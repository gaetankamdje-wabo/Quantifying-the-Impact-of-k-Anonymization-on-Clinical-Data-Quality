# quantify_anonymization_impact()
# ------------------------------------------------------------
# R 4.4 compatible, production-ready implementation that mirrors the
# "Original vs. Anonymized" Assessment
#
#  1. Column presence comparison
#  2. Extended numeric analysis across all common numeric columns
#     (min / max / mean / SD shifts + KS D, sqrt(n_eff) * D, robust p-formatting)
#  3. Extended categorical analysis across all common categorical columns
#     (missing delta, Jaccard overlap ratio, chi-square p over shared levels)
#  4. Ad-hoc single-column comparison with plain-language interpretation and
#     traffic-light verdict / handling recommendation
#
# Corrections vs. prior version (all backward-compatible in behaviour):
#  - Parameter typos fixed: `orgin` → `origin`, `anonimzed` → `anonymized`.
#    Old names retained as deprecated aliases via `...` forwarding.
#  - `cramers_v_tbl`: local variable `c` renamed to `nc` to avoid shadowing
#    base::c().
#  - `sprintf("%d%%", round(...))`: `%d` requires integer type; changed to
#    `sprintf("%.0f%%", ...)` which accepts numeric correctly.
#  - Ad-hoc alignment: raw pre-alignment lengths are now captured before
#    truncation so `N_total_Delta` in the display table is accurate.
#  - `fmt_p` guard: explicit `p == 0L` check added before the
#    `< .Machine$double.xmin` comparison for clarity (both catch underflow,
#    but intent is now unambiguous).
#  - `ks_test_asymp` callers in extended numeric analysis now pass
#    `na.omit()`-cleaned vectors explicitly rather than relying on the
#    internal `is.finite()` guard silently removing NAs.
#  - `safe_stat` guard: `v[is.finite(v)]` already excludes NA/NaN/Inf, so
#    `na.rm = TRUE` is now removed from the inner call to avoid conveying
#    false intent (result is identical; removed for methodological clarity).
#  - Chi-square note: the extended categorical chi-square is intentionally
#    computed over shared levels only (matching app behaviour); categories
#    exclusive to either dataset do not contribute to the test statistic.
#    This is now documented explicitly in the `notes` output slot.
#
# Returns a named list with slots:
#   $column_comparison       data.frame
#   $numeric_analysis        data.frame
#   $categorical_analysis    data.frame
#   $adhoc_column_comparison list | NULL
#   $notes                   list
# ------------------------------------------------------------

quantify_anonymization_impact <- function(
    original_dataframe,
    anonymized_dataframe,
    comparison_column_origin    = NULL,   # corrected spelling
    comparison_column_anonymized = NULL,  # corrected spelling
    # Deprecated aliases kept for backward compatibility
    comparison_column_orgin     = NULL,
    comparison_column_anonimzed = NULL
) {
  
  # ---- Deprecated-alias resolution (backward compatibility) ----------------
  if (is.null(comparison_column_origin)     && !is.null(comparison_column_orgin))
    comparison_column_origin     <- comparison_column_orgin
  if (is.null(comparison_column_anonymized) && !is.null(comparison_column_anonimzed))
    comparison_column_anonymized <- comparison_column_anonimzed
  
  # ==========================================================================
  # HELPERS
  # ==========================================================================
  
  # -- fmt_p -----------------------------------------------------------------
  # Format a p-value to full IEEE 754 double precision (17 significant digits).
  # Explicitly catches exact zero (underflow) as a separate condition for
  # clarity; the numeric guard `p < .Machine$double.xmin` also catches it
  # (since 0 < 2.225e-308 is TRUE), but the explicit check makes intent clear.
  fmt_p <- function(p) {
    if (is.na(p))                                      return("NA")
    if (p == 0 || p < .Machine$double.xmin)
      return(paste0("< ", format(.Machine$double.xmin, scientific = TRUE)))
    format(p, digits = 17, scientific = TRUE, trim = TRUE)
  }
  
  # -- ks_test_asymp ---------------------------------------------------------
  # Two-sample KS test (asymptotic approximation; appropriate for large n and
  # when ties are present).  Returns D, the effective sample size n_eff =
  # n1*n2/(n1+n2), the scaled statistic z = sqrt(n_eff)*D (used to gauge
  # whether the difference is real noise vs. signal), and the asymptotic p.
  #
  # Guard: requires >= 2 distinct finite values in each sample; otherwise all
  # output slots are NA.
  ks_test_asymp <- function(x, y) {
    x <- x[is.finite(x)]
    y <- y[is.finite(y)]
    if (length(unique(x)) < 2L || length(unique(y)) < 2L) {
      return(list(p.value   = NA_real_,
                  statistic = NA_real_,
                  n_eff     = NA_real_,
                  z         = NA_real_))
    }
    res <- tryCatch(stats::ks.test(x, y, exact = FALSE),
                    error = function(e) NULL)
    if (is.null(res)) {
      return(list(p.value   = NA_real_,
                  statistic = NA_real_,
                  n_eff     = NA_real_,
                  z         = NA_real_))
    }
    D     <- unname(res$statistic)
    n_eff <- length(x) * length(y) / (length(x) + length(y))
    z     <- sqrt(n_eff) * D
    list(p.value   = unname(res$p.value),
         statistic = D,
         n_eff     = n_eff,
         z         = z)
  }
  
  # -- cramers_v_tbl ---------------------------------------------------------
  # Cramér's V from a pre-built contingency matrix.  Uses the standard
  # (uncorrected) formula V = sqrt(chi2 / (n * min(r-1, nc-1))).
  # NOTE: local variable renamed from `c` to `nc` to avoid shadowing base::c().
  cramers_v_tbl <- function(tbl) {
    if (length(dim(tbl)) != 2L || any(dim(tbl) < 2L)) {
      return(list(v = NA_real_, p = NA_real_,
                  chi2 = NA_real_, df = NA_real_, n = sum(tbl)))
    }
    chi <- suppressWarnings(
      tryCatch(stats::chisq.test(tbl, correct = FALSE),
               error = function(e) NULL)
    )
    if (is.null(chi)) {
      return(list(v = NA_real_, p = NA_real_,
                  chi2 = NA_real_, df = NA_real_, n = sum(tbl)))
    }
    n  <- sum(tbl)
    r  <- nrow(tbl)
    nc <- ncol(tbl)                      # renamed: was `c`, shadowed base::c()
    v  <- sqrt(as.numeric(chi$statistic) / (n * min(r - 1L, nc - 1L)))
    list(v    = as.numeric(v),
         p    = as.numeric(chi$p.value),
         chi2 = as.numeric(chi$statistic),
         df   = as.numeric(chi$parameter),
         n    = n)
  }
  
  # -- force_num -------------------------------------------------------------
  force_num <- function(v) {
    if (is.numeric(v) || is.integer(v)) return(as.numeric(v))
    suppressWarnings(as.numeric(v))
  }
  
  # -- safe_stat -------------------------------------------------------------
  # Applies `fun` to the finite elements of `v`.
  # `na.rm` is intentionally NOT passed to `fun` because `is.finite()` has
  # already excluded NA/NaN/Inf; passing na.rm = TRUE would be misleading.
  safe_stat <- function(fun, v) {
    v <- v[is.finite(v)]
    if (length(v) == 0L) return(NA_real_)
    fun(v)
  }
  
  # ==========================================================================
  # VALIDATION
  # ==========================================================================
  if (is.null(original_dataframe) || is.null(anonymized_dataframe))
    stop("Both original_dataframe and anonymized_dataframe must be provided.")
  
  df_orig <- as.data.frame(original_dataframe,  stringsAsFactors = FALSE)
  df_anon <- as.data.frame(anonymized_dataframe, stringsAsFactors = FALSE)
  
  # ==========================================================================
  # 1) COLUMN PRESENCE
  # ==========================================================================
  orig_cols <- colnames(df_orig)
  anon_cols <- colnames(df_anon)
  all_cols  <- union(orig_cols, anon_cols)
  
  column_comparison <- data.frame(
    Column        = all_cols,
    In_Original   = ifelse(all_cols %in% orig_cols, "Yes", "No"),
    In_Anonymized = ifelse(all_cols %in% anon_cols, "Yes", "No"),
    stringsAsFactors = FALSE
  )
  
  # ==========================================================================
  # 2) EXTENDED NUMERIC ANALYSIS
  # ==========================================================================
  common_cols  <- intersect(orig_cols, anon_cols)
  numeric_cols <- common_cols[
    vapply(df_orig[, common_cols, drop = FALSE],
           function(x) is.numeric(x) || is.integer(x),
           logical(1L))
  ]
  
  if (length(numeric_cols) > 0L) {
    res_list <- lapply(numeric_cols, function(col) {
      orig <- force_num(df_orig[[col]])
      anon <- force_num(df_anon[[col]])
      
      # Pass finite-only vectors explicitly; ks_test_asymp would guard
      # internally, but explicit na.omit makes the caller's intent clear.
      ks_out <- ks_test_asymp(orig[is.finite(orig)], anon[is.finite(anon)])
      
      data.frame(
        Column         = col,
        Orig_N_total   = length(orig),
        Anon_N_total   = length(anon),
        N_total_Delta  = length(anon) - length(orig),
        
        Orig_Distinct  = length(unique(orig[is.finite(orig)])),
        Anon_Distinct  = length(unique(anon[is.finite(anon)])),
        Distinct_Delta = length(unique(anon[is.finite(anon)])) -
          length(unique(orig[is.finite(orig)])),
        
        Orig_Min       = safe_stat(min,  orig),
        Anon_Min       = safe_stat(min,  anon),
        Min_Shift      = safe_stat(min,  anon) - safe_stat(min,  orig),
        
        Orig_Max       = safe_stat(max,  orig),
        Anon_Max       = safe_stat(max,  anon),
        Max_Shift      = safe_stat(max,  anon) - safe_stat(max,  orig),
        
        Orig_Mean      = safe_stat(mean, orig),
        Anon_Mean      = safe_stat(mean, anon),
        Mean_Shift     = safe_stat(mean, anon) - safe_stat(mean, orig),
        
        Orig_SD        = safe_stat(stats::sd, orig),
        Anon_SD        = safe_stat(stats::sd, anon),
        SD_Shift       = safe_stat(stats::sd, anon) - safe_stat(stats::sd, orig),
        
        KS_D           = ks_out$statistic,
        KS_Z_sqrtNeffD = ks_out$z,
        KS_p_value     = fmt_p(ks_out$p.value),
        stringsAsFactors = FALSE
      )
    })
    numeric_analysis <- do.call(rbind, res_list)
    rownames(numeric_analysis) <- NULL
  } else {
    numeric_analysis <- data.frame(
      Message = "No common numeric columns.",
      stringsAsFactors = FALSE
    )
  }
  
  # ==========================================================================
  # 3) EXTENDED CATEGORICAL ANALYSIS
  # ==========================================================================
  cat_cols <- common_cols[
    vapply(df_orig[, common_cols, drop = FALSE],
           function(x) is.character(x) || is.factor(x),
           logical(1L))
  ]
  
  if (length(cat_cols) > 0L) {
    res_list <- lapply(cat_cols, function(col) {
      orig <- as.character(df_orig[[col]])
      anon <- as.character(df_anon[[col]])
      
      orig_missing <- sum(is.na(orig) | !nzchar(orig, keepNA = FALSE))
      anon_missing <- sum(is.na(anon) | !nzchar(anon, keepNA = FALSE))
      
      orig_levels <- unique(orig[!is.na(orig) & nzchar(orig)])
      anon_levels <- unique(anon[!is.na(anon) & nzchar(anon)])
      
      # Jaccard overlap ratio: |intersect| / |union|
      denom_u      <- length(union(orig_levels, anon_levels))
      overlap_ratio <- if (denom_u == 0L) NA_real_
      else length(intersect(orig_levels, anon_levels)) / denom_u
      
      tab_orig <- table(orig)
      tab_anon <- table(anon)
      
      # Chi-square over shared levels only (mirrors app behaviour).
      # Categories exclusive to one dataset are not tested here; they are
      # captured by `Distinct_Delta` and `Overlap_Ratio` instead.
      common_levels <- intersect(names(tab_orig), names(tab_anon))
      chisq_p <- if (length(common_levels) > 1L) {
        tryCatch(
          stats::chisq.test(
            rbind(tab_orig[common_levels],
                  tab_anon[common_levels]),
            correct = FALSE
          )$p.value,
          error = function(e) NA_real_
        )
      } else {
        NA_real_
      }
      
      data.frame(
        Column         = col,
        Orig_N_total   = length(orig),
        Anon_N_total   = length(anon),
        N_total_Delta  = length(anon) - length(orig),
        Orig_Distinct  = length(orig_levels),
        Anon_Distinct  = length(anon_levels),
        Distinct_Delta = length(anon_levels) - length(orig_levels),
        Orig_Missing   = orig_missing,
        Anon_Missing   = anon_missing,
        Missing_Delta  = anon_missing - orig_missing,
        Overlap_Ratio  = overlap_ratio,
        ChiSq_p_value  = chisq_p,
        stringsAsFactors = FALSE
      )
    })
    categorical_analysis <- do.call(rbind, res_list)
    rownames(categorical_analysis) <- NULL
  } else {
    categorical_analysis <- data.frame(
      Message = "No common categorical columns.",
      stringsAsFactors = FALSE
    )
  }
  
  # ==========================================================================
  # 4) AD-HOC SINGLE-COLUMN COMPARISON + HANDLING RECOMMENDATION
  # ==========================================================================
  adhoc <- NULL
  
  if (!is.null(comparison_column_origin) && !is.null(comparison_column_anonymized)) {
    
    orig_vec <- comparison_column_origin
    anon_vec <- comparison_column_anonymized
    
    # Capture raw lengths BEFORE alignment so the display table can show the
    # true row-count delta (previously these were lost).
    n_orig_raw <- length(orig_vec)
    n_anon_raw <- length(anon_vec)
    n_min      <- min(n_orig_raw, n_anon_raw)
    
    if (n_orig_raw != n_anon_raw) {
      orig_vec <- orig_vec[seq_len(n_min)]
      anon_vec <- anon_vec[seq_len(n_min)]
    }
    
    is_num_o <- is.numeric(orig_vec) || is.integer(orig_vec)
    is_num_a <- is.numeric(anon_vec) || is.integer(anon_vec)
    
    # ------------------------------------------------------------------
    # 4a) NUMERIC ad-hoc
    # ------------------------------------------------------------------
    if (is_num_o && is_num_a) {
      
      x <- as.numeric(orig_vec)
      y <- as.numeric(anon_vec)
      
      ks_out <- ks_test_asymp(stats::na.omit(x), stats::na.omit(y))
      
      n_dist_o <- length(unique(stats::na.omit(x)))
      n_dist_a <- length(unique(stats::na.omit(y)))
      
      # Mixed-type display table: numeric values are coerced to character
      # because the last rows hold formatted strings.  This is intentional
      # (display-only output); numeric precision is preserved in `meta`.
      adhoc_table <- data.frame(
        Metric = c(
          "N total (pre-alignment)", "N total (compared)",
          "Distinct count",
          "Min", "Max", "Mean", "SD",
          "KS D (effect size)", "sqrt(n_eff)*D",
          "KS p-value (asymptotic; ties present)"
        ),
        Original = as.character(c(
          n_orig_raw, n_min,
          n_dist_o,
          min(x, na.rm = TRUE), max(x, na.rm = TRUE),
          mean(x, na.rm = TRUE), stats::sd(x, na.rm = TRUE),
          ks_out$statistic, ks_out$z,
          NA  # placeholder; overwritten below
        )),
        Anonymized = as.character(c(
          n_anon_raw, n_min,
          n_dist_a,
          min(y, na.rm = TRUE), max(y, na.rm = TRUE),
          mean(y, na.rm = TRUE), stats::sd(y, na.rm = TRUE),
          "", "", ""
        )),
        Delta = as.character(c(
          n_anon_raw - n_orig_raw, 0L,
          n_dist_a - n_dist_o,
          min(y, na.rm = TRUE) - min(x, na.rm = TRUE),
          max(y, na.rm = TRUE) - max(x, na.rm = TRUE),
          mean(y, na.rm = TRUE) - mean(x, na.rm = TRUE),
          stats::sd(y, na.rm = TRUE) - stats::sd(x, na.rm = TRUE),
          "", "", ""
        )),
        check.names    = FALSE,
        stringsAsFactors = FALSE
      )
      adhoc_table$Original[10L] <- fmt_p(ks_out$p.value)
      
      meta <- list(
        type      = "numeric",
        mean_orig = mean(x, na.rm = TRUE),
        mean_anon = mean(y, na.rm = TRUE),
        sd_orig   = stats::sd(x, na.rm = TRUE),
        sd_anon   = stats::sd(y, na.rm = TRUE),
        ks_p_raw  = ks_out$p.value,
        ks_D      = ks_out$statistic,
        ks_z      = ks_out$z,
        ks_n_eff  = ks_out$n_eff
      )
      
      # ---- Traffic-light verdict (thresholds mirror the app) ---------------
      D          <- suppressWarnings(as.numeric(meta$ks_D))
      z          <- suppressWarnings(as.numeric(meta$ks_z))
      mean_shift <- suppressWarnings(as.numeric(meta$mean_anon - meta$mean_orig))
      sd_shift   <- suppressWarnings(as.numeric(meta$sd_anon  - meta$sd_orig))
      sd_orig_v  <- suppressWarnings(as.numeric(meta$sd_orig))
      
      rel_mean   <- if (isTRUE(sd_orig_v > 0)) abs(mean_shift) / sd_orig_v
      else NA_real_
      mean_band  <- if (is.na(rel_mean))       "unknown"
      else if (rel_mean <= 0.10)  "tiny"
      else if (rel_mean <= 0.30)  "moderate"
      else                        "big"
      
      spread_pct  <- if (isTRUE(sd_orig_v != 0)) 100 * (sd_shift / sd_orig_v)
      else NA_real_
      spread_band <- if (is.na(spread_pct))        "unknown"
      else if (abs(spread_pct) <= 10) "about the same"
      else if (spread_pct > 10)       "more spread out"
      else                            "more squeezed"
      
      shape_band  <- if (is.na(D))      "unknown"
      else if (D < 0.05) "small"
      else if (D < 0.10) "clear"
      else               "large"
      
      verdict <- if (!is.na(D) && D < 0.05 &&
                     (is.na(rel_mean)   || rel_mean   <= 0.10) &&
                     (is.na(spread_pct) || abs(spread_pct) <= 10)) {
        "\u2705 No meaningful change"
      } else if (!is.na(D) && D < 0.10 &&
                 (is.na(rel_mean) || rel_mean <= 0.30)) {
        "\u26a0\ufe0f Noticeable change"
      } else {
        "\u274c Big change"
      }
      
      real_change_msg <- if (is.na(z))          "n/a"
      else if (z < 1.36)      "Probably just random noise"
      else if (z < 1.63)      "Likely a real difference"
      else                    "Almost certainly a real difference"
      
      handling_recommendation <- list(
        verdict                 = verdict,
        typical_value_shift     = mean_shift,
        typical_value_band      = mean_band,
        spread_shift            = sd_shift,
        spread_band             = spread_band,
        distribution_shape_band = shape_band,
        ks = list(
          D           = D,
          z           = z,
          n_eff       = meta$ks_n_eff,
          p           = meta$ks_p_raw,
          p_formatted = fmt_p(meta$ks_p_raw)
        ),
        decision_support = list(
          "Green (\u2705)"  = paste0("Safe to use as-is for most analyses; ",
                                     "verify critical endpoints if the study is ",
                                     "sensitive to small distributional shifts."),
          "Yellow (\u26a0\ufe0f)" = paste0("Expect some impact; document the shift, ",
                                           "run sensitivity analyses, and consider mitigation ",
                                           "(e.g., calibration, stratification checks)."),
          "Red (\u274c)"    = paste0("High risk of analytic drift; avoid direct ",
                                     "comparability, redesign anonymization, or add ",
                                     "reconstruction/mitigation steps before release.")
        ),
        real_change_or_noise = real_change_msg
      )
      
      adhoc <- list(
        type                    = "numeric",
        comparison_table        = adhoc_table,
        meta                    = meta,
        handling_recommendation = handling_recommendation
      )
      
      # ------------------------------------------------------------------
      # 4b) CATEGORICAL / mixed ad-hoc
      # ------------------------------------------------------------------
    } else {
      
      orig_chr <- as.character(orig_vec)
      anon_chr <- as.character(anon_vec)
      
      orig_missing <- sum(is.na(orig_chr) | !nzchar(orig_chr, keepNA = FALSE))
      anon_missing <- sum(is.na(anon_chr) | !nzchar(anon_chr, keepNA = FALSE))
      
      lv_orig  <- unique(orig_chr[!is.na(orig_chr) & nzchar(orig_chr)])
      lv_anon  <- unique(anon_chr[!is.na(anon_chr) & nzchar(anon_chr)])
      overlap  <- length(intersect(lv_orig, lv_anon)) /
        max(1L, length(union(lv_orig, lv_anon)))
      
      tab_o    <- table(orig_chr, useNA = "no")
      tab_a    <- table(anon_chr, useNA = "no")
      levels_u <- union(names(tab_o), names(tab_a))
      
      # Extend both frequency vectors to the full union of level names,
      # filling absent levels with 0.
      tab_o_u              <- tab_o[levels_u]
      tab_o_u[is.na(tab_o_u)] <- 0L
      tab_a_u              <- tab_a[levels_u]
      tab_a_u[is.na(tab_a_u)] <- 0L
      
      tbl          <- rbind(Original   = as.numeric(tab_o_u),
                            Anonymized = as.numeric(tab_a_u))
      colnames(tbl) <- levels_u
      
      vstats <- cramers_v_tbl(tbl)
      chi_p  <- vstats$p
      v_val  <- vstats$v
      
      extra_rows <- data.frame(
        Statistic = c("N total pre-alignment (orig)",
                      "N total pre-alignment (anon)",
                      "N compared (aligned)",
                      "Distinct levels (orig)",
                      "Distinct levels (anon)"),
        Value     = as.character(c(n_orig_raw, n_anon_raw, n_min,
                                   length(lv_orig), length(lv_anon))),
        check.names      = FALSE,
        stringsAsFactors = FALSE
      )
      
      base_rows <- data.frame(
        Statistic = c(
          "Orig missing", "Anon missing", "Missing \u0394",
          "Level overlap (Jaccard)", "Cram\u00e9r's V (effect size)",
          "Chi-Sq p-value (marginal distribution, shared levels only)"
        ),
        Value = as.character(c(
          orig_missing, anon_missing, anon_missing - orig_missing,
          overlap, v_val, fmt_p(chi_p)
        )),
        check.names      = FALSE,
        stringsAsFactors = FALSE
      )
      
      meta <- list(
        type          = "categorical",
        missing_delta = anon_missing - orig_missing,
        overlap       = as.numeric(overlap),
        chi_p         = as.numeric(chi_p),
        cramer_v      = as.numeric(v_val)
      )
      
      # ---- Traffic-light verdict (thresholds mirror the app) ---------------
      v        <- suppressWarnings(as.numeric(meta$cramer_v))
      p        <- suppressWarnings(as.numeric(meta$chi_p))
      overlapN <- suppressWarnings(as.numeric(meta$overlap))
      miss_d   <- suppressWarnings(as.numeric(meta$missing_delta))
      
      denom    <- length(orig_chr)   # aligned length
      miss_pct <- if (!is.na(denom) && denom > 0L && !is.na(miss_d))
        100 * miss_d / denom
      else NA_real_
      
      mix_band    <- if (is.na(v))       "unknown"
      else if (v < 0.10)  "small"
      else if (v < 0.30)  "medium"
      else                "large"
      
      # FIX: was `sprintf("%d%%", round(...))` — `%d` requires integer;
      # use `%.0f` which correctly accepts numeric.
      overlap_txt <- if (is.na(overlapN)) "n/a"
      else sprintf("%.0f%%", 100 * overlapN)
      
      verdict <- if ((is.na(v) || v < 0.10) &&
                     (is.na(overlapN) || overlapN >= 0.80) &&
                     (is.na(miss_pct) || miss_pct <= 5)) {
        "\u2705 No meaningful change"
      } else if ((is.na(v) || v < 0.30) &&
                 (is.na(overlapN) || overlapN >= 0.60) &&
                 (is.na(miss_pct) || miss_pct <= 10)) {
        "\u26a0\ufe0f Noticeable change"
      } else {
        "\u274c Big change"
      }
      
      real_change_msg <- if (is.na(p))       "n/a"
      else if (p < 0.01)  "Almost certainly a real difference"
      else if (p < 0.05)  "Likely a real difference"
      else                "Probably just random noise"
      
      handling_recommendation <- list(
        verdict            = verdict,
        category_mix_band  = mix_band,
        shared_labels_kept = overlap_txt,
        extra_blanks_after_anonymization = list(
          delta               = miss_d,
          approx_percent_points = miss_pct
        ),
        chi_square = list(
          p           = p,
          p_formatted = fmt_p(p)
        ),
        cramer_v = v,
        decision_support = list(
          "Green (\u2705)"  = paste0("Category distribution is stable; downstream ",
                                     "categorical analyses should remain comparable."),
          "Yellow (\u26a0\ufe0f)" = paste0("Expect some drift; document changes, check ",
                                           "key strata, and run sensitivity analyses."),
          "Red (\u274c)"    = paste0("Large distortion; treat as non-comparable or apply ",
                                     "mitigation (e.g., harmonize categories, revise ",
                                     "generalization/suppression rules).")
        ),
        real_change_or_noise = real_change_msg
      )
      
      adhoc <- list(
        type                    = "categorical",
        comparison_table        = rbind(extra_rows, base_rows),
        meta                    = meta,
        handling_recommendation = handling_recommendation
      )
    }
  }
  
  # ==========================================================================
  # RETURN BUNDLE
  # ==========================================================================
  list(
    column_comparison       = column_comparison,
    numeric_analysis        = numeric_analysis,
    categorical_analysis    = categorical_analysis,
    adhoc_column_comparison = adhoc,
    notes = list(
      numeric = paste0(
        "Numeric comparison: min/max/mean/SD shifts + KS D (effect size), ",
        "sqrt(n_eff)*D (signal vs. noise gauge), and full-precision asymptotic p. ",
        "Distinct counts and N-delta use finite values only."
      ),
      categorical = paste0(
        "Categorical comparison: missing delta, Jaccard overlap ratio ",
        "(|intersect|/|union|), and chi-square p. ",
        "IMPORTANT: the chi-square is computed over shared level names only; ",
        "categories exclusive to one dataset are NOT included in the test ",
        "statistic but are reflected in Distinct_Delta and Overlap_Ratio. ",
        "Ad-hoc additionally provides Cram\u00e9r's V over the full union of levels."
      ),
      alignment = paste0(
        "When comparison vectors differ in length, both are trimmed to the ",
        "shorter length before statistics are computed. Pre-alignment lengths ",
        "are reported in the ad-hoc display table."
      )
    )
  )
}

# ==========================================================================
# Usage example
# ==========================================================================
# result <- quantify_anonymization_impact(
#   original_dataframe       = demo_original_data,
#   anonymized_dataframe     = demo_anonymized_data,
#   comparison_column_origin    = demo_original_data$icd,
#   comparison_column_anonymized = demo_anonymized_data$icd
# )
#
# result$column_comparison
# result$numeric_analysis
# result$categorical_analysis
# result$adhoc_column_comparison$handling_recommendation$verdict