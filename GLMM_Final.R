# ============================================================================
# GLMM_v3.R
# ----------------------------------------------------------------------------
# TRACE Luquillo soil CO2 warming experiment.
# Final consolidated v3 pipeline. Supersedes all prior v2 / v3 scripts.
#
# DATA
#   T_VWC.csv      hourly soil temperature & VWC (temp/VWC primary)
#   FLUX.csv       hourly soil CO2 flux         (flux primary + mechanistic)
#   CLIMATE.csv    monthly cwd, prec, binary Dry/Wet season
#   Date filter:   >= 2018-09-01 applied to both pipelines.
#
# MODELS  (all share the random structure; Gaussian REML; nlminb optimiser,
#          with dispersed-start retries if the Hessian is not certified)
#
#   Temperature   temp_mean ~ Treatment * Year_f + season
#                          + (1|Plot) + ar1(time_ou + 0 | Plot_Year)
#                 dispformula = ~ 1
#
#   VWC           vwc_pct ~ Treatment * Year_f + season
#                        + (1|Plot) + ar1(time_ou + 0 | Plot_Year)
#                 dispformula = ~ 1; vwc_pct = vwc_mean * 100
#
#   Flux          log_flux ~ Treatment * Year_f + Treatment * season
#                         + (1|Plot) + ar1(time_ou + 0 | Plot_Year)
#                 dispformula = ~ 1
#
#   Mechanistic   log_flux ~ Treatment + Year_f + season + temp_c + vwc_c
#                         + I(vwc_c^2) + Treatment:(temp_c + vwc_c + I(vwc_c^2))
#                         + (1|Plot) + ar1(time_ou + 0 | Plot_Year)
#                 dispformula = ~ 1
#
# DESIGN RATIONALE  (empirically validated against raw per-year differences)
#
#   AR(1) within Plot_Year. ar1(time_ou + 0 | Plot) spanning the full record
#   absorbs year-level fixed-effect variation, distorting per-year contrasts.
#   Restricting the AR(1) to within Plot_Year decouples year-level signal
#   from the autocorrelation process; per-year model contrasts then track
#   raw year-by-year warmed-control differences to within 0.5 unit.
#
#   No Treatment:season for temperature and VWC. The Hessian goes non-PD when
#   Treatment:season is included, because at this sample size the warming-
#   by-season interaction is not robustly identifiable for these two
#   responses. Flux retains Treatment:season as the primary scientific
#   interest is the treatment effect on flux.
#
#   Mechanistic dispformula = ~ 1. The richer ~ Treatment + season is non-
#   identifiable here because temperature and VWC absorb the seasonal and
#   treatment mean structure that the dispformula otherwise models in the
#   variance.
#
#   Default nlminb is used for all fits. If it returns a finite fit with a
#   non-positive-definite Hessian, the script retries from dispersed random-
#   effect starts; all reported coefficient extractors assume the documented
#   default treatment coding.
#
# CONTRASTS
#   adjust = "none" (unadjusted p-values; comparison count stated in methods)
#   Marginal warmed - control                       (1 per response)
#   Per-Year warmed - control                       (8 per response)
#   Per-season warmed - control                     (flux only: 2)
#   Per-cell sample sizes attached; empty cells flagged.
#
# DERIVED QUANTITIES FROM MECHANISTIC MODEL
#   Q10 by treatment      from Treatment:temp_mean slopes, Q10 = exp(10*slope)
#   VWC optimum by trt    from quadratic peak x* = -b1 / (2*b2)
#   95% CIs via delta method on the relevant coefficient vector.
#
# OUTPUTS
#   OUTPUT_GLMM_V3/tables/         CSVs (coefficients, contrasts, coverage)
#   OUTPUT_GLMM_V3/objects/        RDS with all fits and contrasts
#   OUTPUT_GLMM_V3/REPORT.txt      comprehensive plain-text report
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(readr)
  library(glmmTMB)
  library(emmeans)
  library(purrr)
  library(tibble)
  library(stringr)
})

set.seed(42)

# ---------------------------- configuration ---------------------------------
# Run from the repository/project directory; do not override the caller's working
# directory with a machine-specific path.
tv_path     <- "T_VWC.csv"
flux_path   <- "FLUX.csv"
clim_path   <- "CLIMATE.csv"
date_filter <- as.Date("2018-09-01")
out_dir     <- "OUTPUT_GLMM_FINAL"
report_path <- file.path(out_dir, "REPORT.txt")

MIN_PLOT_YEAR_DAYS    <- 30L   # coverage warning threshold
EMPTY_CELL_THRESHOLD  <- 0L    # cells with min(n_ctrl, n_warm) <= threshold flagged
MAX_RETRY             <- 10L   # dispersed-start retries to certify pdHess (vwc_inter needs 1)

stopifnot(file.exists(tv_path), file.exists(flux_path), file.exists(clim_path))
dir.create(out_dir,                          recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "tables"),     showWarnings = FALSE)
dir.create(file.path(out_dir, "objects"),    showWarnings = FALSE)

# Optimiser controls for the default nlminb fit.
ctrl_default <- glmmTMBControl(optCtrl = list(iter.max = 2000, eval.max = 2000))

# ----------------------------- helpers --------------------------------------
`%||%` <- function(x, y) if (is.null(x)) y else x

is_failed_glmm <- function(x) {
  if (inherits(x, "error") || inherits(x, "try-error")) return(TRUE)
  if (!inherits(x, "glmmTMB")) return(TRUE)
  ll <- tryCatch(as.numeric(logLik(x)), error = function(e) NA_real_)
  !is.finite(ll)
}
is_pdHess <- function(x) {
  if (is_failed_glmm(x)) return(FALSE)
  tryCatch(isTRUE(x$sdr$pdHess), error = function(e) FALSE)
}
safe_AIC <- function(model) {
  if (is_failed_glmm(model)) return(NA_real_)
  tryCatch(as.numeric(AIC(model)), error = function(e) NA_real_)
}

# Extract AR(1) phi, tau (decorrelation in observation steps), Plot sigma.
extract_random_summary <- function(model) {
  out <- list(phi = NA_real_, tau = NA_real_,
              sigma_plot = NA_real_, sigma_plot_year = NA_real_)
  if (is_failed_glmm(model)) return(out)
  vc <- tryCatch(VarCorr(model)$cond, error = function(e) NULL)
  if (is.null(vc) || length(vc) == 0) return(out)
  for (i in seq_along(vc)) {
    nm     <- names(vc)[i]
    sd_val <- as.numeric(attr(vc[[i]], "stddev"))[1]
    cm     <- attr(vc[[i]], "correlation")
    if (!is.null(cm) && is.matrix(cm) && nrow(cm) >= 2) {
      out$phi             <- as.numeric(cm[1, 2])
      out$sigma_plot_year <- sd_val
      if (!is.na(out$phi) && out$phi > 0 && out$phi < 1)
        out$tau <- -1 / log(out$phi)
    } else {
      out$sigma_plot <- sd_val
    }
  }
  out
}

# Fit with dispersed-start retry. The retry is seeded from a finite, but
# uncertified, default fit and perturbs only the random-effect parameters.
fit_with_rescue <- function(formula, dispformula, data, label) {
  fit_std <- tryCatch(glmmTMB(
    formula = formula, dispformula = dispformula,
    data = data, family = gaussian(), REML = TRUE, control = ctrl_default
  ), error = function(e) e)

  if (!is_failed_glmm(fit_std) && is_pdHess(fit_std)) {
    return(list(fit = fit_std, tactic = "standard"))
  }
  seed_ok <- inherits(fit_std, "glmmTMB") &&
    !is.null(tryCatch(getME(fit_std, "theta"), error = function(e) NULL))
  if (!seed_ok) {
    return(list(fit = NULL, tactic = "FAILED (no usable fit to seed retry)"))
  }

  th0 <- getME(fit_std, "theta")
  be0 <- fixef(fit_std)$cond
  bd0 <- fixef(fit_std)$disp
  set.seed(101)
  for (i in seq_len(MAX_RETRY)) {
    st <- list(beta = be0, betadisp = bd0,
               theta = th0 + rnorm(length(th0), 0, 0.5))
    fit_r <- tryCatch(glmmTMB(
      formula = formula, dispformula = dispformula,
      data = data, family = gaussian(), REML = TRUE,
      start = st, control = ctrl_default
    ), error = function(e) e)
    if (!is_failed_glmm(fit_r) && is_pdHess(fit_r)) {
      return(list(fit = fit_r, tactic = sprintf("dispersed start %d", i)))
    }
  }
  list(fit = fit_std,
       tactic = sprintf("standard (pdHess caveat; %d retries exhausted)", MAX_RETRY))
}

# Per-cell sample sizes (n_control, n_warmed, n_min) used for empty-cell
# flagging in per-(Year x season) contrast tables.
build_cell_n <- function(data, by_vars = c("Year_f", "season")) {
  # Complete treatment levels before widening so a missing arm is reported as
  # zero rather than causing n_control/n_warmed to be absent.
  data %>%
    count(across(all_of(by_vars)), Treatment, name = "n_obs") %>%
    tidyr::complete(!!!rlang::syms(by_vars), Treatment,
                    fill = list(n_obs = 0L)) %>%
    pivot_wider(names_from = Treatment, values_from = n_obs,
                values_fill = 0L, names_prefix = "n_") %>%
    mutate(n_min = pmin(n_control, n_warmed))
}

# Unadjusted-p contrast extraction. If cell_n + by-keys provided, joins
# per-cell sample sizes and adds an empty_cell flag.
extract_contrast <- function(model, specs, cell_n = NULL, by_keys = NULL) {
  if (is_failed_glmm(model)) return(tibble(status = "model failed"))
  ct <- tryCatch(
    as_tibble(summary(
      contrast(emmeans(model, specs = specs), method = "revpairwise",
               adjust = "none"),
      infer = c(TRUE, TRUE))),
    error = function(e) tibble(status = paste0("emmeans failed: ", e$message))
  )
  if ("status" %in% names(ct)) return(ct)
  required_cols <- c("contrast", "estimate", "SE", "lower.CL", "upper.CL", "p.value")
  missing_cols <- setdiff(required_cols, names(ct))
  if (length(missing_cols) > 0L) {
    return(tibble(status = paste0("emmeans output missing columns: ",
                                  paste(missing_cols, collapse = ", "))))
  }
  if (xor(is.null(cell_n), is.null(by_keys))) {
    return(tibble(status = "cell_n and by_keys must be supplied together"))
  }
  if (!is.null(cell_n)) {
    missing_ct_keys <- setdiff(by_keys, names(ct))
    missing_n_keys <- setdiff(by_keys, names(cell_n))
    if (length(missing_ct_keys) > 0L || length(missing_n_keys) > 0L ||
        !"n_min" %in% names(cell_n)) {
      return(tibble(status = "contrast coverage join is missing required keys"))
    }
    ct <- left_join(ct, cell_n, by = by_keys)
    if (anyNA(ct$n_min)) {
      return(tibble(status = "contrast coverage join did not match every contrast"))
    }
    ct <- ct %>%
      mutate(empty_cell = n_min <= EMPTY_CELL_THRESHOLD,
             notes = if_else(empty_cell, "[empty/low-info cell]", ""))
  }
  ct
}

# Q10 by treatment under the CENTRED-INTERACTION mechanistic parameterisation.
# control slope = temp_c ; warmed slope = temp_c + Treatmentwarmed:temp_c.
# Slopes are centring-invariant, so Q10 = exp(10*slope) is on the real scale.
derive_q10 <- function(model) {
  if (is_failed_glmm(model))
    return(tibble(treatment = c("control","warmed"),
                  slope = NA, SE = NA, Q10 = NA, Q10_lo = NA, Q10_hi = NA,
                  status = "model failed"))
  fe <- fixef(model)$cond; ve <- vcov(model)$cond
  nm_m <- "temp_c"; nm_i <- "Treatmentwarmed:temp_c"
  if (!nm_m %in% names(fe))
    return(tibble(treatment = c("control","warmed"),
                  slope = NA, SE = NA, Q10 = NA, Q10_lo = NA, Q10_hi = NA,
                  status = "temp_c not found"))
  if (!nm_i %in% names(fe))
    return(tibble(treatment = c("control", "warmed"),
                  slope = NA, SE = NA, Q10 = NA, Q10_lo = NA, Q10_hi = NA,
                  status = "treatment-specific temperature slope not found"))
  b_m <- fe[nm_m]; b_i <- fe[nm_i]
  slope_c <- as.numeric(b_m); slope_w <- as.numeric(b_m + b_i)
  s_c <- sqrt(ve[nm_m, nm_m])
  s_w <- sqrt(max(ve[nm_m,nm_m] + ve[nm_i,nm_i] + 2*ve[nm_m,nm_i], 0))
  tibble(
    treatment = c("control", "warmed"),
    slope     = c(slope_c, slope_w),
    SE        = c(as.numeric(s_c), as.numeric(s_w)),
    Q10       = c(exp(10*slope_c), exp(10*slope_w)),
    Q10_lo    = c(exp(10*(slope_c - 1.96*s_c)), exp(10*(slope_w - 1.96*s_w))),
    Q10_hi    = c(exp(10*(slope_c + 1.96*s_c)), exp(10*(slope_w + 1.96*s_w)))
  )
}

# VWC optimum by treatment under the CENTRED-INTERACTION parameterisation.
# control: b1=vwc_c, b2=I(vwc_c^2); warmed: b1=vwc_c + Tw:vwc_c,
# b2=I(vwc_c^2)+Tw:I(vwc_c^2). Vertex x* = -b1/(2 b2) is on the CENTRED scale;
# add VWC_CENTRE to return the optimum on the real vwc_mean (fraction) scale.
# 95% CI via delta method: control over (b1,b2); warmed over all four coefs.
# VWC_CENTRE is passed in from the fitting frame (grand mean of vwc_mean).
derive_vwc_optimum <- function(model, VWC_CENTRE) {
  if (is_failed_glmm(model))
    return(tibble(treatment = c("control","warmed"),
                  optimum = NA, SE = NA, opt_lo = NA, opt_hi = NA,
                  note = "model failed"))
  fe <- fixef(model)$cond; ve <- vcov(model)$cond
  b1m <- "vwc_c"; b2m <- "I(vwc_c^2)"
  b1i <- "Treatmentwarmed:vwc_c"; b2i <- "Treatmentwarmed:I(vwc_c^2)"
  if (!all(c(b1m, b2m) %in% names(fe)))
    return(tibble(treatment = c("control","warmed"),
                  optimum = NA, SE = NA, opt_lo = NA, opt_hi = NA,
                  note = "main quadratic coefficients missing"))
  has_i <- all(c(b1i, b2i) %in% names(fe))
  B1m <- fe[b1m]; B2m <- fe[b2m]

  ctrl <- local({
    if (B2m >= 0) return(list(opt = NA, se = NA,
                              note = sprintf("non-concave control (b2 = %.3f >= 0)", B2m)))
    v <- as.numeric(-B1m / (2 * B2m))
    g <- c(-1/(2*B2m), B1m/(2*B2m^2))
    V <- ve[c(b1m, b2m), c(b1m, b2m)]
    list(opt = v + VWC_CENTRE, se = sqrt(max(as.numeric(t(g) %*% V %*% g), 0)), note = "")
  })
  warm <- local({
    if (!has_i) return(list(opt = NA, se = NA, note = "warmed interaction terms missing"))
    B1i <- fe[b1i]; B2i <- fe[b2i]; B1 <- B1m + B1i; B2 <- B2m + B2i
    if (B2 >= 0) return(list(opt = NA, se = NA,
                             note = sprintf("non-concave warmed (b2 = %.3f >= 0)", B2)))
    v <- as.numeric(-B1 / (2 * B2))
    g <- c(-1/(2*B2), -1/(2*B2), B1/(2*B2^2), B1/(2*B2^2))  # order b1m,b1i,b2m,b2i
    V <- ve[c(b1m, b1i, b2m, b2i), c(b1m, b1i, b2m, b2i)]
    list(opt = v + VWC_CENTRE, se = sqrt(max(as.numeric(t(g) %*% V %*% g), 0)), note = "")
  })
  tibble(
    treatment = c("control", "warmed"),
    optimum   = c(ctrl$opt, warm$opt),
    SE        = c(ctrl$se,  warm$se),
    opt_lo    = c(ctrl$opt - 1.96 * ctrl$se, warm$opt - 1.96 * warm$se),
    opt_hi    = c(ctrl$opt + 1.96 * ctrl$se, warm$opt + 1.96 * warm$se),
    note      = c(ctrl$note, warm$note)
  )
}

# Formatting for the plain-text report.
fmt   <- function(x, d = 3) if (is.na(x)) "NA" else
                            formatC(x, format = "f", digits = d, big.mark = ",")
fmt_p <- function(p) if (is.na(p)) "NA" else
                     if (p < 0.0001) "<0.0001" else
                     formatC(p, format = "f", digits = 4)

contrast_to_text <- function(ct, group_col = NULL, value_label = "estimate",
                              show_n = FALSE) {
  if (is.null(ct) || nrow(ct) == 0) return("  [no contrasts]\n")
  if ("status" %in% names(ct))
    return(paste0("  [", paste(ct$status, collapse = "; "), "]\n"))
  has_group <- !is.null(group_col) && group_col %in% names(ct)
  notes <- if ("notes" %in% names(ct)) ct$notes else rep("", nrow(ct))

  if (has_group && show_n) {
    header <- sprintf("  %-12s %10s %10s %22s %10s  %6s %6s  %s",
                      group_col, value_label, "SE", "95% CI", "p",
                      "n_ctrl", "n_warm", "notes")
  } else if (has_group) {
    header <- sprintf("  %-12s %10s %10s %22s %10s  %s",
                      group_col, value_label, "SE", "95% CI", "p", "notes")
  } else {
    header <- sprintf("  %-25s %10s %10s %22s %10s",
                      "contrast", value_label, "SE", "95% CI", "p")
  }
  sep <- paste0("  ", strrep("-", nchar(header) - 2))
  body <- vapply(seq_len(nrow(ct)), function(i) {
    ci <- sprintf("[%s, %s]", fmt(ct$lower.CL[i]), fmt(ct$upper.CL[i]))
    label <- if (has_group) as.character(ct[[group_col]][i]) else
              as.character(ct$contrast[i])
    if (has_group && show_n) {
      sprintf("  %-12s %10s %10s %22s %10s  %6s %6s  %s",
              label, fmt(ct$estimate[i]), fmt(ct$SE[i]),
              ci, fmt_p(ct$p.value[i]),
              as.character(ct$n_control[i] %||% "NA"),
              as.character(ct$n_warmed[i] %||% "NA"),
              as.character(notes[i] %||% ""))
    } else if (has_group) {
      sprintf("  %-12s %10s %10s %22s %10s  %s",
              label, fmt(ct$estimate[i]), fmt(ct$SE[i]),
              ci, fmt_p(ct$p.value[i]),
              as.character(notes[i] %||% ""))
    } else {
      sprintf("  %-25s %10s %10s %22s %10s",
              label, fmt(ct$estimate[i]), fmt(ct$SE[i]),
              ci, fmt_p(ct$p.value[i]))
    }
  }, character(1))
  paste(c(header, sep, body, ""), collapse = "\n")
}

# ============================================================================
# 1. DATA IMPORT
# ============================================================================
cat(">>> Loading data ...\n")

# CLIMATE.csv: re-level season to Dry-then-Wet for chronological ordering.
clim <- read_csv(clim_path, show_col_types = FALSE) %>%
  mutate(
    Year   = as.integer(Year),
    Month  = as.integer(Month),
    cwd    = as.numeric(cwd),
    prec   = as.numeric(prec),
    season = factor(season, levels = c("Dry", "Wet"))
  )

# T_VWC.csv: source-of-truth for temperature and VWC. Flux column may be
# absent here; if present it is ignored. No flux-based filter.
tv_raw <- read_csv(tv_path, show_col_types = FALSE) %>%
  mutate(
    DayHour     = ymd_hms(DayHour, quiet = TRUE, tz = "UTC"),
    Date        = as.Date(DayHour),
    Plot        = factor(Plot),
    Treatment   = factor(Treatment, levels = c("control", "warmed")),
    Temperature = as.numeric(Temperature),
    VWC         = as.numeric(VWC)
  )

# FLUX.csv: source for the flux pipeline. Flux > 0 filter applied as standard
# QC for log-transformed flux.
flux_raw <- read_csv(flux_path, show_col_types = FALSE) %>%
  mutate(
    DayHour     = ymd_hms(DayHour, quiet = TRUE, tz = "UTC"),
    Date        = as.Date(DayHour),
    Plot        = factor(Plot),
    Treatment   = factor(Treatment, levels = c("control", "warmed")),
    Flux        = as.numeric(Flux),
    Temperature = as.numeric(Temperature),
    VWC         = as.numeric(VWC)
  )

# ============================================================================
# 2. DAILY AGGREGATION
# ============================================================================
# Temp/VWC daily aggregate: filter on temp/VWC presence, n_hours >= 12.
# No Flux filter — T_VWC.csv is the source-of-truth.
cat(">>> Building daily temp/VWC and flux datasets ...\n")

d_tv <- tv_raw %>%
  filter(!is.na(Date), !is.na(Temperature), !is.na(VWC)) %>%
  group_by(Plot, Treatment, Date) %>%
  summarise(
    Year      = year(first(Date)),
    Month     = month(first(Date)),
    temp_mean = mean(Temperature, na.rm = TRUE),
    vwc_mean  = mean(VWC,         na.rm = TRUE),
    n_hours   = dplyr::n(),
    .groups   = "drop"
  ) %>%
  filter(n_hours >= 12) %>%
  left_join(clim, by = c("Year", "Month")) %>%
  filter(!is.na(season)) %>%
  filter(Date >= date_filter) %>%
  arrange(Plot, Date) %>%
  mutate(Year_f = factor(Year))

# Flux daily aggregate: Flux > 0 retained for log transformation.
d_fl <- flux_raw %>%
  filter(!is.na(Date), !is.na(Flux), Flux > 0) %>%
  group_by(Plot, Treatment, Date) %>%
  summarise(
    Year      = year(first(Date)),
    Month     = month(first(Date)),
    flux_mean = mean(Flux,        na.rm = TRUE),
    temp_mean = mean(Temperature, na.rm = TRUE),
    vwc_mean  = mean(VWC,         na.rm = TRUE),
    n_hours   = dplyr::n(),
    .groups   = "drop"
  ) %>%
  filter(n_hours >= 12) %>%
  left_join(clim, by = c("Year", "Month")) %>%
  filter(!is.na(season)) %>%
  filter(Date >= date_filter) %>%
  arrange(Plot, Date) %>%
  mutate(log_flux = log(flux_mean),
         Year_f   = factor(Year))

# ============================================================================
# 3. AR(1) COORDINATE AND PLOT_YEAR GROUPING
# ============================================================================
# time_ou = numFactor(day_num) indexes the AR(1) process. Plot_Year restricts
# the AR(1) to within-year so it cannot drift across years and absorb year-
# level fixed-effect variation. Day numbering anchored at the earliest date
# across both datasets so flux and temp/VWC share a coordinate system.
min_date_all <- min(c(d_tv$Date, d_fl$Date), na.rm = TRUE)

augment_for_glmm <- function(df) {
  df %>%
    mutate(
      day_num   = as.integer(Date - min_date_all) + 1L,
      time_ou   = glmmTMB::numFactor(day_num),
      Plot_Year = factor(interaction(Plot, Year_f, drop = TRUE)),
      season    = factor(season, levels = c("Dry", "Wet"))
    ) %>%
    droplevels()
}
d_analysis <- augment_for_glmm(d_tv) %>%
  mutate(vwc_pct = vwc_mean * 100)   # VWC scaled for numerical stability
d_flux     <- augment_for_glmm(d_fl)

validate_analysis_frame <- function(data, label) {
  if (nrow(data) == 0L) stop(label, " has no eligible plot-days after filtering.")
  if (n_distinct(data$Plot) < 2L) stop(label, " has fewer than two plots.")
  if (!all(c("control", "warmed") %in% as.character(unique(data$Treatment)))) {
    stop(label, " must contain both control and warmed observations.")
  }
  invisible(data)
}
validate_analysis_frame(d_analysis, "Temp/VWC analysis data")
validate_analysis_frame(d_flux, "Flux analysis data")

# Grand-mean centring constants for the mechanistic covariates (computed on the
# flux analysis frame, the data the mechanistic models are fitted to). Retained
# so derived quantities can be back-transformed to real units:
#   Q10 uses slopes (centring-invariant); VWC optimum vertex + VWC_CENTRE.
TEMP_CENTRE <- mean(d_flux$temp_mean, na.rm = TRUE)
VWC_CENTRE  <- mean(d_flux$vwc_mean,  na.rm = TRUE)
d_flux <- d_flux %>%
  mutate(temp_c = temp_mean - TEMP_CENTRE,
         vwc_c  = vwc_mean  - VWC_CENTRE)

# ============================================================================
# 4. COVERAGE DIAGNOSTICS
# ============================================================================
coverage_tv_year_season <- d_analysis %>%
  count(Plot, Treatment, Year, season, name = "n_days")
coverage_fl_year_season <- d_flux %>%
  count(Plot, Treatment, Year, season, name = "n_days")

coverage_tv_year <- d_analysis %>%
  count(Plot, Treatment, Year, name = "n_days") %>%
  mutate(low_coverage = n_days < MIN_PLOT_YEAR_DAYS)
coverage_fl_year <- d_flux %>%
  count(Plot, Treatment, Year, name = "n_days") %>%
  mutate(low_coverage = n_days < MIN_PLOT_YEAR_DAYS)

cell_n_tv <- build_cell_n(d_analysis)
cell_n_fl <- build_cell_n(d_flux)
cell_n_tv_year <- build_cell_n(d_analysis, "Year_f")
cell_n_fl_year <- build_cell_n(d_flux, "Year_f")
cell_n_tv_season <- build_cell_n(d_analysis, "season")
cell_n_fl_season <- build_cell_n(d_flux, "season")

write_csv(coverage_tv_year_season,
          file.path(out_dir, "tables", "coverage_tv_year_season.csv"))
write_csv(coverage_fl_year_season,
          file.path(out_dir, "tables", "coverage_flux_year_season.csv"))
write_csv(coverage_tv_year,
          file.path(out_dir, "tables", "coverage_tv_year.csv"))
write_csv(coverage_fl_year,
          file.path(out_dir, "tables", "coverage_flux_year.csv"))
write_csv(cell_n_tv, file.path(out_dir, "tables", "cell_n_tv.csv"))
write_csv(cell_n_fl, file.path(out_dir, "tables", "cell_n_flux.csv"))
write_csv(cell_n_tv_year, file.path(out_dir, "tables", "cell_n_tv_year.csv"))
write_csv(cell_n_fl_year, file.path(out_dir, "tables", "cell_n_flux_year.csv"))
write_csv(cell_n_tv_season, file.path(out_dir, "tables", "cell_n_tv_season.csv"))
write_csv(cell_n_fl_season, file.path(out_dir, "tables", "cell_n_flux_season.csv"))

# ============================================================================
# 5. MODEL FITTING
# ============================================================================
cat(">>> Fitting primary models (dispformula ~1; dispersed-start retry) ...\n")

# Common random-structure expression.
rs_term <- "(1 | Plot) + ar1(time_ou + 0 | Plot_Year)"

# All models: dispformula ~1 (uniform; dispersion does not affect the fixed
# effects and ~1 gives certified convergence across the set). Observed-level
# ar1 via time_ou. Pooled models give the balanced long-term marginal warming
# effect (Treatment averaged over Year and season); interaction models give the
# year- and season-resolved effects. Mechanistic pathway (common slopes) gives
# the flux-environment associations; mechanistic modification (Treatment x
# covariate) tests whether warming shifts the response surface (Q10, VWC optimum).

# --- Temperature ---
cat("    temperature (pooled) ...\n")
temp_pooled_res <- fit_with_rescue(
  formula     = as.formula(paste("temp_mean ~ Treatment + Year_f + season +", rs_term)),
  dispformula = ~ 1, data = d_analysis, label = "temp_pooled")
cat("    temperature (interaction) ...\n")
temp_inter_res <- fit_with_rescue(
  formula     = as.formula(paste("temp_mean ~ Treatment * Year_f + season +", rs_term)),
  dispformula = ~ 1, data = d_analysis, label = "temp_inter")

# --- VWC ---
cat("    VWC (pooled) ...\n")
vwc_pooled_res <- fit_with_rescue(
  formula     = as.formula(paste("vwc_pct ~ Treatment + Year_f + season +", rs_term)),
  dispformula = ~ 1, data = d_analysis, label = "vwc_pooled")
cat("    VWC (interaction) ...\n")
vwc_inter_res <- fit_with_rescue(
  formula     = as.formula(paste("vwc_pct ~ Treatment * Year_f + season +", rs_term)),
  dispformula = ~ 1, data = d_analysis, label = "vwc_inter")

# --- Flux ---
cat("    flux (pooled) ...\n")
flux_pooled_res <- fit_with_rescue(
  formula     = as.formula(paste("log_flux ~ Treatment + Year_f + season +", rs_term)),
  dispformula = ~ 1, data = d_flux, label = "flux_pooled")
cat("    flux (interaction) ...\n")
flux_inter_res <- fit_with_rescue(
  formula     = as.formula(paste("log_flux ~ Treatment * Year_f + Treatment * season +", rs_term)),
  dispformula = ~ 1, data = d_flux, label = "flux_inter")

# --- Mechanistic (centred covariates) ---
cat("    mechanistic pathway (common slopes) ...\n")
mech_path_res <- fit_with_rescue(
  formula = as.formula(paste(
    "log_flux ~ Treatment + Year_f + season + temp_c + vwc_c + I(vwc_c^2) +", rs_term)),
  dispformula = ~ 1, data = d_flux, label = "mech_path")
cat("    mechanistic modification (Treatment x covariate) ...\n")
mech_mod_res <- fit_with_rescue(
  formula = as.formula(paste(
    "log_flux ~ Treatment + Year_f + season + temp_c + vwc_c + I(vwc_c^2)",
    "+ Treatment:temp_c + Treatment:vwc_c + Treatment:I(vwc_c^2) +", rs_term)),
  dispformula = ~ 1, data = d_flux, label = "mech_mod")

# ============================================================================
# 6. CONTRASTS
# ============================================================================
cat(">>> Extracting contrasts (unadjusted p-values) ...\n")

# Overall long-term marginal warmed - control: from the POOLED models
# (Treatment averaged over Year and season = the balanced long-term effect).
temp_marg <- extract_contrast(temp_pooled_res$fit, ~ Treatment)
vwc_marg  <- extract_contrast(vwc_pooled_res$fit,  ~ Treatment)
flux_marg <- extract_contrast(flux_pooled_res$fit, ~ Treatment)

# Per-Year and per-season warmed - control: from the INTERACTION models
# (the year- and season-resolved effects the interactions exist to estimate).
temp_by_year <- extract_contrast(temp_inter_res$fit, ~ Treatment | Year_f,
                                 cell_n_tv_year, "Year_f")
vwc_by_year <- extract_contrast(vwc_inter_res$fit, ~ Treatment | Year_f,
                                cell_n_tv_year, "Year_f")
flux_by_year <- extract_contrast(flux_inter_res$fit, ~ Treatment | Year_f,
                                 cell_n_fl_year, "Year_f")
# Temperature and VWC do not include Treatment:season, so no model-derived
# treatment-by-season contrasts are reported for those responses.
flux_by_season <- extract_contrast(flux_inter_res$fit, ~ Treatment | season,
                                   cell_n_fl_season, "season")

# Write contrast tables.
write_csv(temp_marg,      file.path(out_dir, "tables", "contrast_temp_marginal.csv"))
write_csv(vwc_marg,       file.path(out_dir, "tables", "contrast_vwc_marginal.csv"))
write_csv(flux_marg,      file.path(out_dir, "tables", "contrast_flux_marginal.csv"))
write_csv(temp_by_year,   file.path(out_dir, "tables", "contrast_temp_by_year.csv"))
write_csv(vwc_by_year,    file.path(out_dir, "tables", "contrast_vwc_by_year.csv"))
write_csv(flux_by_year,   file.path(out_dir, "tables", "contrast_flux_by_year.csv"))
write_csv(flux_by_season, file.path(out_dir, "tables", "contrast_flux_by_season.csv"))

# ============================================================================
# 7. Q10 AND VWC OPTIMUM FROM MECHANISTIC MODELS
# ============================================================================
cat(">>> Deriving Q10 and VWC optima from mechanistic fits ...\n")

# Q10 and VWC optima (treatment-specific response surface) come from the
# MODIFICATION model (Treatment x covariate). VWC optimum is back-transformed
# to the real fraction scale using VWC_CENTRE.
q10_tbl     <- derive_q10(mech_mod_res$fit)
vwc_opt_tbl <- derive_vwc_optimum(mech_mod_res$fit, VWC_CENTRE)

# Mechanistic residual warming effect (warming effect after temperature and
# moisture are accounted for, common slopes) comes from the PATHWAY model's
# Treatment marginal. On log_flux, so exp(contrast) is a geometric-mean ratio.
mech_marg   <- extract_contrast(mech_path_res$fit, ~ Treatment)
if (!"status" %in% names(mech_marg)) {
  mech_marg <- mech_marg %>%
    mutate(ratio        = exp(estimate),
           ratio_lower  = exp(lower.CL),
           ratio_upper  = exp(upper.CL))
}

write_csv(q10_tbl,     file.path(out_dir, "tables", "q10_by_treatment.csv"))
write_csv(vwc_opt_tbl, file.path(out_dir, "tables", "vwc_optimum_by_treatment.csv"))
write_csv(mech_marg,   file.path(out_dir, "tables", "mech_residual_marginal.csv"))

# ============================================================================
# 7B. DERIVED CONTRAST TESTS AND VARIABILITY
# ----------------------------------------------------------------------------
# Pre-specified follow-on tests not directly produced by emmeans:
#   - Q10 contrast: warmed vs control temperature-slope difference (from the
#     modification model's Treatmentwarmed:temp_c interaction)
#   - VWC optimum contrast: warmed vs control quadratic-peak difference (delta
#     method over the four moisture coefficients; centring cancels in the diff)
#   - Treatment effect on variability: empirical per-(Plot x Year x season) SD
#     of the response modelled as log(SD) ~ Treatment + season + (1|Plot);
#     exp(Treatment coef) = warmed/control SD ratio. Model-independent; replaces
#     the former dispformula-based variance test (dispformula is now ~1).
# ============================================================================
cat(">>> Computing derived contrasts and variability ...\n")

# Q10 contrast: the slope difference IS the Treatmentwarmed:temp_c interaction.
derive_q10_contrast <- function(model) {
  if (is_failed_glmm(model)) return(tibble(status = "model failed"))
  fe <- fixef(model)$cond; ve <- vcov(model)$cond
  nm_i <- "Treatmentwarmed:temp_c"
  if (!nm_i %in% names(fe)) return(tibble(status = "interaction slope not found"))
  d  <- as.numeric(fe[nm_i]); se_d <- sqrt(ve[nm_i, nm_i]); zv <- d / se_d
  tibble(
    contrast      = "warmed - control (slope of log_flux on temperature)",
    delta_slope   = d,
    SE            = as.numeric(se_d),
    delta_lo      = d - 1.96 * se_d,
    delta_hi      = d + 1.96 * se_d,
    Q10_ratio     = exp(10 * d),
    Q10_ratio_lo  = exp(10 * (d - 1.96 * se_d)),
    Q10_ratio_hi  = exp(10 * (d + 1.96 * se_d)),
    z_value       = as.numeric(zv),
    p_value       = 2 * pnorm(-abs(zv))
  )
}

# VWC optimum contrast under the centred-interaction parameterisation.
# control: B1m=vwc_c, B2m=I(vwc_c^2). warmed: B1w=B1m+Tw:vwc_c, B2w=B2m+Tw:I(vwc_c^2).
# delta = opt_w - opt_c (centring constant cancels). Gradient over
# (vwc_c, I(vwc_c^2), Tw:vwc_c, Tw:I(vwc_c^2)); verified vs numerical grad.
derive_vwc_opt_contrast <- function(model) {
  if (is_failed_glmm(model)) return(tibble(status = "model failed"))
  fe <- fixef(model)$cond; ve <- vcov(model)$cond
  b1m <- "vwc_c"; b2m <- "I(vwc_c^2)"
  b1i <- "Treatmentwarmed:vwc_c"; b2i <- "Treatmentwarmed:I(vwc_c^2)"
  if (!all(c(b1m, b2m, b1i, b2i) %in% names(fe)))
    return(tibble(status = "coefficient names missing"))
  B1m <- fe[b1m]; B2m <- fe[b2m]; B1i <- fe[b1i]; B2i <- fe[b2i]
  B1w <- B1m + B1i; B2w <- B2m + B2i
  if (B2m >= 0 || B2w >= 0)
    return(tibble(status = sprintf("non-concave (b2c = %.3f, b2w = %.3f); contrast undefined",
                                   B2m, B2w)))
  opt_c <- as.numeric(-B1m / (2 * B2m))
  opt_w <- as.numeric(-B1w / (2 * B2w))
  d <- opt_w - opt_c
  g <- c(-1/(2*B2w) + 1/(2*B2m),                 # d/d vwc_c
         B1w/(2*B2w^2) - B1m/(2*B2m^2),          # d/d I(vwc_c^2)
         -1/(2*B2w),                             # d/d Tw:vwc_c
         B1w/(2*B2w^2))                          # d/d Tw:I(vwc_c^2)
  V  <- ve[c(b1m, b2m, b1i, b2i), c(b1m, b2m, b1i, b2i)]
  se_d <- sqrt(max(as.numeric(t(g) %*% V %*% g), 0)); zv <- d / se_d
  tibble(
    contrast    = "warmed - control (VWC optimum, fraction scale)",
    opt_control = opt_c + VWC_CENTRE, opt_warmed = opt_w + VWC_CENTRE,
    delta_opt   = d, SE = as.numeric(se_d),
    delta_lo    = d - 1.96 * se_d, delta_hi = d + 1.96 * se_d,
    z_value     = as.numeric(zv), p_value = 2 * pnorm(-abs(zv)),
    note        = ""
  )
}

# Exploratory empirical-SD variability. The finite-sample expectation of log(S)
# depends on n_days; correct that bias before comparing cells. This does not
# remove uncertainty from estimating each cell SD, so results remain descriptive.
empirical_sd_variability <- function(data, response, label, MIN_CELL_N = 5L) {
  cells <- data %>%
    group_by(Plot, Treatment, Year_f, season) %>%
    summarise(cell_sd = sd(.data[[response]], na.rm = TRUE),
              n_days  = sum(!is.na(.data[[response]])), .groups = "drop") %>%
    filter(n_days >= MIN_CELL_N, is.finite(cell_sd), cell_sd > 0)
  if (nrow(cells) < 4 || n_distinct(cells$Treatment) < 2)
    return(tibble(model = label, status = "insufficient cells"))
  cells <- cells %>%
    mutate(log_sd_bias = 0.5 * (digamma((n_days - 1) / 2) + log(2) - log(n_days - 1)),
           log_sd = log(cell_sd) - log_sd_bias)
  m <- tryCatch(glmmTMB(log_sd ~ Treatment + season + (1 | Plot),
                        data = cells, family = gaussian(), REML = TRUE,
                        control = ctrl_default), error = function(e) e)
  if (is_failed_glmm(m)) return(tibble(model = label, status = "SD model failed"))
  s <- summary(m)$coefficients$cond; rn <- rownames(s)
  tw <- rn[grepl("^Treatmentwarmed$", rn)]
  if (length(tw) == 0) return(tibble(model = label, status = "Treatment coef not found"))
  est <- s[tw, "Estimate"]; se <- s[tw, "Std. Error"]; p <- s[tw, "Pr(>|z|)"]
  gm <- cells %>% group_by(Treatment) %>%
    summarise(gm_sd = exp(mean(log_sd)), .groups = "drop")
  raw_ratio <- gm$gm_sd[gm$Treatment == "warmed"] / gm$gm_sd[gm$Treatment == "control"]
  tibble(
    model           = label,
    n_cells         = nrow(cells),
    logSD_diff      = as.numeric(est),
    SE              = as.numeric(se),
    SD_ratio        = exp(as.numeric(est)),
    SD_ratio_lo     = exp(as.numeric(est - 1.96 * se)),
    SD_ratio_hi     = exp(as.numeric(est + 1.96 * se)),
    p_value         = as.numeric(p),
    raw_gm_SD_ratio = as.numeric(raw_ratio),
    status          = ""
  )
}

q10_contrast     <- derive_q10_contrast(mech_mod_res$fit)
vwc_opt_contrast <- derive_vwc_opt_contrast(mech_mod_res$fit)
variability_summary <- bind_rows(
  empirical_sd_variability(d_analysis, "temp_mean", "temperature"),
  empirical_sd_variability(d_analysis, "vwc_pct",   "vwc"),
  empirical_sd_variability(d_flux,     "log_flux",  "flux (log)")
)

write_csv(q10_contrast,        file.path(out_dir, "tables", "q10_contrast.csv"))
write_csv(vwc_opt_contrast,    file.path(out_dir, "tables", "vwc_optimum_contrast.csv"))
write_csv(variability_summary, file.path(out_dir, "tables", "variability_summary.csv"))

# ============================================================================
# 8. COEFFICIENT TABLES
# ============================================================================
compact_coef <- function(fit, label) {
  if (is_failed_glmm(fit)) return(tibble(model = label, status = "failed"))
  s <- summary(fit)
  cond <- as.data.frame(s$coefficients$cond) %>%
    rownames_to_column("term") %>%
    mutate(component = "conditional", model = label, .before = 1)
  disp <- if (!is.null(s$coefficients$disp))
    as.data.frame(s$coefficients$disp) %>%
      rownames_to_column("term") %>%
      mutate(component = "dispersion", model = label, .before = 1) else NULL
  bind_rows(cond, disp)
}
write_csv(compact_coef(temp_pooled_res$fit, "temp_pooled"),
          file.path(out_dir, "tables", "coef_temp_pooled.csv"))
write_csv(compact_coef(temp_inter_res$fit,  "temp_inter"),
          file.path(out_dir, "tables", "coef_temp_inter.csv"))
write_csv(compact_coef(vwc_pooled_res$fit,  "vwc_pooled"),
          file.path(out_dir, "tables", "coef_vwc_pooled.csv"))
write_csv(compact_coef(vwc_inter_res$fit,   "vwc_inter"),
          file.path(out_dir, "tables", "coef_vwc_inter.csv"))
write_csv(compact_coef(flux_pooled_res$fit, "flux_pooled"),
          file.path(out_dir, "tables", "coef_flux_pooled.csv"))
write_csv(compact_coef(flux_inter_res$fit,  "flux_inter"),
          file.path(out_dir, "tables", "coef_flux_inter.csv"))
write_csv(compact_coef(mech_path_res$fit,   "mech_pathway"),
          file.path(out_dir, "tables", "coef_mech_pathway.csv"))
write_csv(compact_coef(mech_mod_res$fit,    "mech_modification"),
          file.path(out_dir, "tables", "coef_mech_modification.csv"))

# ============================================================================
# 9. PLAIN-TEXT REPORT
# ============================================================================
cat(">>> Writing report ...\n")

# Per-year / per-season comparison counts (dynamic; reported in headings).
n_yr_temp <- if (!"status" %in% names(temp_by_year))   nrow(temp_by_year)   else 0L
n_yr_vwc  <- if (!"status" %in% names(vwc_by_year))    nrow(vwc_by_year)    else 0L
n_yr_flux <- if (!"status" %in% names(flux_by_year))   nrow(flux_by_year)   else 0L
n_se_flux <- if (!"status" %in% names(flux_by_season)) nrow(flux_by_season) else 0L

rule <- function(c = "=", n = 75) strrep(c, n)
fit_block <- function(res, label) {
  if (is.null(res$fit)) {
    return(c(sprintf("  Status:           FAILED (%s)", res$tactic), ""))
  }
  rs <- extract_random_summary(res$fit)
  c(
    sprintf("  Status:           CONVERGED (%s)", res$tactic),
    sprintf("  AIC:              %s", fmt(safe_AIC(res$fit), 1)),
    sprintf("  pdHess:           %s",
            if (is_pdHess(res$fit)) "TRUE" else "FALSE (caveat: SEs approximate)"),
    sprintf("  sigma(Plot):      %s", fmt(rs$sigma_plot,      4)),
    sprintf("  sigma(Plot_Year): %s", fmt(rs$sigma_plot_year, 4)),
    sprintf("  phi (lag-1):      %s", fmt(rs$phi,             4)),
    sprintf("  tau (steps):      %s", fmt(rs$tau,             2)),
    ""
  )
}

report_lines <- character(0)
add <- function(...) report_lines <<- c(report_lines, ...)

add(
  rule("="),
  "TRACE LUQUILLO SOIL CO2 WARMING - GLMM v3 FINAL REPORT",
  rule("="),
  sprintf("Generated:    %s", Sys.time()),
  sprintf("Working dir:  %s", getwd()),
  sprintf("Output:       %s", out_dir),
  sprintf("Date filter:  Date >= %s", date_filter),
  ""
)

# Data summary
add(
  rule("-"),
  "DATA SUMMARY",
  rule("-"),
  sprintf("  T_VWC.csv rows (hourly):  %s", format(nrow(tv_raw),  big.mark=",")),
  sprintf("  FLUX.csv  rows (hourly):  %s", format(nrow(flux_raw), big.mark=",")),
  sprintf("  CLIMATE.csv rows:         %s", format(nrow(clim),    big.mark=",")),
  "",
  sprintf("  d_analysis (temp/VWC) plot-days: %s",
          format(nrow(d_analysis), big.mark=",")),
  sprintf("  d_flux               plot-days:  %s",
          format(nrow(d_flux),     big.mark=",")),
  sprintf("  Year range (temp/VWC):  %d - %d",
          min(d_analysis$Year), max(d_analysis$Year)),
  sprintf("  Year range (flux):      %d - %d",
          min(d_flux$Year),     max(d_flux$Year)),
  sprintf("  Plots:                  %s",
          paste(sort(unique(as.character(d_analysis$Plot))), collapse = ", ")),
  ""
)

# Coverage tables
coverage_block <- function(cell_n, label) {
  pv <- cell_n %>%
    mutate(yr_seas = paste(Year_f, season, sep = "/")) %>%
    select(yr_seas, n_control, n_warmed, n_min) %>%
    arrange(yr_seas)
  out <- c(sprintf("  %s coverage by (Year x season) cell:", label),
           sprintf("  %-10s %10s %10s %10s",
                   "yr_seas", "n_control", "n_warmed", "n_min"),
           paste0("  ", strrep("-", 44)))
  for (i in seq_len(nrow(pv))) {
    out <- c(out, sprintf("  %-10s %10d %10d %10d",
                          pv$yr_seas[i], pv$n_control[i],
                          pv$n_warmed[i], pv$n_min[i]))
  }
  out <- c(out, "")
  out
}
add(rule("-"), "COVERAGE", rule("-"))
add(coverage_block(cell_n_tv, "Temp/VWC"))
add(coverage_block(cell_n_fl, "Flux"))

# Per-response sections
add(rule("-"),
    "TEMPERATURE (response = temp_mean, deg C)",
    rule("-"),
    "  Pooled:      temp_mean ~ Treatment + Year_f + season",
    "  Interaction: temp_mean ~ Treatment * Year_f + season",
    "               + (1|Plot) + ar1(time_ou + 0 | Plot_Year); dispformula = ~ 1",
    "")
add("  [Pooled model - overall long-term effect]")
add(fit_block(temp_pooled_res, "temp_pooled"))
add("  Marginal warmed - control (deg C):")
add(contrast_to_text(temp_marg, value_label = "deg C"))
add("  [Interaction model - year/season resolution]")
add(fit_block(temp_inter_res, "temp_inter"))
add(sprintf("  Per-year warmed - control (%d comparisons; p unadjusted):", n_yr_temp))
add(contrast_to_text(temp_by_year, group_col = "Year_f",
                     value_label = "deg C", show_n = FALSE))

add(rule("-"),
    "VWC (response = vwc_pct = vwc_mean * 100, percent)",
    rule("-"),
    "  Pooled:      vwc_pct ~ Treatment + Year_f + season",
    "  Interaction: vwc_pct ~ Treatment * Year_f + season",
    "               + (1|Plot) + ar1(time_ou + 0 | Plot_Year); dispformula = ~ 1",
    "",
    "  Note: VWC scaled to percent for numerical stability. Divide by 100",
    "        for original VWC fraction units.",
    "")
add("  [Pooled model - overall long-term effect]")
add(fit_block(vwc_pooled_res, "vwc_pooled"))
add("  Marginal warmed - control (vwc_pct):")
add(contrast_to_text(vwc_marg, value_label = "vwc_pct"))
add("  [Interaction model - year/season resolution]")
add(fit_block(vwc_inter_res, "vwc_inter"))
add(sprintf("  Per-year warmed - control (vwc_pct; %d comparisons; p unadjusted):", n_yr_vwc))
add(contrast_to_text(vwc_by_year, group_col = "Year_f",
                     value_label = "vwc_pct"))

add(rule("-"),
    "FLUX (response = log_flux)",
    rule("-"),
    "  Pooled:      log_flux ~ Treatment + Year_f + season   [long-term summary;",
    "               a single pooled effect averages over the episodic year",
    "               structure - interpret alongside the interaction model]",
    "  Interaction: log_flux ~ Treatment * Year_f + Treatment * season",
    "               + (1|Plot) + ar1(time_ou + 0 | Plot_Year); dispformula = ~ 1",
    "")
add("  [Pooled model - overall long-term effect]")
add(fit_block(flux_pooled_res, "flux_pooled"))
add("  [Interaction model - year/season resolution]")
add(fit_block(flux_inter_res, "flux_inter"))

# Flux contrasts: report log-scale and geometric-mean ratios.
add_flux_contrast <- function(ct, group_col, value_label) {
  if ("status" %in% names(ct)) {
    add(paste0("  [", paste(ct$status, collapse = "; "), "]\n"))
    return(invisible())
  }
  ct2 <- ct %>%
    mutate(ratio_w_c   = exp(estimate),
           ratio_lo    = exp(lower.CL),
           ratio_hi    = exp(upper.CL))
  add("  Log-scale contrasts (warmed - control on log_flux):")
  add(contrast_to_text(ct, group_col = group_col, value_label = value_label))
  add("  Back-transformed geometric-mean ratios (warmed / control):")
  if (is.null(group_col)) {
    add(sprintf("    ratio = %s   95%% CI [%s, %s]",
                fmt(ct2$ratio_w_c, 3), fmt(ct2$ratio_lo, 3),
                fmt(ct2$ratio_hi, 3)))
  } else {
    for (i in seq_len(nrow(ct2))) {
      add(sprintf("    %s : ratio = %s   95%% CI [%s, %s]",
                  ct2[[group_col]][i],
                  fmt(ct2$ratio_w_c[i], 3),
                  fmt(ct2$ratio_lo[i], 3),
                  fmt(ct2$ratio_hi[i], 3)))
    }
  }
  add("")
}
add("  Marginal warmed - control:")
add_flux_contrast(flux_marg, NULL, "log_flux")
add(sprintf("  Per-year warmed - control (%d comparisons; p unadjusted):", n_yr_flux))
add_flux_contrast(flux_by_year, "Year_f", "log_flux")
add(sprintf("  Per-season warmed - control (%d comparisons; p unadjusted):", n_se_flux))
add_flux_contrast(flux_by_season, "season", "log_flux")

add(rule("-"),
    "MECHANISTIC FLUX MODELS (Q10 + VWC optimum)",
    rule("-"),
    "  Pathway:      log_flux ~ Treatment + Year_f + season",
    "                + temp_c + vwc_c + I(vwc_c^2)   [common slopes]",
    "  Modification: + Treatment:temp_c + Treatment:vwc_c + Treatment:I(vwc_c^2)",
    "                + (1|Plot) + ar1(time_ou + 0 | Plot_Year); dispformula = ~ 1",
    sprintf("  Covariates centred: temp_c on %s deg C, vwc_c on %s (fraction).",
            fmt(TEMP_CENTRE, 3), fmt(VWC_CENTRE, 4)),
    "  Q10 and VWC optima are from the MODIFICATION model (treatment-specific",
    "  surface); the residual warming effect is from the PATHWAY model.",
    "")
add("  [Modification model - treatment-specific response surface]")
add(fit_block(mech_mod_res, "mech_modification"))

# Q10 results
add("  Q10 by treatment (Q10 = exp(10 * slope) of log_flux on temperature):")
if ("status" %in% names(q10_tbl)) {
  add(paste0("  [", paste(unique(q10_tbl$status), collapse = "; "), "]"))
} else {
  add(sprintf("  %-10s %10s %10s %10s %10s",
              "treatment", "slope", "Q10", "Q10_lo", "Q10_hi"))
  add(paste0("  ", strrep("-", 55)))
  for (i in seq_len(nrow(q10_tbl))) {
    add(sprintf("  %-10s %10s %10s %10s %10s",
                q10_tbl$treatment[i],
                fmt(q10_tbl$slope[i],  5),
                fmt(q10_tbl$Q10[i],    3),
                fmt(q10_tbl$Q10_lo[i], 3),
                fmt(q10_tbl$Q10_hi[i], 3)))
  }
}
add("")

# VWC optimum (back-transformed to real fraction scale via VWC_CENTRE)
add("  VWC optimum by treatment (peak of quadratic, real VWC fraction scale):")
add(sprintf("  %-10s %10s %10s %10s  %s",
            "treatment", "optimum", "opt_lo", "opt_hi", "note"))
add(paste0("  ", strrep("-", 58)))
for (i in seq_len(nrow(vwc_opt_tbl))) {
  add(sprintf("  %-10s %10s %10s %10s  %s",
              vwc_opt_tbl$treatment[i],
              fmt(vwc_opt_tbl$optimum[i], 4),
              fmt(vwc_opt_tbl$opt_lo[i],  4),
              fmt(vwc_opt_tbl$opt_hi[i],  4),
              if ("note" %in% names(vwc_opt_tbl)) vwc_opt_tbl$note[i] else ""))
}
add("")

# Mechanistic residual (from PATHWAY model: warming effect after T+VWC absorbed)
add("  Residual warmed-control on log_flux, from PATHWAY model",
    "  (warming effect after temperature and moisture absorbed, common slopes;",
    "  exp(contrast) is a geometric-mean ratio):")
if (!"status" %in% names(mech_marg)) {
  add(sprintf("    estimate = %s   SE = %s   95%% CI [%s, %s]   p = %s",
              fmt(mech_marg$estimate[1]),  fmt(mech_marg$SE[1]),
              fmt(mech_marg$lower.CL[1]),  fmt(mech_marg$upper.CL[1]),
              fmt_p(mech_marg$p.value[1])))
  add(sprintf("    ratio    = %s   95%% CI [%s, %s]",
              fmt(mech_marg$ratio[1]),
              fmt(mech_marg$ratio_lower[1]),
              fmt(mech_marg$ratio_upper[1])))
}
add("")

# Q10 contrast
add(rule("-"),
    "DERIVED CONTRAST TESTS",
    rule("-"))
add("  Q10 contrast (warmed - control temperature slope on log_flux):")
if (!"status" %in% names(q10_contrast)) {
  add(sprintf("    delta_slope = %s   SE = %s   95%% CI [%s, %s]   p = %s",
              fmt(q10_contrast$delta_slope, 5),
              fmt(q10_contrast$SE,          5),
              fmt(q10_contrast$delta_lo,    5),
              fmt(q10_contrast$delta_hi,    5),
              fmt_p(q10_contrast$p_value)))
  add(sprintf("    Q10 ratio (warmed/control) = %s   95%% CI [%s, %s]",
              fmt(q10_contrast$Q10_ratio,    3),
              fmt(q10_contrast$Q10_ratio_lo, 3),
              fmt(q10_contrast$Q10_ratio_hi, 3)))
}
add("")
add("  VWC optimum contrast (warmed - control):")
if (!"status" %in% names(vwc_opt_contrast)) {
  add(sprintf("    opt_control = %s   opt_warmed = %s",
              fmt(vwc_opt_contrast$opt_control, 4),
              fmt(vwc_opt_contrast$opt_warmed,  4)))
  add(sprintf("    delta_opt   = %s   SE = %s   95%% CI [%s, %s]   p = %s",
              fmt(vwc_opt_contrast$delta_opt, 4),
              fmt(vwc_opt_contrast$SE,        4),
              fmt(vwc_opt_contrast$delta_lo,  4),
              fmt(vwc_opt_contrast$delta_hi,  4),
              fmt_p(vwc_opt_contrast$p_value)))
}
add("")

# Variability (empirical per-cell SD; model-independent)
add(rule("-"),
    "EXPLORATORY TREATMENT EFFECT ON VARIABILITY (EMPIRICAL PER-CELL SD)",
    rule("-"),
    "  Per-(Plot x Year x season) SD of the response (scatter around each",
    "  cell mean), with a finite-sample log-SD bias correction, then modelled",
    "  as log(SD) ~ Treatment + season + (1|Plot). This is exploratory because",
    "  cell-SD estimation uncertainty and temporal dependence are not modelled.",
    "  SD_ratio = exp(Treatment coef) = warmed/control SD ratio (adjusted).",
    "  raw_gm_SD_ratio = unadjusted geometric-mean SD ratio. Model-independent",
    "  (does not use the GLMM dispersion structure).",
    "")
add(sprintf("  %-12s %7s %10s %10s %9s %9s %9s %9s",
            "model", "n_cells", "logSD_diff", "SD_ratio", "SD_lo", "SD_hi", "p", "raw_ratio"))
add(paste0("  ", strrep("-", 88)))
for (i in seq_len(nrow(variability_summary))) {
  if ("status" %in% names(variability_summary) &&
      !is.na(variability_summary$status[i]) && variability_summary$status[i] != "") {
    add(sprintf("  %-12s  [%s]", variability_summary$model[i], variability_summary$status[i]))
    next
  }
  add(sprintf("  %-12s %7d %10s %10s %9s %9s %9s %9s",
              variability_summary$model[i],
              variability_summary$n_cells[i],
              fmt(variability_summary$logSD_diff[i],      4),
              fmt(variability_summary$SD_ratio[i],        3),
              fmt(variability_summary$SD_ratio_lo[i],     3),
              fmt(variability_summary$SD_ratio_hi[i],     3),
              fmt_p(variability_summary$p_value[i]),
              fmt(variability_summary$raw_gm_SD_ratio[i], 3)))
}
add("")

# Methodological notes
add(rule("-"),
    "METHODOLOGICAL NOTES",
    rule("-"),
    "  Multiple comparisons: per-year and per-season (flux only) warmed-",
    "  control contrasts are pre-specified questions, each indexing",
    "  a distinct ecological context (year x season cell). p-values are",
    "  reported unadjusted. Comparison counts are stated in the headings of",
    "  the contrast tables.",
    "",
    "  AR(1) within Plot_Year. ar1(time_ou + 0 | Plot) spanning the full",
    "  multi-year record absorbs year-level fixed-effect variation, which",
    "  was empirically shown to distort per-year warming contrasts.",
    "  Restricting the AR(1) to within Plot_Year eliminates this absorption",
    "  while preserving within-year serial correlation handling.",
    "",
    "  Temporal-correlation units. ar1() treats consecutive retained plot-days",
    "  as equally spaced (lag 1) regardless of calendar gaps, so the AR(1)",
    "  decorrelation time tau is expressed in observation steps, not days.",
    "  ou(), which uses true numeric distances and would give tau in days,",
    "  produced non-positive-definite Hessians for the temperature and flux",
    "  primaries under every optimiser tested and was therefore not adopted;",
    "  a uniform ar1() structure is retained across all four models. Treatment",
    "  contrasts were robust to this choice (see model-vs-raw validation).",
    "",
    "  No Treatment:season for temperature and VWC. Including a Treatment x",
    "  season interaction in these two primaries produces a non-positive-",
    "  definite Hessian; the data does not robustly identify the warming-",
    "  by-season interaction for these responses at this sample size. Per-",
    "  season warmed-control contrasts are therefore not reported for temp",
    "  and VWC. Flux retains Treatment * season as primary scientific focus.",
    "",
    "  Mechanistic dispformula = ~ 1. The richer ~ Treatment + season was",
    "  non-identifiable here because temperature and VWC absorb the seasonal",
    "  and treatment mean structure that the dispformula otherwise models",
    "  in the variance.",
    "",
    "  Optimisation: a positive-definite Hessian is required for a certified",
    "  fit. If the standard nlminb fit is finite but has a non-positive-",
    "  nlminb fit is finite but has a non-positive-definite Hessian, the script",
    "  retries from dispersed random-effect starting values. Each fit block",
    "  reports the tactic actually used; an uncertified fit is explicitly",
    "  flagged and its standard errors should not be used for inference.",
    "",
    "  Q10 derivation: Treatment-specific temperature slopes from the",
    "  mechanistic model give Q10 = exp(10 * slope) per treatment. 95% CIs",
    "  from the slope SEs on the log scale, back-transformed.",
    "",
    "  VWC optimum derivation: per-treatment quadratic peak x* = -b1 / (2*b2)",
    "  from Treatment:vwc_mean and Treatment:I(vwc_mean^2). 95% CIs via",
    "  delta method on the relevant 2x2 coefficient sub-block of vcov.",
    "",
    rule("="),
    "END OF REPORT",
    rule("="))

writeLines(report_lines, report_path)

# ============================================================================
# 10. SAVE OBJECTS
# ============================================================================
saveRDS(
  list(
    d_analysis = d_analysis,
    d_flux     = d_flux,
    clim       = clim,
    coverage   = list(
      tv_year_season   = coverage_tv_year_season,
      fl_year_season   = coverage_fl_year_season,
      tv_year          = coverage_tv_year,
      fl_year          = coverage_fl_year,
      cell_n_tv        = cell_n_tv,
      cell_n_fl        = cell_n_fl,
      cell_n_tv_year   = cell_n_tv_year,
      cell_n_fl_year   = cell_n_fl_year,
      cell_n_tv_season = cell_n_tv_season,
      cell_n_fl_season = cell_n_fl_season
    ),
    fits = list(
      temp_pooled = temp_pooled_res,
      temp_inter  = temp_inter_res,
      vwc_pooled  = vwc_pooled_res,
      vwc_inter   = vwc_inter_res,
      flux_pooled = flux_pooled_res,
      flux_inter  = flux_inter_res,
      mech_path   = mech_path_res,
      mech_mod    = mech_mod_res
    ),
    contrasts = list(
      temp_marg      = temp_marg,
      vwc_marg       = vwc_marg,
      flux_marg      = flux_marg,
      mech_marg      = mech_marg,
      temp_by_year   = temp_by_year,
      vwc_by_year    = vwc_by_year,
      flux_by_year   = flux_by_year,
      flux_by_season = flux_by_season
    ),
    derived = list(
      q10_by_treatment         = q10_tbl,
      vwc_optimum_by_treatment = vwc_opt_tbl,
      q10_contrast             = q10_contrast,
      vwc_optimum_contrast     = vwc_opt_contrast,
      variability_summary      = variability_summary
    ),
    config = list(
      date_filter          = date_filter,
      MIN_PLOT_YEAR_DAYS   = MIN_PLOT_YEAR_DAYS,
      EMPTY_CELL_THRESHOLD = EMPTY_CELL_THRESHOLD,
      MAX_RETRY            = MAX_RETRY,
      TEMP_CENTRE          = TEMP_CENTRE,
      VWC_CENTRE           = VWC_CENTRE,
      sessionInfo          = sessionInfo()
    )
  ),
  file.path(out_dir, "objects", "glmm_v3_objects.rds")
)

cat(">>> Done.\n")
cat(sprintf("    Report:  %s\n", normalizePath(report_path, mustWork = FALSE)))
cat(sprintf("    Tables:  %s\n",
            normalizePath(file.path(out_dir, "tables"), mustWork = FALSE)))
cat(sprintf("    Objects: %s\n",
            normalizePath(file.path(out_dir, "objects"), mustWork = FALSE)))
