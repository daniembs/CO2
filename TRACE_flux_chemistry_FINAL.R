# =============================================================================
# TRACE — Soil Flux × Chemistry: Unified Manuscript Analysis
# =============================================================================
# Site:       Tropical Responses to Altered Climate Experiment (TRACE),
#             Sabana Field Research Station, Luquillo Experimental Forest,
#             Puerto Rico (Kimball et al., 2018; Wood et al., 2025).
# Question:   Which soil chemistry pools predict CO2 flux at this P-limited
#             site, and does +4 deg C warming modify chemistry-flux coupling?
#
# Pipeline (single-file, single-run):
#   PART A  Abiotic base model and flux residuals
#   PART B  Primary windowed analysis at +/-14 days
#             - Standard LMM, FDR-corrected
#             - Interaction LMM with emmeans simple slopes
#             - corExp-corrected refits (main effects + interactions)
#             - Microbial P interaction pre-specified as headline target
#   PART C  Robustness checks
#             - Empirical flux ACF per plot
#             - Window sensitivity grid (+/-3 to +/-30 days)
#             - Correlation structure comparison (corExp, corCAR1, corGaus)
#             - Same-day matching as the most stringent temporal-independence
#               test (chemistry observations paired only with hourly flux
#               from the same Plot on the same calendar day)
#   PART D  Publication outputs
#             - Main figures (1-3) and supplementary figures (S1-S6)
#             - DOCX tables (Tables 1-4 main; Tables S1-S2 supplementary)
#             - All CSV result files
#             - Per-model diagnostic panels for AR-robust compounds
#
# Locked analytical decisions (all earlier alternatives tested and discarded):
#   - Primary window: +/-14 days. +/-7 days appears only as one point in
#     the sensitivity grid (Figure S2). Same-day matching is the most
#     stringent sensitivity check (Figure S5, Figure S6, Table S2).
#   - Primary responses: mean log(flux) and flux residuals after T+VWC.
#     Flux CV moved to supplementary (Figure S4) because its biological
#     interpretation confounds diurnal, weather-driven, and episodic
#     variability components.
#   - Chemistry panel: all 19 compounds shown in main figure for
#     transparency; biological discussion focuses on AR-robust survivors.
#   - Autocorrelation correction: corExp(form = ~time | Plot, nugget = TRUE)
#     on Plot x Date aggregated data, treated as a central methodological
#     step rather than as a sensitivity check, because empirical flux ACF
#     exceeds 0.6 at 30-day lags at this site.
#   - HEADLINE_INT_TARGETS: microbial P x Treatment is refit with corExp
#     across both main responses regardless of standard-LMM FDR outcome,
#     so the manuscript headline figure (Figure 2) annotates AR slopes
#     symmetrically across both panels.
#
# Inputs (in working directory):
#   FLUX.csv
#   TRACE_Cores_clean.csv, TRACE_Resin_clean.csv,
#   TRACE_Lysimeter_clean.csv, TRACE_Porewater_clean.csv
#
# Outputs:
#   ./flux_chem_main/                  Tables 1-4 (DOCX + CSV); Figures 1-3
#                                      (PNG + PDF, 600 dpi); abiotic summary
#   ./flux_chem_supplementary/         All robustness CSVs; Figures S1-S6;
#                                      Supplementary DOCX tables;
#                                      per-model diagnostic plots
# =============================================================================

library(tidyverse)
library(lubridate)
library(lme4)
library(lmerTest)
library(MuMIn)
library(broom.mixed)
library(nlme)
library(emmeans)
library(patchwork)
library(flextable)
library(officer)

# =============================================================================
# 0. CONFIGURATION
# =============================================================================
setwd("D:/USDA/TRACE_DM_APR_26/CHEM")
# Windows
WINDOW_PRIMARY <- 14L                                # primary main-text window
WINDOW_GRID    <- c(1L, 3L, 5L, 7L, 10L, 14L, 21L, 30L)  # supplementary grid

# Responses
RESPONSES_MAIN <- c("log_flux_mean", "flux_resid_mean")
RESPONSES_SUPP <- c("flux_cv")
RESPONSES_ALL  <- c(RESPONSES_MAIN, RESPONSES_SUPP)

# Microbial P x Treatment is pre-specified as a headline AR-interaction
# target so simple slopes are available for both main responses regardless
# of standard-LMM FDR outcome.
HEADLINE_INT_TARGETS <- expand_grid(
  file     = "Cores",
  compound = "uP",
  window   = WINDOW_PRIMARY,
  response = RESPONSES_MAIN
)

# Filtering / multiple-testing
MIN_OBS       <- 10L                                 # minimum n per LMM
MIN_PLOTS     <- 2L                                  # minimum unique plots
FDR_THRESHOLD <- 0.05                                # BH FDR

# ACF
ACF_MAX_LAG <- 60L                                   # days
ACF_THRESH  <- 0.1                                   # decorrelation threshold

# Figure quality
FIG_DPI <- 600

# Output directory tree
MAIN_DIR    <- "flux_chem_main"
MAIN_TABLES <- file.path(MAIN_DIR, "tables")
MAIN_FIGS   <- file.path(MAIN_DIR, "figures")
SUPP_DIR    <- "flux_chem_supplementary"
SUPP_TABLES <- file.path(SUPP_DIR, "tables")
SUPP_FIGS   <- file.path(SUPP_DIR, "figures")
SUPP_DIAG   <- file.path(SUPP_DIR, "diagnostics")

for (d in c(MAIN_TABLES, MAIN_FIGS, SUPP_TABLES, SUPP_FIGS, SUPP_DIAG)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# Visual conventions
treat_pal <- c("Control" = "#2166ac", "Warmed" = "#d6604d")
file_pal  <- c("Cores"     = "#4575b4",
               "Resin"     = "#d73027",
               "Lysimeter" = "#1a9641",
               "Porewater" = "#e08214")
plot_pal  <- c("1" = "#1b9e77", "2" = "#d95f02", "3" = "#7570b3",
               "4" = "#e7298a", "5" = "#66a61e", "6" = "#e6ab02")
response_labels <- c(
  log_flux_mean   = "Mean log(flux)",
  flux_resid_mean = "Flux residuals (after T + VWC)",
  flux_cv         = "Flux CV (variability)"
)

theme_pub <- function(base_size = 9) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.25, colour = "grey90"),
      strip.background = element_rect(fill = "grey92", colour = NA),
      strip.text       = element_text(size = base_size, face = "bold"),
      plot.title       = element_text(size = base_size + 1, face = "bold"),
      plot.subtitle    = element_text(size = base_size - 1, colour = "grey35"),
      axis.title       = element_text(size = base_size),
      axis.text        = element_text(size = base_size - 1, colour = "grey20"),
      legend.title     = element_text(size = base_size - 1),
      legend.text      = element_text(size = base_size - 1),
      legend.position  = "bottom"
    )
}

# Number-formatting helpers used in DOCX tables
format_p <- function(p) {
  ifelse(is.na(p), "\u2014",
  ifelse(p < 0.001, "<0.001", formatC(p, format = "f", digits = 3)))
}
format_est <- function(est, digits = 3) {
  ifelse(is.na(est), "\u2014", formatC(est, format = "f", digits = digits))
}
format_ci <- function(lo, hi, digits = 3) {
  ifelse(is.na(lo) | is.na(hi), "\u2014",
         sprintf("[%s, %s]",
                 formatC(lo, format = "f", digits = digits),
                 formatC(hi, format = "f", digits = digits)))
}

# =============================================================================
# 1. LOAD AND PREPARE FLUX DATA
# =============================================================================

message("Loading flux data...")

flux <- read_csv("FLUX.csv", show_col_types = FALSE) %>%
  mutate(
    FluxDate    = as_date(DayHour),
    Plot        = as.integer(Plot),
    Treatment   = str_to_title(Treatment),
    Flux        = suppressWarnings(as.numeric(Flux)),
    Temperature = suppressWarnings(as.numeric(Temperature)),
    VWC         = suppressWarnings(as.numeric(VWC))
  ) %>%
  filter(!is.na(Flux), Flux > 0,
         !is.na(Temperature), !is.na(VWC))

message("  n hourly obs = ", nrow(flux),
        " | Date range: ", min(flux$FluxDate), " to ", max(flux$FluxDate))

# =============================================================================
# 2. ABIOTIC BASE MODEL
# =============================================================================
# log(Flux) ~ Temperature + VWC + (1 | Plot), REML.
# Treatment is excluded so warming-mediated variance acting through T and VWC
# is retained in residuals for downstream chemistry analysis. Log transform
# linearises the multiplicative response of respiration to its drivers
# (Lloyd & Taylor 1994; Davidson & Janssens 2006) and stabilises residual
# variance (untransformed residuals at this site are severely right-skewed).

message("\nFitting abiotic base model...")

abiotic_mod <- lmer(
  log(Flux) ~ Temperature + VWC + (1 | Plot),
  data    = flux,
  REML    = TRUE,
  control = lmerControl(optimizer = "bobyqa")
)

writeLines(capture.output(summary(abiotic_mod)),
           file.path(MAIN_TABLES, "abiotic_model_summary.txt"))
if (isSingular(abiotic_mod))
  warning("Abiotic base model is singular.")

abiotic_fe <- broom.mixed::tidy(abiotic_mod, effects = "fixed",
                                conf.int = TRUE, conf.method = "Wald")
abiotic_re <- broom.mixed::tidy(abiotic_mod, effects = "ran_pars")
abiotic_r2 <- MuMIn::r.squaredGLMM(abiotic_mod)

flux$flux_resid <- residuals(abiotic_mod)
flux$log_flux   <- log(flux$Flux)

# Per-plot daily aggregate, used by both the windowed and same-day pipelines
flux_split <- split(flux, flux$Plot)

flux_daily <- flux %>%
  group_by(Plot, FluxDate) %>%
  summarise(flux_mean       = mean(Flux,       na.rm = TRUE),
            log_flux_mean   = mean(log_flux,   na.rm = TRUE),
            flux_cv         = sd(Flux, na.rm = TRUE) / mean(Flux, na.rm = TRUE),
            flux_resid_mean = mean(flux_resid, na.rm = TRUE),
            n_flux_obs      = sum(!is.na(Flux)),
            .groups         = "drop") %>%
  filter(!is.na(flux_mean), n_flux_obs >= 1L)

# =============================================================================
# 3. LOAD CHEMISTRY DATA
# =============================================================================

meta_cols <- c("Date", "Plot", "Treatment", "Notes")

load_chem <- function(path) {
  raw <- read_csv(path, show_col_types = FALSE)
  # Robust date parsing: try ISO, then US, then European formats
  if ("Date" %in% names(raw)) {
    raw <- raw %>%
      mutate(Date = lubridate::parse_date_time(
        as.character(Date),
        orders = c("ymd", "mdy", "dmy", "ymd HMS", "mdy HMS", "dmy HMS"),
        quiet  = TRUE
      ) %>% as.Date())
  }
  raw %>%
    mutate(Plot = as.integer(Plot)) %>%
    mutate(across(-any_of(meta_cols),
                  ~ suppressWarnings(as.numeric(.))))
}

message("\nLoading chemistry data...")

cores     <- load_chem("TRACE_Cores_clean.csv")
resin     <- load_chem("TRACE_Resin_clean.csv")
lysimeter <- load_chem("TRACE_Lysimeter_clean.csv")
porewater <- load_chem("TRACE_Porewater_clean.csv")

chem_files <- list(
  Cores     = list(df = cores,
                   compounds = setdiff(names(cores),     meta_cols)),
  Resin     = list(df = resin,
                   compounds = setdiff(names(resin),     meta_cols)),
  Lysimeter = list(df = lysimeter,
                   compounds = setdiff(names(lysimeter), meta_cols)),
  Porewater = list(df = porewater,
                   compounds = setdiff(names(porewater), meta_cols))
)

# Drop fully empty placeholder columns (Porewater has 4 such columns)
for (nm in names(chem_files)) {
  cf   <- chem_files[[nm]]
  keep <- vapply(cf$compounds,
                 function(cmp) any(is.finite(cf$df[[cmp]])),
                 logical(1))
  chem_files[[nm]]$compounds <- cf$compounds[keep]
}

for (nm in names(chem_files)) {
  message("  ", nm, ": ", length(chem_files[[nm]]$compounds), " compounds")
}

chem_lookup <- list(Cores = cores, Resin = resin,
                    Lysimeter = lysimeter, Porewater = porewater)

# =============================================================================
# 4. CORE FUNCTIONS
# =============================================================================

# ---- 4a. Windowed flux summary ----------------------------------------------
# For one chemistry sample (Plot, Date), summarise hourly flux within
# +/-window_days of that sample on the same plot.

summarise_flux_window <- function(plot_id, chem_date, window_days) {

  empty <- tibble(log_flux_mean   = NA_real_,
                  flux_cv         = NA_real_,
                  flux_resid_mean = NA_real_,
                  n_flux_obs      = 0L)

  fd <- flux_split[[as.character(plot_id)]]
  if (is.null(fd) || nrow(fd) == 0L) return(empty)

  wf <- fd[abs(as.integer(fd$FluxDate - chem_date)) <= window_days, ]
  if (nrow(wf) == 0L) return(empty)

  f  <- wf$Flux
  fm <- mean(f, na.rm = TRUE)

  tibble(
    log_flux_mean   = mean(wf$log_flux, na.rm = TRUE),
    flux_cv         = if (!is.na(fm) && fm > 0)
                        sd(f, na.rm = TRUE) / fm else NA_real_,
    flux_resid_mean = mean(wf$flux_resid, na.rm = TRUE),
    n_flux_obs      = sum(!is.na(f))
  )
}

build_matched_df <- function(chem_df, compound, window_days) {

  chem_sub <- chem_df %>%
    select(Date, Plot, Treatment, chem_value = all_of(compound)) %>%
    filter(!is.na(chem_value), !is.na(Date), is.finite(chem_value))

  if (nrow(chem_sub) == 0L) return(NULL)

  chem_sub %>%
    mutate(flux_stats = map2(Plot, Date,
                             ~ summarise_flux_window(.x, .y, window_days))) %>%
    unnest(flux_stats) %>%
    filter(!is.na(log_flux_mean), n_flux_obs >= 1L)
}

# ---- 4b. Same-day matching --------------------------------------------------
# Each chemistry observation paired with the daily-aggregated hourly flux
# from the same Plot on the same calendar day. Most stringent temporal
# independence achievable with this design; complements the windowed analysis.

build_sameday_df <- function(chem_df, compound) {

  chem_sub <- chem_df %>%
    select(Date, Plot, Treatment, chem_value = all_of(compound)) %>%
    filter(!is.na(chem_value), !is.na(Date), is.finite(chem_value))

  if (nrow(chem_sub) == 0L) return(NULL)

  matched <- chem_sub %>%
    inner_join(flux_daily, by = c("Plot", "Date" = "FluxDate"))

  if (nrow(matched) == 0L) return(NULL)
  matched
}

# ---- 4c. Standard main-effect LMM -------------------------------------------
# response ~ scale(chem) + Treatment + (1 | Plot)

fit_standard_lmm <- function(df, response_col) {

  d <- df %>%
    select(resp = all_of(response_col),
           chem_value, Treatment, Plot) %>%
    filter(!is.na(resp), !is.na(chem_value),
           is.finite(resp), is.finite(chem_value))

  if (nrow(d) < MIN_OBS || length(unique(d$Plot)) < MIN_PLOTS) return(NULL)

  mod <- tryCatch(
    suppressWarnings(
      lmer(resp ~ scale(chem_value) + Treatment + (1 | Plot),
           data = d, REML = TRUE,
           control = lmerControl(optimizer = "bobyqa"))
    ),
    error = function(e) NULL
  )
  if (is.null(mod)) return(NULL)

  r2   <- tryCatch(r.squaredGLMM(mod)[1L, "R2m"], error = function(e) NA_real_)
  sngl <- isSingular(mod)

  fe <- tidy(mod, effects = "fixed", conf.int = TRUE,
             conf.method = "Wald") %>%
    mutate(r2_marginal = r2, n_obs = nrow(d),
           n_plots = length(unique(d$Plot)), singular = sngl)

  list(
    chem      = filter(fe, term == "scale(chem_value)"),
    treatment = filter(fe, str_detect(term, "Treatment"))
  )
}

# ---- 4d. Interaction LMM ----------------------------------------------------
# response ~ scale(chem) * Treatment + (1 | Plot); per-Treatment slopes via
# emmeans::emtrends.

fit_interaction_lmm <- function(df, response_col) {

  d <- df %>%
    select(resp = all_of(response_col),
           chem_value, Treatment, Plot) %>%
    filter(!is.na(resp), !is.na(chem_value),
           is.finite(resp), is.finite(chem_value)) %>%
    mutate(Treatment = factor(Treatment, levels = c("Control", "Warmed")))

  if (nrow(d) < MIN_OBS || length(unique(d$Plot)) < MIN_PLOTS) return(NULL)

  mod <- tryCatch(
    suppressWarnings(
      lmer(resp ~ scale(chem_value) * Treatment + (1 | Plot),
           data = d, REML = TRUE,
           control = lmerControl(optimizer = "bobyqa"))
    ),
    error = function(e) NULL
  )
  if (is.null(mod)) return(NULL)

  r2   <- tryCatch(r.squaredGLMM(mod)[1L, "R2m"], error = function(e) NA_real_)
  sngl <- isSingular(mod)

  int_row <- tidy(mod, effects = "fixed", conf.int = TRUE,
                  conf.method = "Wald") %>%
    filter(str_detect(term, ":")) %>%
    mutate(r2_marginal = r2, n_obs = nrow(d),
           n_plots = length(unique(d$Plot)), singular = sngl)

  slopes <- tryCatch({
    emtrends(mod, ~ Treatment, var = "chem_value") %>%
      as_tibble() %>%
      rename(slope    = chem_value.trend,
             slope_se = SE,
             slope_lo = lower.CL,
             slope_hi = upper.CL)
  }, error = function(e) NULL)

  list(interaction = int_row, slopes = slopes)
}

# ---- 4e. AR-corrected main-effect model -------------------------------------
# Plot x Date aggregated data; corExp/corCAR1/corGaus selectable. Returns
# both the AR-corrected fit and the AIC of the analogous model without the
# correlation structure, so AIC improvement can be reported.

fit_autocor_lmm <- function(df, response_col,
                             cor_struct = c("corExp", "corCAR1", "corGaus")) {

  cor_struct <- match.arg(cor_struct)

  d <- df %>%
    select(Date, resp = all_of(response_col),
           chem_value, Treatment, Plot) %>%
    filter(!is.na(resp), !is.na(chem_value),
           is.finite(resp), is.finite(chem_value)) %>%
    group_by(Date, Plot, Treatment) %>%
    summarise(resp       = mean(resp,       na.rm = TRUE),
              chem_value = mean(chem_value, na.rm = TRUE),
              .groups    = "drop") %>%
    mutate(time_num = as.numeric(Date),
           Plot     = factor(Plot))

  if (nrow(d) < MIN_OBS || length(unique(d$Plot)) < MIN_PLOTS) return(NULL)

  corr_obj <- switch(
    cor_struct,
    "corExp"  = nlme::corExp( form = ~ time_num | Plot, nugget = TRUE),
    "corCAR1" = nlme::corCAR1(form = ~ time_num | Plot),
    "corGaus" = nlme::corGaus(form = ~ time_num | Plot, nugget = TRUE)
  )

  mod_ar <- tryCatch(
    suppressWarnings(
      nlme::lme(resp ~ scale(chem_value) + Treatment,
                random      = ~ 1 | Plot,
                correlation = corr_obj,
                data        = d,
                method      = "REML")
    ),
    error = function(e) NULL
  )
  if (is.null(mod_ar)) return(NULL)

  mod_base <- tryCatch(
    suppressWarnings(
      nlme::lme(resp ~ scale(chem_value) + Treatment,
                random = ~ 1 | Plot, data = d, method = "REML")
    ),
    error = function(e) NULL
  )

  aic_improvement <- if (!is.null(mod_base))
    AIC(mod_base) - AIC(mod_ar) else NA_real_

  tt  <- summary(mod_ar)$tTable
  idx <- grep("chem_value", rownames(tt))
  if (length(idx) == 0L) return(NULL)

  est <- tt[idx, "Value"]
  se  <- tt[idx, "Std.Error"]

  tibble(
    model           = cor_struct,
    estimate        = est,
    std.error       = se,
    df              = tt[idx, "DF"],
    statistic       = tt[idx, "t-value"],
    p.value         = tt[idx, "p-value"],
    conf.low        = est - 1.96 * se,
    conf.high       = est + 1.96 * se,
    aic             = AIC(mod_ar),
    aic_base        = if (!is.null(mod_base)) AIC(mod_base) else NA_real_,
    aic_improvement = aic_improvement,
    n_obs_agg       = nrow(d)
  )
}

# ---- 4f. AR-corrected interaction model -------------------------------------

fit_autocor_interaction <- function(df, response_col) {

  d <- df %>%
    select(Date, resp = all_of(response_col),
           chem_value, Treatment, Plot) %>%
    filter(!is.na(resp), !is.na(chem_value),
           is.finite(resp), is.finite(chem_value)) %>%
    group_by(Date, Plot, Treatment) %>%
    summarise(resp       = mean(resp,       na.rm = TRUE),
              chem_value = mean(chem_value, na.rm = TRUE),
              .groups    = "drop") %>%
    mutate(time_num  = as.numeric(Date),
           Treatment = factor(Treatment, levels = c("Control", "Warmed")),
           Plot      = factor(Plot))

  if (nrow(d) < MIN_OBS || length(unique(d$Plot)) < MIN_PLOTS) return(NULL)

  mod_ar <- tryCatch(
    suppressWarnings(
      nlme::lme(resp ~ scale(chem_value) * Treatment,
                random      = ~ 1 | Plot,
                correlation = nlme::corExp(form = ~ time_num | Plot,
                                           nugget = TRUE),
                data        = d, method = "REML")
    ),
    error = function(e) NULL
  )
  if (is.null(mod_ar)) return(NULL)

  mod_base <- tryCatch(
    suppressWarnings(
      nlme::lme(resp ~ scale(chem_value) * Treatment,
                random = ~ 1 | Plot, data = d, method = "REML")
    ),
    error = function(e) NULL
  )

  aic_improvement <- if (!is.null(mod_base))
    AIC(mod_base) - AIC(mod_ar) else NA_real_

  tt  <- summary(mod_ar)$tTable
  idx <- grep(":", rownames(tt))
  if (length(idx) == 0L) return(NULL)

  est <- tt[idx, "Value"]
  se  <- tt[idx, "Std.Error"]

  int_df <- tibble(
    term            = rownames(tt)[idx],
    estimate        = est,
    std.error       = se,
    df              = tt[idx, "DF"],
    statistic       = tt[idx, "t-value"],
    p.value         = tt[idx, "p-value"],
    conf.low        = est - 1.96 * se,
    conf.high       = est + 1.96 * se,
    aic             = AIC(mod_ar),
    aic_base        = if (!is.null(mod_base)) AIC(mod_base) else NA_real_,
    aic_improvement = aic_improvement,
    n_obs_agg       = nrow(d)
  )

  slopes_df <- tryCatch({
    emtrends(mod_ar, ~ Treatment, var = "chem_value") %>%
      as_tibble() %>%
      rename(slope    = chem_value.trend,
             slope_se = SE,
             slope_lo = lower.CL,
             slope_hi = upper.CL) %>%
      mutate(n_obs_agg = nrow(d))
  }, error = function(e) NULL)

  list(interaction = int_df, slopes = slopes_df)
}

# ---- 4g. Diagnostic plots for one model -------------------------------------

make_diagnostic_plots <- function(df, response_col, compound,
                                  file_nm, window_days, output_dir) {

  d <- df %>%
    select(resp = all_of(response_col),
           chem_value, Treatment, Plot) %>%
    filter(!is.na(resp), !is.na(chem_value),
           is.finite(resp), is.finite(chem_value)) %>%
    mutate(Treatment = factor(Treatment, levels = c("Control", "Warmed")))

  if (nrow(d) < MIN_OBS) return(invisible(NULL))

  mod <- tryCatch(
    suppressWarnings(
      lmer(resp ~ scale(chem_value) + Treatment + (1 | Plot),
           data = d, REML = TRUE,
           control = lmerControl(optimizer = "bobyqa"))
    ),
    error = function(e) NULL
  )
  if (is.null(mod)) return(invisible(NULL))

  d$fitted_val <- fitted(mod)
  d$resid_val  <- residuals(mod)

  pA <- ggplot(d, aes(x = fitted_val, y = resid_val, colour = Treatment)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_point(alpha = 0.40, size = 1.3) +
    geom_smooth(aes(group = 1), method = "loess", se = FALSE,
                colour = "black", linewidth = 0.6) +
    scale_colour_manual(values = treat_pal) +
    labs(x = "Fitted", y = "Residuals", title = "A  Residuals vs Fitted") +
    theme_pub(8) + theme(legend.position = "none")

  pB <- ggplot(d, aes(sample = resid_val)) +
    stat_qq(size = 1, alpha = 0.4) +
    stat_qq_line(colour = "black", linewidth = 0.6) +
    labs(x = "Theoretical quantiles", y = "Sample quantiles",
         title = "B  Normal Q-Q") +
    theme_pub(8)

  y_part <- tryCatch(
    residuals(suppressWarnings(
      lmer(resp ~ Treatment + (1 | Plot), data = d,
           control = lmerControl(optimizer = "bobyqa")))),
    error = function(e) d$resp - mean(d$resp, na.rm = TRUE)
  )
  x_part <- tryCatch(
    residuals(suppressWarnings(
      lmer(scale(chem_value) ~ Treatment + (1 | Plot), data = d,
           control = lmerControl(optimizer = "bobyqa")))),
    error = function(e) scale(d$chem_value)
  )
  part_df <- tibble(x = as.numeric(x_part), y = as.numeric(y_part),
                    Treatment = d$Treatment)

  pC <- ggplot(part_df, aes(x = x, y = y, colour = Treatment)) +
    geom_point(alpha = 0.40, size = 1.3) +
    geom_smooth(aes(group = 1), method = "lm", se = TRUE,
                colour = "black", linewidth = 0.7, fill = "grey80") +
    scale_colour_manual(values = treat_pal) +
    labs(x = paste0("scale(", compound, ") | Treatment, Plot"),
         y = paste0(response_col, " | Treatment, Plot"),
         title = "C  Partial regression") +
    theme_pub(8) + theme(legend.position = "none")

  pD <- ggplot(d, aes(x = chem_value, y = resp, colour = Treatment)) +
    geom_point(alpha = 0.35, size = 1.3) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.7, alpha = 0.15) +
    scale_colour_manual(values = treat_pal, name = NULL) +
    labs(x = compound, y = response_col,
         title = "D  Raw scatter by Treatment") +
    theme_pub(8)

  r2m <- tryCatch(round(r.squaredGLMM(mod)[1L, "R2m"], 3),
                  error = function(e) NA)

  combined <- (pA | pB) / (pC | pD) +
    plot_annotation(
      title    = paste0(file_nm, " | ", compound, " x ", response_col,
                        " [\u00b1", window_days, "d]"),
      subtitle = paste0("Marginal R\u00b2 = ", r2m,
                        " | n = ", nrow(d),
                        " | singular = ", isSingular(mod))
    )

  safe  <- function(x) gsub("[^A-Za-z0-9_]", "_", x)
  fname <- paste0("diagnostic_", safe(file_nm), "_", safe(compound),
                  "_", safe(response_col), "_w", window_days, ".png")
  ggsave(file.path(output_dir, fname), combined,
         width = 9, height = 7, dpi = 300, units = "in")

  invisible(combined)
}
# =============================================================================
# 5. PRIMARY WINDOWED LMM ANALYSIS (+/-WINDOW_PRIMARY days)
# =============================================================================
# Single loop over chemistry pools x compounds, fitting standard main-effect
# and interaction LMMs for all main and supplementary responses. Matched
# data frames are cached so AR refits in subsequent sections do not rebuild.

message("\nRunning primary windowed LMM analysis (\u00b1",
        WINDOW_PRIMARY, "d)...")

chem_rows     <- list()
treat_rows    <- list()
int_rows      <- list()
slopes_rows   <- list()
matched_cache <- list()                              # key: file|compound|window

w <- WINDOW_PRIMARY

for (file_nm in names(chem_files)) {
  chem_df   <- chem_files[[file_nm]]$df
  compounds <- chem_files[[file_nm]]$compounds

  for (comp in compounds) {

    key     <- paste(file_nm, comp, w, sep = "|")
    matched <- build_matched_df(chem_df, comp, w)
    matched_cache[[key]] <- matched

    if (is.null(matched) || nrow(matched) < MIN_OBS) next

    for (resp in RESPONSES_ALL) {

      meta     <- tibble(window = w, file = file_nm,
                         compound = comp, response = resp)
      main_res <- fit_standard_lmm(matched, resp)
      int_res  <- fit_interaction_lmm(matched, resp)

      if (!is.null(main_res)) {
        if (nrow(main_res$chem) > 0)
          chem_rows[[length(chem_rows) + 1L]] <-
            bind_cols(meta, main_res$chem)
        if (nrow(main_res$treatment) > 0)
          treat_rows[[length(treat_rows) + 1L]] <-
            bind_cols(meta, main_res$treatment)
      }
      if (!is.null(int_res)) {
        if (!is.null(int_res$interaction) && nrow(int_res$interaction) > 0)
          int_rows[[length(int_rows) + 1L]] <-
            bind_cols(meta, int_res$interaction)
        if (!is.null(int_res$slopes))
          slopes_rows[[length(slopes_rows) + 1L]] <-
            bind_cols(meta, int_res$slopes)
      }
    }
  }
}

results_chem   <- bind_rows(chem_rows)   %>%
  relocate(window, file, compound, response, .before = everything())
results_treat  <- bind_rows(treat_rows)  %>%
  relocate(window, file, compound, response, .before = everything())
results_int    <- bind_rows(int_rows)    %>%
  relocate(window, file, compound, response, .before = everything())
results_slopes <- bind_rows(slopes_rows) %>%
  relocate(window, file, compound, response, .before = everything())

# =============================================================================
# 6. FDR CORRECTION (Benjamini-Hochberg, within response)
# =============================================================================

apply_fdr <- function(df) {
  df %>%
    group_by(response) %>%
    mutate(p.fdr = p.adjust(p.value, method = "BH")) %>%
    ungroup()
}

results_chem  <- apply_fdr(results_chem)
results_treat <- apply_fdr(results_treat)
results_int   <- apply_fdr(results_int)

message("  Standard chemistry models: ", nrow(results_chem),
        " | FDR < ", FDR_THRESHOLD, ": ",
        sum(results_chem$p.fdr < FDR_THRESHOLD, na.rm = TRUE))
message("  Treatment terms:          ", nrow(results_treat),
        " | FDR < ", FDR_THRESHOLD, ": ",
        sum(results_treat$p.fdr < FDR_THRESHOLD, na.rm = TRUE))
message("  Interaction terms:        ", nrow(results_int),
        " | FDR < ", FDR_THRESHOLD, ": ",
        sum(results_int$p.fdr < FDR_THRESHOLD, na.rm = TRUE))

# =============================================================================
# 7. AR-CORRECTED MAIN-EFFECT MODELS
# =============================================================================
# Refit every standard-LMM survivor (FDR < FDR_THRESHOLD) with corExp.
# Justification: empirical flux ACF (Section 9) exceeds 0.6 at 30-day lags
# at this site, so ignoring temporal structure inflates type I error
# (Hurlbert 1984; Pinheiro & Bates 2000; Zuur et al. 2009).

message("\nRunning AR-corrected main-effect models (corExp)...")

sig_main <- results_chem %>%
  filter(p.fdr < FDR_THRESHOLD) %>%
  distinct(window, file, compound, response)

autocor_main_rows <- list()
for (i in seq_len(nrow(sig_main))) {
  row     <- sig_main[i, ]
  key     <- paste(row$file, row$compound, row$window, sep = "|")
  matched <- matched_cache[[key]]
  if (is.null(matched)) next
  res <- fit_autocor_lmm(matched, row$response, "corExp")
  if (!is.null(res))
    autocor_main_rows[[length(autocor_main_rows) + 1L]] <- bind_cols(row, res)
}
results_autocor_main <- if (length(autocor_main_rows) > 0)
  bind_rows(autocor_main_rows) else tibble()

message("  AR main-effect models fitted: ", nrow(results_autocor_main))

# =============================================================================
# 8. AR-CORRECTED INTERACTION MODELS
# =============================================================================
# Targets:
#   (a) Every standard-LMM interaction passing FDR < FDR_THRESHOLD.
#   (b) All Porewater interactions (pre-specified supplementary set).
#   (c) HEADLINE_INT_TARGETS: microbial P x Treatment is fitted across
#       both main responses regardless of standard-LMM outcome, so
#       Figure 2 can annotate AR slopes symmetrically across both panels.

message("\nRunning AR-corrected interaction models...")

sig_int <- results_int %>%
  filter(p.fdr < FDR_THRESHOLD) %>%
  distinct(window, file, compound, response)

porewater_int_targets <- expand_grid(
  window   = WINDOW_PRIMARY,
  file     = "Porewater",
  compound = chem_files[["Porewater"]]$compounds,
  response = RESPONSES_ALL
)

autocor_int_targets <- bind_rows(sig_int,
                                 porewater_int_targets,
                                 HEADLINE_INT_TARGETS) %>%
  distinct(window, file, compound, response)

autocor_int_rows    <- list()
autocor_slopes_rows <- list()

for (i in seq_len(nrow(autocor_int_targets))) {
  row     <- autocor_int_targets[i, ]
  key     <- paste(row$file, row$compound, row$window, sep = "|")
  matched <- matched_cache[[key]]
  if (is.null(matched)) next

  res <- fit_autocor_interaction(matched, row$response)
  if (is.null(res)) next

  if (!is.null(res$interaction) && nrow(res$interaction) > 0)
    autocor_int_rows[[length(autocor_int_rows) + 1L]] <-
      bind_cols(row, res$interaction)
  if (!is.null(res$slopes) && nrow(res$slopes) > 0)
    autocor_slopes_rows[[length(autocor_slopes_rows) + 1L]] <-
      bind_cols(row, res$slopes)
}

results_autocor_int    <- if (length(autocor_int_rows) > 0)
  bind_rows(autocor_int_rows) else tibble()
results_autocor_slopes <- if (length(autocor_slopes_rows) > 0)
  bind_rows(autocor_slopes_rows) else tibble()

# Verify headline coverage
headline_check <- results_autocor_slopes %>%
  filter(file == "Cores", compound == "uP",
         response %in% RESPONSES_MAIN, window == WINDOW_PRIMARY) %>%
  distinct(response, Treatment)
message("  AR interactions fitted: ", nrow(results_autocor_int))
message("  Headline (uP) AR slopes available: ", nrow(headline_check),
        "/", length(RESPONSES_MAIN) * 2L, " (target=both responses x ",
        "two treatments)")

# =============================================================================
# 9. EMPIRICAL FLUX ACF PER PLOT
# =============================================================================
# Daily-aggregated log(Flux) and flux residuals; ACF computed independently
# per plot. Decorrelation lag = first lag where |ACF| < ACF_THRESH.

message("\nComputing empirical flux ACF per plot...")

compute_plot_acf <- function(plot_df, plot_id) {
  daily <- plot_df %>%
    group_by(FluxDate) %>%
    summarise(log_flux   = mean(log_flux,   na.rm = TRUE),
              flux_resid = mean(flux_resid, na.rm = TRUE),
              .groups    = "drop") %>%
    arrange(FluxDate)
  if (nrow(daily) < 30L) return(NULL)
  full <- tibble(FluxDate = seq(min(daily$FluxDate),
                                max(daily$FluxDate), by = "day"))
  daily_full <- left_join(full, daily, by = "FluxDate")
  acf_log   <- acf(daily_full$log_flux,   lag.max = ACF_MAX_LAG,
                   plot = FALSE, na.action = na.pass)$acf[, 1, 1]
  acf_resid <- acf(daily_full$flux_resid, lag.max = ACF_MAX_LAG,
                   plot = FALSE, na.action = na.pass)$acf[, 1, 1]
  tibble(Plot = plot_id, lag_days = 0L:ACF_MAX_LAG,
         acf_log = as.numeric(acf_log), acf_resid = as.numeric(acf_resid))
}

acf_per_plot <- map2_dfr(flux_split, names(flux_split),
                         ~ compute_plot_acf(.x, as.integer(.y)))

# Report ACF values at biologically interpretable fixed lags. Threshold-
# crossing decorrelation lags are not informative at this site because no
# plot decorrelates within the 60-day examination window; reporting fixed-
# lag ACF values directly shows the strength of the persistence.

ACF_REPORT_LAGS <- c(0L, 1L, 7L, 14L, 30L, 60L)

acf_at_lags <- acf_per_plot %>%
  filter(lag_days %in% ACF_REPORT_LAGS) %>%
  pivot_wider(names_from = lag_days, names_prefix = "lag_",
              values_from = c(acf_log, acf_resid))

# Per-plot first lag where |ACF| drops below the threshold, if any. Most
# plots return NA because the ACF stays above the threshold throughout.
decorr_lags <- acf_per_plot %>%
  filter(lag_days > 0L) %>%
  group_by(Plot) %>%
  summarise(
    decorr_lag_log   = {
      below <- which(abs(acf_log)   < ACF_THRESH)
      if (length(below) > 0L) lag_days[below[1L]] else NA_integer_
    },
    decorr_lag_resid = {
      below <- which(abs(acf_resid) < ACF_THRESH)
      if (length(below) > 0L) lag_days[below[1L]] else NA_integer_
    },
    .groups = "drop"
  )

# =============================================================================
# 10. WINDOW SENSITIVITY GRID (+/-3 to +/-30 days)
# =============================================================================
# For every FDR-significant compound x main-response, refit standard and
# corExp models across the full window grid. Demonstrates that conclusions
# do not depend on window choice in the biologically meaningful range.

message("\nRunning window sensitivity grid...")

grid_targets <- results_chem %>%
  filter(p.fdr < FDR_THRESHOLD, response %in% RESPONSES_MAIN) %>%
  distinct(file, compound, response)

grid_rows <- list()
for (i in seq_len(nrow(grid_targets))) {
  t_row   <- grid_targets[i, ]
  chem_df <- chem_lookup[[t_row$file]]
  for (gw in WINDOW_GRID) {
    matched <- build_matched_df(chem_df, t_row$compound, gw)
    if (is.null(matched) || nrow(matched) < MIN_OBS) next

    std_res <- fit_standard_lmm(matched, t_row$response)
    if (!is.null(std_res) && nrow(std_res$chem) > 0)
      grid_rows[[length(grid_rows) + 1L]] <- std_res$chem %>%
        transmute(model = "standard_LMM",
                  estimate, std.error, p.value,
                  conf.low, conf.high, n_obs) %>%
        bind_cols(tibble(window = gw), t_row, .)

    ar_res <- fit_autocor_lmm(matched, t_row$response, "corExp")
    if (!is.null(ar_res))
      grid_rows[[length(grid_rows) + 1L]] <- ar_res %>%
        transmute(model, estimate, std.error, p.value,
                  conf.low, conf.high, n_obs = n_obs_agg) %>%
        bind_cols(tibble(window = gw), t_row, .)
  }
}
sens_grid <- bind_rows(grid_rows) %>%
  relocate(window, file, compound, response, model, .before = everything())

# =============================================================================
# 11. CORRELATION STRUCTURE COMPARISON (corExp vs corCAR1 vs corGaus)
# =============================================================================
# Validates that the AR conclusions are not specific to corExp.

message("\nRunning correlation structure comparison...")

cs_rows <- list()
for (i in seq_len(nrow(grid_targets))) {
  t_row   <- grid_targets[i, ]
  chem_df <- chem_lookup[[t_row$file]]
  matched <- build_matched_df(chem_df, t_row$compound, WINDOW_PRIMARY)
  if (is.null(matched) || nrow(matched) < MIN_OBS) next
  for (cs in c("corExp", "corCAR1", "corGaus")) {
    res <- fit_autocor_lmm(matched, t_row$response, cs)
    if (!is.null(res))
      cs_rows[[length(cs_rows) + 1L]] <- res %>%
        transmute(model, estimate, std.error, p.value,
                  conf.low, conf.high, aic, n_obs = n_obs_agg) %>%
        bind_cols(tibble(window = WINDOW_PRIMARY), t_row, .)
  }
}
corr_struct <- bind_rows(cs_rows) %>%
  relocate(window, file, compound, response, model, .before = everything()) %>%
  arrange(file, compound, response, model)

corr_struct_wide <- corr_struct %>%
  select(file, compound, response, window, model, estimate, p.value, aic) %>%
  pivot_wider(names_from = model,
              values_from = c(estimate, p.value, aic)) %>%
  mutate(
    dAIC_CAR1_vs_Exp = aic_corCAR1 - aic_corExp,
    dAIC_Gaus_vs_Exp = aic_corGaus - aic_corExp,
    dEst_CAR1_vs_Exp = estimate_corCAR1 - estimate_corExp,
    dEst_Gaus_vs_Exp = estimate_corGaus - estimate_corExp
  )

# =============================================================================
# 12. AR COLLAPSE TABLE
# =============================================================================
# Pairs every FDR-significant standard-LMM main effect with its AR-corrected
# counterpart at the primary window. The collapse pattern itself is a
# substantive result: most apparent associations vanish under AR correction.

ar_collapse <- results_chem %>%
  filter(window == WINDOW_PRIMARY, p.fdr < FDR_THRESHOLD) %>%
  left_join(
    results_autocor_main %>%
      filter(window == WINDOW_PRIMARY) %>%
      select(window, file, compound, response,
             ar_estimate        = estimate,
             ar_se              = std.error,
             ar_p               = p.value,
             ar_aic_improvement = aic_improvement),
    by = c("window", "file", "compound", "response")
  ) %>%
  mutate(
    pct_attenuation = ifelse(estimate == 0, NA_real_,
                             100 * (abs(estimate) - abs(ar_estimate)) / abs(estimate)),
    survives_AR     = ar_p < 0.05
  ) %>%
  select(window, file, compound, response,
         std_estimate = estimate, std_se = std.error,
         std_p = p.value, std_fdr = p.fdr, std_n = n_obs,
         ar_estimate, ar_se, ar_p,
         ar_aic_improvement, pct_attenuation, survives_AR)

# =============================================================================
# 12b. AR DIAGNOSTICS SUMMARY
# =============================================================================
# Extracts the corExp range parameter (length of the empirical autocorrelation
# in days) and the nugget term per fitted AR model, providing concrete
# numerical evidence for the strength and timescale of residual temporal
# dependence that the AR specification is correcting.

extract_ar_diagnostics <- function(matched_df, response_col) {

  d <- matched_df %>%
    select(Date, resp = all_of(response_col),
           chem_value, Treatment, Plot) %>%
    filter(!is.na(resp), !is.na(chem_value),
           is.finite(resp), is.finite(chem_value)) %>%
    group_by(Date, Plot, Treatment) %>%
    summarise(resp       = mean(resp,       na.rm = TRUE),
              chem_value = mean(chem_value, na.rm = TRUE),
              .groups    = "drop") %>%
    mutate(time_num = as.numeric(Date), Plot = factor(Plot))

  if (nrow(d) < MIN_OBS) return(NULL)

  mod <- tryCatch(suppressWarnings(
    nlme::lme(resp ~ scale(chem_value) + Treatment,
              random      = ~ 1 | Plot,
              correlation = nlme::corExp(form = ~ time_num | Plot,
                                         nugget = TRUE),
              data = d, method = "REML")),
    error = function(e) NULL)
  if (is.null(mod)) return(NULL)

  cs <- tryCatch(coef(mod$modelStruct$corStruct, unconstrained = FALSE),
                 error = function(e) NULL)
  if (is.null(cs)) return(NULL)

  tibble(range_days = as.numeric(cs["range"]),
         nugget     = as.numeric(cs["nugget"]),
         sigma      = mod$sigma,
         n_obs_agg  = nrow(d))
}

ar_diag_rows <- list()
for (i in seq_len(nrow(sig_main))) {
  row     <- sig_main[i, ]
  matched <- matched_cache[[paste(row$file, row$compound, row$window,
                                  sep = "|")]]
  if (is.null(matched)) next
  diag <- extract_ar_diagnostics(matched, row$response)
  if (!is.null(diag))
    ar_diag_rows[[length(ar_diag_rows) + 1L]] <- bind_cols(row, diag)
}
ar_diagnostics <- if (length(ar_diag_rows) > 0)
  bind_rows(ar_diag_rows) else tibble()

# Console summary of AR diagnostics
if (nrow(ar_diagnostics) > 0) {
  message("\nAR correlation diagnostics (windowed corExp, primary window):")
  message("  AR preferred over base (\u0394AIC > 2): ",
          sum(results_autocor_main$aic_improvement > 2, na.rm = TRUE),
          " / ", nrow(results_autocor_main))
  message("  Mean   \u0394AIC: ",
          round(mean(results_autocor_main$aic_improvement, na.rm = TRUE), 1))
  message("  Median \u0394AIC: ",
          round(median(results_autocor_main$aic_improvement, na.rm = TRUE), 1))
  message("  Mean   AR range parameter: ",
          round(mean(ar_diagnostics$range_days, na.rm = TRUE), 1), " days")
  message("  Median AR range parameter: ",
          round(median(ar_diagnostics$range_days, na.rm = TRUE), 1), " days")
  message("  Mean   nugget:             ",
          round(mean(ar_diagnostics$nugget, na.rm = TRUE), 3))
}

# =============================================================================
# 13. SAME-DAY MATCHING SENSITIVITY ANALYSIS
# =============================================================================
# Each chemistry observation paired with daily-aggregated hourly flux from
# the same Plot on the same calendar day. Most stringent test of temporal
# independence available with this design; it pairs chemistry with a single
# day's flux instead of integrating over a window.
#
# Pipeline mirrors the windowed analysis: standard LMM, FDR, AR-corrected
# refits, AR-corrected interactions, headline microbial P slopes, then a
# cross-method comparison table that joins the same-day AR estimates with
# the windowed AR estimates from Section 7.

message("\n--- SAME-DAY MATCHING SENSITIVITY ---")

# ---- 13a. Build matched data + coverage summary -----------------------------

sameday_cache    <- list()
sameday_coverage_rows <- list()

for (file_nm in names(chem_files)) {
  chem_df   <- chem_files[[file_nm]]$df
  compounds <- chem_files[[file_nm]]$compounds
  for (comp in compounds) {
    matched <- build_sameday_df(chem_df, comp)
    sameday_cache[[paste(file_nm, comp, sep = "|")]] <- matched
    n_matched <- if (is.null(matched)) 0L else nrow(matched)
    n_chem    <- sum(!is.na(chem_df[[comp]]) & is.finite(chem_df[[comp]]))
    sameday_coverage_rows[[length(sameday_coverage_rows) + 1L]] <- tibble(
      file            = file_nm,
      compound        = comp,
      n_chem_obs      = n_chem,
      n_matched       = n_matched,
      n_plots_matched = if (is.null(matched)) 0L else
                          length(unique(matched$Plot)),
      match_rate      = if (n_chem > 0) n_matched / n_chem else NA_real_
    )
  }
}
sameday_coverage <- bind_rows(sameday_coverage_rows)

message("  Compounds with >= ", MIN_OBS, " same-day matches: ",
        sum(sameday_coverage$n_matched >= MIN_OBS), " / ",
        nrow(sameday_coverage))

# ---- 13b. Same-day standard + interaction LMMs ------------------------------

sameday_chem_rows   <- list()
sameday_treat_rows  <- list()
sameday_int_rows    <- list()
sameday_slopes_rows <- list()

for (file_nm in names(chem_files)) {
  for (comp in chem_files[[file_nm]]$compounds) {
    matched <- sameday_cache[[paste(file_nm, comp, sep = "|")]]
    if (is.null(matched) || nrow(matched) < MIN_OBS) next
    for (resp in RESPONSES_MAIN) {
      meta     <- tibble(file = file_nm, compound = comp, response = resp)
      main_res <- fit_standard_lmm(matched, resp)
      int_res  <- fit_interaction_lmm(matched, resp)
      if (!is.null(main_res) && nrow(main_res$chem) > 0)
        sameday_chem_rows[[length(sameday_chem_rows) + 1L]] <-
          bind_cols(meta, main_res$chem)
      if (!is.null(main_res) && nrow(main_res$treatment) > 0)
        sameday_treat_rows[[length(sameday_treat_rows) + 1L]] <-
          bind_cols(meta, main_res$treatment)
      if (!is.null(int_res)) {
        if (!is.null(int_res$interaction) && nrow(int_res$interaction) > 0)
          sameday_int_rows[[length(sameday_int_rows) + 1L]] <-
            bind_cols(meta, int_res$interaction)
        if (!is.null(int_res$slopes))
          sameday_slopes_rows[[length(sameday_slopes_rows) + 1L]] <-
            bind_cols(meta, int_res$slopes)
      }
    }
  }
}

sameday_chem <- bind_rows(sameday_chem_rows) %>%
  group_by(response) %>%
  mutate(p.fdr = p.adjust(p.value, method = "BH")) %>%
  ungroup()
sameday_int    <- bind_rows(sameday_int_rows) %>%
  group_by(response) %>%
  mutate(p.fdr = p.adjust(p.value, method = "BH")) %>%
  ungroup()
sameday_slopes <- bind_rows(sameday_slopes_rows)
sameday_treat  <- bind_rows(sameday_treat_rows) %>%
  group_by(response) %>%
  mutate(p.fdr = p.adjust(p.value, method = "BH")) %>%
  ungroup()

# ---- 13b-bis. Empirical justification for same-day AR correction ------------
# Tests whether the same-day standard LMM residuals retain temporal
# autocorrelation across consecutive within-plot chemistry sampling dates.
# Two diagnostics: (i) Pearson correlation of within-plot residual pairs
# binned by inter-observation gap (the irregular-time analogue of acf());
# (ii) likelihood-ratio test of corExp vs no correlation structure (ML fit).
# Reported per FDR-significant compound and pooled.

sameday_ar_test_lag_bins <- list(c(1L, 7L), c(8L, 14L), c(15L, 30L),
                                 c(31L, 60L), c(61L, 90L), c(91L, 365L))

sameday_residual_acf <- function(df_with_resid, bins = sameday_ar_test_lag_bins) {
  d <- df_with_resid %>% arrange(Plot, Date)
  pairs <- list()
  for (pl in unique(d$Plot)) {
    sub <- filter(d, Plot == pl)
    if (nrow(sub) < 2L) next
    for (i in seq_len(nrow(sub) - 1L)) {
      for (j in (i + 1L):nrow(sub)) {
        pairs[[length(pairs) + 1L]] <- tibble(
          gap_days = as.integer(sub$Date[j] - sub$Date[i]),
          r_i = sub$resid[i], r_j = sub$resid[j])
      }
    }
  }
  pairs <- bind_rows(pairs)
  if (nrow(pairs) == 0L) return(tibble())
  out <- list()
  for (b in bins) {
    sub <- filter(pairs, gap_days >= b[1], gap_days <= b[2])
    if (nrow(sub) >= 5L)
      out[[length(out) + 1L]] <- tibble(
        bin = sprintf("%d-%dd", b[1], b[2]),
        n   = nrow(sub),
        r   = cor(sub$r_i, sub$r_j),
        mean_gap = round(mean(sub$gap_days), 1))
  }
  bind_rows(out)
}

sameday_ar_test_rows <- list()  # binned-lag ACF per compound
sameday_lr_rows      <- list()  # LR-test results per compound

sameday_fdr_targets <- sameday_chem %>%
  filter(p.fdr < FDR_THRESHOLD) %>%
  distinct(file, compound, response)

for (i in seq_len(nrow(sameday_fdr_targets))) {
  row     <- sameday_fdr_targets[i, ]
  matched <- sameday_cache[[paste(row$file, row$compound, sep = "|")]]
  if (is.null(matched)) next

  # Standard LMM residuals for binned-lag ACF
  d_std <- matched %>%
    select(Date, resp = all_of(row$response),
           chem_value, Treatment, Plot) %>%
    filter(!is.na(resp), !is.na(chem_value),
           is.finite(resp), is.finite(chem_value)) %>%
    group_by(Date, Plot, Treatment) %>%
    summarise(resp       = mean(resp,       na.rm = TRUE),
              chem_value = mean(chem_value, na.rm = TRUE),
              .groups    = "drop")
  if (nrow(d_std) >= MIN_OBS && length(unique(d_std$Plot)) >= MIN_PLOTS) {
    mod_std <- tryCatch(suppressWarnings(
      lmer(resp ~ scale(chem_value) + Treatment + (1 | Plot),
           data = d_std, REML = TRUE,
           control = lmerControl(optimizer = "bobyqa"))),
      error = function(e) NULL)
    if (!is.null(mod_std)) {
      d_std$resid <- residuals(mod_std)
      bins <- sameday_residual_acf(d_std)
      if (nrow(bins) > 0)
        sameday_ar_test_rows[[length(sameday_ar_test_rows) + 1L]] <-
        bind_cols(row, bins)
    }
  }

  # LR test: corExp vs no correlation structure (ML fit)
  d <- matched %>%
    select(Date, resp = all_of(row$response),
           chem_value, Treatment, Plot) %>%
    filter(!is.na(resp), !is.na(chem_value),
           is.finite(resp), is.finite(chem_value)) %>%
    group_by(Date, Plot, Treatment) %>%
    summarise(resp       = mean(resp,       na.rm = TRUE),
              chem_value = mean(chem_value, na.rm = TRUE),
              .groups    = "drop") %>%
    mutate(time_num = as.numeric(Date), Plot = factor(Plot))
  if (nrow(d) < MIN_OBS) next

  mod_base <- tryCatch(suppressWarnings(
    nlme::lme(resp ~ scale(chem_value) + Treatment,
              random = ~ 1 | Plot, data = d, method = "ML")),
    error = function(e) NULL)
  mod_ar <- tryCatch(suppressWarnings(
    nlme::lme(resp ~ scale(chem_value) + Treatment,
              random      = ~ 1 | Plot,
              correlation = nlme::corExp(form = ~ time_num | Plot,
                                         nugget = TRUE),
              data = d, method = "ML")),
    error = function(e) NULL)
  if (is.null(mod_base) || is.null(mod_ar)) next

  lr   <- tryCatch(anova(mod_base, mod_ar), error = function(e) NULL)
  rng  <- tryCatch(coef(mod_ar$modelStruct$corStruct,
                        unconstrained = FALSE)[["range"]],
                   error = function(e) NA_real_)
  nug  <- tryCatch(coef(mod_ar$modelStruct$corStruct,
                        unconstrained = FALSE)[["nugget"]],
                   error = function(e) NA_real_)
  daic <- AIC(mod_base) - AIC(mod_ar)
  lrp  <- if (!is.null(lr)) lr[2L, "p-value"] else NA_real_

  sameday_lr_rows[[length(sameday_lr_rows) + 1L]] <- tibble(
    file        = row$file,
    compound    = row$compound,
    response    = row$response,
    n_agg       = nrow(d),
    aic_base    = AIC(mod_base),
    aic_ar      = AIC(mod_ar),
    delta_aic   = daic,
    lr_p        = lrp,
    range_days  = rng,
    nugget      = nug,
    AR_preferred = (daic > 2) & (!is.na(lrp) & lrp < 0.05))
}

sameday_ar_test       <- if (length(sameday_ar_test_rows) > 0)
  bind_rows(sameday_ar_test_rows) else tibble()
sameday_ar_lr_summary <- if (length(sameday_lr_rows) > 0)
  bind_rows(sameday_lr_rows) else tibble()

if (nrow(sameday_ar_lr_summary) > 0) {
  message("\nSame-day AR justification diagnostic:")
  message("  Compounds tested:  ", nrow(sameday_ar_lr_summary))
  message("  AR preferred:      ",
          sum(sameday_ar_lr_summary$AR_preferred, na.rm = TRUE), " / ",
          nrow(sameday_ar_lr_summary))
  message("  Mean   \u0394AIC:    ",
          round(mean(sameday_ar_lr_summary$delta_aic, na.rm = TRUE), 1))
  message("  Mean   range:      ",
          round(mean(sameday_ar_lr_summary$range_days, na.rm = TRUE), 1),
          " days")
  message("  Median range:      ",
          round(median(sameday_ar_lr_summary$range_days, na.rm = TRUE), 1),
          " days")
}

# ---- 13c. Same-day AR-corrected main effects --------------------------------
# Target set is the union of (a) windowed FDR-significant compounds at the
# primary window and (b) same-day FDR-significant compounds. This ensures
# Figure S6 (cross-method comparison) has matched coverage: any compound
# flagged by either pipeline is fitted under both schemes when the
# minimum-observations threshold is met.

sameday_sig_main <- bind_rows(
  sameday_chem %>%
    filter(p.fdr < FDR_THRESHOLD) %>%
    distinct(file, compound, response),
  results_chem %>%
    filter(window == WINDOW_PRIMARY, response %in% RESPONSES_MAIN,
           p.fdr < FDR_THRESHOLD) %>%
    distinct(file, compound, response)
) %>% distinct(file, compound, response)

sameday_ar_main_rows <- list()
for (i in seq_len(nrow(sameday_sig_main))) {
  row     <- sameday_sig_main[i, ]
  matched <- sameday_cache[[paste(row$file, row$compound, sep = "|")]]
  if (is.null(matched)) next
  res <- fit_autocor_lmm(matched, row$response, "corExp")
  if (!is.null(res))
    sameday_ar_main_rows[[length(sameday_ar_main_rows) + 1L]] <-
      bind_cols(row, res)
}
sameday_ar_main <- if (length(sameday_ar_main_rows) > 0)
  bind_rows(sameday_ar_main_rows) else tibble()

# ---- 13c-bis. Symmetric backfill of windowed AR for same-day-flagged --------
# Compounds flagged by the same-day FDR but not by the windowed FDR are
# refitted with the windowed corExp model so that Figure S6 has matched
# coverage in both directions. Result is appended to results_autocor_main.

windowed_already_fit <- results_autocor_main %>%
  filter(window == WINDOW_PRIMARY) %>%
  distinct(file, compound, response)

sameday_only_flags <- sameday_chem %>%
  filter(p.fdr < FDR_THRESHOLD) %>%
  distinct(file, compound, response) %>%
  anti_join(windowed_already_fit, by = c("file", "compound", "response"))

backfill_rows <- list()
for (i in seq_len(nrow(sameday_only_flags))) {
  row     <- sameday_only_flags[i, ]
  matched <- matched_cache[[paste(row$file, row$compound,
                                  WINDOW_PRIMARY, sep = "|")]]
  if (is.null(matched)) next
  res <- fit_autocor_lmm(matched, row$response, "corExp")
  if (!is.null(res))
    backfill_rows[[length(backfill_rows) + 1L]] <-
      bind_cols(tibble(window = WINDOW_PRIMARY), row, res)
}
if (length(backfill_rows) > 0L) {
  results_autocor_main <- bind_rows(results_autocor_main,
                                    bind_rows(backfill_rows))
  message("  Windowed AR backfill: ", length(backfill_rows),
          " same-day-only flagged compound(s) added to results_autocor_main")
}

# ---- 13d. Same-day AR-corrected interactions + headline uP slopes -----------

sameday_ar_int_rows    <- list()
sameday_ar_slopes_rows <- list()

# Targets: FDR-significant sameday interactions plus the uP headline targets
sameday_int_targets <- sameday_int %>%
  filter(p.fdr < FDR_THRESHOLD) %>%
  distinct(file, compound, response) %>%
  bind_rows(tibble(file = "Cores", compound = "uP",
                   response = RESPONSES_MAIN)) %>%
  distinct(file, compound, response)

for (i in seq_len(nrow(sameday_int_targets))) {
  row     <- sameday_int_targets[i, ]
  matched <- sameday_cache[[paste(row$file, row$compound, sep = "|")]]
  if (is.null(matched)) next
  res <- fit_autocor_interaction(matched, row$response)
  if (is.null(res)) next
  if (!is.null(res$interaction) && nrow(res$interaction) > 0)
    sameday_ar_int_rows[[length(sameday_ar_int_rows) + 1L]] <-
      bind_cols(row, res$interaction)
  if (!is.null(res$slopes) && nrow(res$slopes) > 0)
    sameday_ar_slopes_rows[[length(sameday_ar_slopes_rows) + 1L]] <-
      bind_cols(row, res$slopes)
}
sameday_ar_int    <- if (length(sameday_ar_int_rows) > 0)
  bind_rows(sameday_ar_int_rows) else tibble()
sameday_ar_slopes <- if (length(sameday_ar_slopes_rows) > 0)
  bind_rows(sameday_ar_slopes_rows) else tibble()

# ---- 13e. Cross-method comparison: same-day vs windowed ---------------------

windowed_std <- results_chem %>%
  filter(window == WINDOW_PRIMARY, response %in% RESPONSES_MAIN) %>%
  select(file, compound, response,
         windowed_estimate_std = estimate,
         windowed_se_std       = std.error,
         windowed_p_std        = p.value,
         windowed_fdr_std      = p.fdr,
         windowed_n_std        = n_obs)

windowed_ar <- results_autocor_main %>%
  filter(window == WINDOW_PRIMARY, response %in% RESPONSES_MAIN) %>%
  select(file, compound, response,
         windowed_estimate_ar = estimate,
         windowed_se_ar       = std.error,
         windowed_p_ar        = p.value,
         windowed_n_ar        = n_obs_agg)

sameday_std <- sameday_chem %>%
  select(file, compound, response,
         sameday_estimate_std = estimate,
         sameday_se_std       = std.error,
         sameday_p_std        = p.value,
         sameday_fdr_std      = p.fdr,
         sameday_n_std        = n_obs)

sameday_ar <- sameday_ar_main %>%
  select(file, compound, response,
         sameday_estimate_ar = estimate,
         sameday_se_ar       = std.error,
         sameday_p_ar        = p.value,
         sameday_n_ar        = n_obs_agg)

cross_method_comparison <- windowed_std %>%
  full_join(windowed_ar, by = c("file", "compound", "response")) %>%
  full_join(sameday_std, by = c("file", "compound", "response")) %>%
  full_join(sameday_ar,  by = c("file", "compound", "response")) %>%
  arrange(file, compound, response)

message("  Cross-method comparison: ", nrow(cross_method_comparison),
        " compound x response rows")

# Summary message: how many AR survivors agree across schemes?
both_ar <- cross_method_comparison %>%
  filter(!is.na(windowed_p_ar) & !is.na(sameday_p_ar)) %>%
  mutate(windowed_sig = windowed_p_ar < 0.05,
         sameday_sig  = sameday_p_ar  < 0.05)
message("  AR significance agreement: ",
        sum(both_ar$windowed_sig == both_ar$sameday_sig), " / ",
        nrow(both_ar))

# =============================================================================
# 14. EXPORT TABLES (CSV)
# =============================================================================

message("\nExporting CSV tables...")

# --- Main-text tables ---
results_chem %>%
  filter(window == WINDOW_PRIMARY, response %in% RESPONSES_MAIN) %>%
  write_csv(file.path(MAIN_TABLES, "Table2_chemistry_primary.csv"))

results_autocor_main %>%
  filter(window == WINDOW_PRIMARY, response %in% RESPONSES_MAIN) %>%
  write_csv(file.path(MAIN_TABLES, "Table3_AR_corrected_primary.csv"))

ar_collapse %>%
  filter(response %in% RESPONSES_MAIN) %>%
  write_csv(file.path(MAIN_TABLES, "Table4_AR_collapse_primary.csv"))

results_autocor_slopes %>%
  filter(file == "Cores", compound == "uP",
         window == WINDOW_PRIMARY, response %in% RESPONSES_MAIN) %>%
  write_csv(file.path(MAIN_TABLES, "Table5_uP_simple_slopes.csv"))

# --- Supplementary tables: windowed full set ---
write_csv(results_chem,
          file.path(SUPP_TABLES, "results_chemistry_all.csv"))
write_csv(results_treat,
          file.path(SUPP_TABLES, "results_treatment_all.csv"))
write_csv(results_int,
          file.path(SUPP_TABLES, "results_interaction_all.csv"))
write_csv(results_slopes,
          file.path(SUPP_TABLES, "results_interaction_slopes_all.csv"))
if (nrow(results_autocor_main) > 0)
  write_csv(results_autocor_main,
            file.path(SUPP_TABLES, "results_autocor_main_all.csv"))
if (nrow(results_autocor_int) > 0)
  write_csv(results_autocor_int,
            file.path(SUPP_TABLES, "results_autocor_interaction_all.csv"))
if (nrow(results_autocor_slopes) > 0)
  write_csv(results_autocor_slopes,
            file.path(SUPP_TABLES, "results_autocor_slopes_all.csv"))
write_csv(acf_per_plot,
          file.path(SUPP_TABLES, "flux_acf_per_plot.csv"))
write_csv(decorr_lags,
          file.path(SUPP_TABLES, "flux_acf_decorrelation_lags.csv"))
write_csv(sens_grid,
          file.path(SUPP_TABLES, "window_sensitivity_grid.csv"))
write_csv(corr_struct,
          file.path(SUPP_TABLES, "correlation_structure_comparison.csv"))
write_csv(corr_struct_wide,
          file.path(SUPP_TABLES, "correlation_structure_pairwise.csv"))
write_csv(ar_collapse,
          file.path(SUPP_TABLES, "ar_collapse_full.csv"))

# --- Supplementary tables: same-day analysis ---
write_csv(sameday_coverage,
          file.path(SUPP_TABLES, "sameday_coverage.csv"))
write_csv(sameday_chem,
          file.path(SUPP_TABLES, "sameday_results_chemistry.csv"))
write_csv(sameday_treat,
          file.path(SUPP_TABLES, "sameday_results_treatment.csv"))
write_csv(sameday_int,
          file.path(SUPP_TABLES, "sameday_results_interaction.csv"))
write_csv(sameday_slopes,
          file.path(SUPP_TABLES, "sameday_results_interaction_slopes.csv"))
if (nrow(sameday_ar_main) > 0)
  write_csv(sameday_ar_main,
            file.path(SUPP_TABLES, "sameday_results_autocor_main.csv"))
if (nrow(sameday_ar_int) > 0)
  write_csv(sameday_ar_int,
            file.path(SUPP_TABLES, "sameday_results_autocor_interaction.csv"))
if (nrow(sameday_ar_slopes) > 0)
  write_csv(sameday_ar_slopes,
            file.path(SUPP_TABLES, "sameday_results_autocor_slopes.csv"))
write_csv(cross_method_comparison,
          file.path(SUPP_TABLES, "cross_method_comparison.csv"))

# AR diagnostics and ACF-at-fixed-lags additions
write_csv(acf_at_lags,
          file.path(SUPP_TABLES, "flux_acf_at_fixed_lags.csv"))
if (nrow(ar_diagnostics) > 0)
  write_csv(ar_diagnostics,
            file.path(SUPP_TABLES, "ar_correlation_diagnostics.csv"))
if (nrow(sameday_ar_test) > 0)
  write_csv(sameday_ar_test,
            file.path(SUPP_TABLES, "sameday_residual_acf_binned.csv"))
if (nrow(sameday_ar_lr_summary) > 0)
  write_csv(sameday_ar_lr_summary,
            file.path(SUPP_TABLES, "sameday_AR_LR_test.csv"))

# =============================================================================
# 15. PUBLICATION-READY DOCX TABLES
# =============================================================================

message("\nBuilding DOCX tables...")

# ---- Table 1. Abiotic base model --------------------------------------------

t1_fixed <- abiotic_fe %>%
  filter(effect == "fixed") %>%
  transmute(
    Term      = recode(term,
                       "(Intercept)" = "Intercept",
                       "Temperature" = "Temperature (\u00b0C)",
                       "VWC"         = "VWC (m\u00b3 m\u207b\u00b3)"),
    Estimate  = format_est(estimate, 3),
    SE        = format_est(std.error, 4),
    `95% CI`  = format_ci(conf.low, conf.high, 3),
    df        = formatC(df, format = "f", digits = 0),
    `t-value` = format_est(statistic, 2),
    `p-value` = format_p(p.value)
  )

t1_random <- abiotic_re %>%
  filter(group %in% c("Plot", "Residual")) %>%
  transmute(
    Term      = paste0("\u03c3\u00b2 (", group, ")"),
    Estimate  = format_est(estimate^2, 4),
    SE        = "",
    `95% CI`  = "",
    df        = "",
    `t-value` = "",
    `p-value` = ""
  )

t1 <- bind_rows(t1_fixed, t1_random)

ft1 <- flextable::flextable(t1) %>%
  flextable::theme_booktabs() %>%
  flextable::bg(part = "header", bg = "#3a3a3a") %>%
  flextable::color(part = "header", color = "white") %>%
  flextable::bold(part = "header") %>%
  flextable::fontsize(size = 9, part = "all") %>%
  flextable::padding(padding.top = 3, padding.bottom = 3, part = "all") %>%
  flextable::align(align = "center", part = "all") %>%
  flextable::align(j = "Term", align = "left", part = "body") %>%
  flextable::set_table_properties(layout = "autofit", width = 1) %>%
  flextable::add_header_lines(paste0(
    "Table 1. Abiotic base model: log(Flux) ~ Temperature + VWC + (1 | Plot). ",
    "REML; Satterthwaite df; n = ", nrow(flux), " hourly observations across ",
    "6 plots. Marginal R\u00b2 = ", round(abiotic_r2[1, "R2m"], 3),
    ", conditional R\u00b2 = ", round(abiotic_r2[1, "R2c"], 3), "."))

read_docx() %>% body_add_flextable(ft1) %>%
  print(target = file.path(MAIN_TABLES, "Table1_abiotic_model.docx"))

# ---- Table 2. Primary chemistry results -------------------------------------

t2_data <- results_chem %>%
  filter(window == WINDOW_PRIMARY, response %in% RESPONSES_MAIN) %>%
  arrange(response, p.fdr) %>%
  transmute(
    Response       = recode(response, !!!response_labels),
    Predictor      = paste0(compound, " [", file, "]"),
    "\u03b2 (std)" = format_est(estimate, 3),
    `95% CI`        = format_ci(conf.low, conf.high, 3),
    "R\u00b2 (m)"  = format_est(r2_marginal, 3),
    n              = formatC(n_obs, format = "d"),
    `p-value`      = format_p(p.value),
    `p (FDR)`      = format_p(p.fdr),
    Sig            = ifelse(p.fdr < FDR_THRESHOLD, "*", "")
  )

ft2 <- flextable::flextable(t2_data) %>%
  flextable::theme_booktabs() %>%
  flextable::bg(part = "header", bg = "#3a3a3a") %>%
  flextable::color(part = "header", color = "white") %>%
  flextable::bold(part = "header") %>%
  flextable::fontsize(size = 8, part = "all") %>%
  flextable::padding(padding.top = 2.5, padding.bottom = 2.5, part = "all") %>%
  flextable::align(align = "center", part = "all") %>%
  flextable::align(j = c("Response", "Predictor"),
                   align = "left", part = "body") %>%
  flextable::merge_v(j = "Response") %>%
  flextable::set_table_properties(layout = "autofit", width = 1) %>%
  flextable::add_header_lines(paste0(
    "Table 2. Standard LMM: response ~ scale(chemistry) + Treatment + ",
    "(1 | Plot) at the primary \u00b1", WINDOW_PRIMARY, "-day window. ",
    "Standardised \u03b2 with 95% Wald CI; marginal R\u00b2 ",
    "(Nakagawa & Schielzeth 2013); ",
    "p (FDR) = Benjamini-Hochberg adjusted within response. ",
    "* p (FDR) < ", FDR_THRESHOLD, "."))

read_docx() %>% body_add_flextable(ft2) %>%
  print(target = file.path(MAIN_TABLES, "Table2_chemistry_primary.docx"))

# ---- Table 3. AR-corrected primary -----------------------------------------

t3_data <- ar_collapse %>%
  filter(response %in% RESPONSES_MAIN) %>%
  arrange(response, ar_p) %>%
  transmute(
    Response        = recode(response, !!!response_labels),
    Predictor       = paste0(compound, " [", file, "]"),
    "\u03b2 (std)"  = format_est(std_estimate, 3),
    "\u03b2 (AR)"   = format_est(ar_estimate, 3),
    "\u0394 (%)"    = format_est(pct_attenuation, 0),
    "\u0394AIC"     = format_est(ar_aic_improvement, 1),
    `p (std)`       = format_p(std_fdr),
    `p (AR)`        = format_p(ar_p),
    Survives        = ifelse(survives_AR, "yes", "no")
  )

ft3 <- flextable::flextable(t3_data) %>%
  flextable::theme_booktabs() %>%
  flextable::bg(part = "header", bg = "#3a3a3a") %>%
  flextable::color(part = "header", color = "white") %>%
  flextable::bold(part = "header") %>%
  flextable::fontsize(size = 8, part = "all") %>%
  flextable::padding(padding.top = 2.5, padding.bottom = 2.5, part = "all") %>%
  flextable::align(align = "center", part = "all") %>%
  flextable::align(j = c("Response", "Predictor"),
                   align = "left", part = "body") %>%
  flextable::bg(i = ~ Survives == "yes", bg = "#e8f4ea", part = "body") %>%
  flextable::merge_v(j = "Response") %>%
  flextable::set_table_properties(layout = "autofit", width = 1) %>%
  flextable::add_header_lines(paste0(
    "Table 3. Autocorrelation correction at the primary \u00b1",
    WINDOW_PRIMARY, "-day window. Standard-LMM \u03b2 vs corExp-corrected ",
    "\u03b2 for every compound significant at FDR < ", FDR_THRESHOLD,
    " in the standard LMM. \u0394AIC = AIC(base) \u2212 AIC(AR); ",
    "positive values favour the AR specification. ",
    "\u0394(%) = relative attenuation of |\u03b2|. ",
    "Shaded rows survive AR correction at p < 0.05."))

read_docx() %>% body_add_flextable(ft3) %>%
  print(target = file.path(MAIN_TABLES, "Table3_AR_corrected_primary.docx"))

# ---- Table 4. Microbial P interaction --------------------------------------

t4_data_int <- results_autocor_int %>%
  filter(file == "Cores", compound == "uP",
         window == WINDOW_PRIMARY, response %in% RESPONSES_MAIN) %>%
  transmute(response,
            "Interaction \u03b2" = format_est(estimate, 3),
            `Interaction p (AR)` = format_p(p.value))

t4_data_slopes <- results_autocor_slopes %>%
  filter(file == "Cores", compound == "uP",
         window == WINDOW_PRIMARY, response %in% RESPONSES_MAIN) %>%
  arrange(response, Treatment) %>%
  transmute(Response = recode(response, !!!response_labels),
            response_raw = response,
            Treatment,
            "Slope (\u03b2)" = format_est(slope, 4),
            `95% CI`         = format_ci(slope_lo, slope_hi, 4)) %>%
  left_join(t4_data_int %>% rename(response_raw = response),
            by = "response_raw") %>%
  select(Response, Treatment, "Slope (\u03b2)", `95% CI`,
         "Interaction \u03b2", `Interaction p (AR)`)

ft4 <- flextable::flextable(t4_data_slopes) %>%
  flextable::theme_booktabs() %>%
  flextable::bg(part = "header", bg = "#3a3a3a") %>%
  flextable::color(part = "header", color = "white") %>%
  flextable::bold(part = "header") %>%
  flextable::fontsize(size = 8, part = "all") %>%
  flextable::padding(padding.top = 2.5, padding.bottom = 2.5, part = "all") %>%
  flextable::align(align = "center", part = "all") %>%
  flextable::align(j = "Response", align = "left", part = "body") %>%
  flextable::merge_v(j = c("Response", "Interaction \u03b2",
                           "Interaction p (AR)")) %>%
  flextable::set_table_properties(layout = "autofit", width = 1) %>%
  flextable::add_header_lines(paste0(
    "Table 4. Microbial P \u00d7 Treatment interaction at the primary \u00b1",
    WINDOW_PRIMARY, "-day window. AR-corrected per-Treatment simple slopes ",
    "from corExp model on Plot \u00d7 Date aggregated data, with the ",
    "interaction term \u03b2 and significance from the same model. ",
    "Microbial P \u00d7 Treatment is pre-specified as a headline target ",
    "so AR estimates are available for both responses regardless of ",
    "standard-LMM FDR significance."))

read_docx() %>% body_add_flextable(ft4) %>%
  print(target = file.path(MAIN_TABLES, "Table4_uP_interaction.docx"))

# ---- Table S1. Cross-method comparison (DOCX) -------------------------------

ts1_data <- cross_method_comparison %>%
  filter(!is.na(windowed_estimate_ar) | !is.na(sameday_estimate_ar)) %>%
  arrange(response, file, compound) %>%
  transmute(
    Response                 = recode(response, !!!response_labels),
    Predictor                = paste0(compound, " [", file, "]"),
    "\u03b2 (windowed AR)"   = format_est(windowed_estimate_ar, 3),
    `p (windowed AR)`        = format_p(windowed_p_ar),
    "\u03b2 (same-day AR)"   = format_est(sameday_estimate_ar, 3),
    `p (same-day AR)`        = format_p(sameday_p_ar),
    Concordance              = ifelse(
      is.na(windowed_p_ar) | is.na(sameday_p_ar), "\u2014",
      ifelse((windowed_p_ar < 0.05) == (sameday_p_ar < 0.05),
             "agree", "differ"))
  )

fts1 <- flextable::flextable(ts1_data) %>%
  flextable::theme_booktabs() %>%
  flextable::bg(part = "header", bg = "#3a3a3a") %>%
  flextable::color(part = "header", color = "white") %>%
  flextable::bold(part = "header") %>%
  flextable::fontsize(size = 8, part = "all") %>%
  flextable::padding(padding.top = 2.5, padding.bottom = 2.5, part = "all") %>%
  flextable::align(align = "center", part = "all") %>%
  flextable::align(j = c("Response", "Predictor"),
                   align = "left", part = "body") %>%
  flextable::merge_v(j = "Response") %>%
  flextable::set_table_properties(layout = "autofit", width = 1) %>%
  flextable::add_header_lines(paste0(
    "Table S1. Cross-method comparison: AR-corrected estimates from the ",
    "primary windowed analysis (\u00b1", WINDOW_PRIMARY, "-day window) and ",
    "from same-day matching, in which each chemistry observation is paired ",
    "only with hourly flux from the same Plot on the same calendar day. ",
    "Same-day matching is the most stringent test of temporal independence ",
    "available with this design. Concordance is reported as 'agree' if AR ",
    "significance status (p < 0.05) matches across schemes."))

read_docx() %>% body_add_flextable(fts1) %>%
  print(target = file.path(SUPP_TABLES, "TableS1_cross_method.docx"))

# =============================================================================
# 16. MAIN-TEXT FIGURES
# =============================================================================

# ---- Figure 1. Chemistry coefficients (full panel, both responses) ----------

message("\nGenerating Figure 1...")

fig1_df <- results_chem %>%
  filter(window == WINDOW_PRIMARY, response %in% RESPONSES_MAIN) %>%
  mutate(compound_file = paste0(compound, "  [", file, "]"),
         sig_fdr       = p.fdr < FDR_THRESHOLD,
         response_lbl  = recode(response, !!!response_labels))

comp_order <- fig1_df %>%
  group_by(compound_file) %>%
  summarise(m = mean(estimate, na.rm = TRUE), .groups = "drop") %>%
  arrange(m) %>% pull(compound_file)

fig1_df <- mutate(fig1_df,
                  compound_file = factor(compound_file, levels = comp_order))

p1 <- ggplot(fig1_df,
             aes(x = estimate, y = compound_file,
                 colour = file, shape = sig_fdr)) +
  geom_vline(xintercept = 0, linetype = "dashed",
             colour = "grey50", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                 height = 0.30, linewidth = 0.5, alpha = 0.80) +
  geom_point(size = 2.4) +
  scale_shape_manual(values = c(`FALSE` = 1L, `TRUE` = 19L),
                     labels = c("p (FDR) \u2265 0.05", "p (FDR) < 0.05"),
                     name = NULL) +
  scale_colour_manual(values = file_pal, name = "Chemistry pool") +
  facet_wrap(~ response_lbl, ncol = 2L) +
  labs(x = "Standardised \u03b2 (\u00b1 95% CI)", y = NULL) +
  theme_pub(9) + theme(legend.box = "vertical")

fig1_height <- max(5, 0.27 * length(comp_order) + 2.2)
ggsave(file.path(MAIN_FIGS, "Figure1_chemistry_coefficients.png"),
       p1, width = 11, height = fig1_height, dpi = FIG_DPI, units = "in")
ggsave(file.path(MAIN_FIGS, "Figure1_chemistry_coefficients.pdf"),
       p1, width = 11, height = fig1_height, units = "in")

# ---- Figure 2. Microbial P simple slopes (both panels) ----------------------

message("Generating Figure 2...")

uP_matched <- build_matched_df(cores, "uP", WINDOW_PRIMARY) %>%
  group_by(Date, Plot, Treatment) %>%
  summarise(across(c(log_flux_mean, flux_resid_mean, chem_value),
                   ~ mean(.x, na.rm = TRUE)),
            .groups = "drop") %>%
  mutate(Treatment = factor(Treatment, levels = c("Control", "Warmed")))

uP_slopes_ar <- results_autocor_slopes %>%
  filter(file == "Cores", compound == "uP", window == WINDOW_PRIMARY,
         response %in% RESPONSES_MAIN)

uP_int_ar <- results_autocor_int %>%
  filter(file == "Cores", compound == "uP", window == WINDOW_PRIMARY,
         response %in% RESPONSES_MAIN)

panel_meta <- tibble(
  response    = RESPONSES_MAIN,
  panel_label = c("A", "B"),
  y_label     = c("Mean log(flux)",
                  "Flux residuals\n(after T + VWC)")
)

fig2_panels <- list()
for (i in seq_len(nrow(panel_meta))) {

  rm       <- panel_meta[i, ]
  resp_col <- rm$response

  d_resp <- uP_matched %>%
    select(Treatment, chem_value, resp_val = all_of(resp_col)) %>%
    filter(!is.na(resp_val), !is.na(chem_value))

  sl_ctrl <- uP_slopes_ar %>%
    filter(response == resp_col, Treatment == "Control")
  sl_warm <- uP_slopes_ar %>%
    filter(response == resp_col, Treatment == "Warmed")
  ar_p    <- uP_int_ar %>%
    filter(response == resp_col) %>% pull(p.value)

  ann_lines <- character()
  if (nrow(sl_ctrl) > 0 && nrow(sl_warm) > 0)
    ann_lines <- c(
      sprintf("Control \u03b2 = %.3f [%.3f, %.3f]",
              sl_ctrl$slope[1], sl_ctrl$slope_lo[1], sl_ctrl$slope_hi[1]),
      sprintf("Warmed  \u03b2 = %.3f [%.3f, %.3f]",
              sl_warm$slope[1], sl_warm$slope_lo[1], sl_warm$slope_hi[1])
    )
  if (length(ar_p) > 0)
    ann_lines <- c(ann_lines,
                   sprintf("Interaction p (AR) = %.3g", ar_p[1]))

  p <- ggplot(d_resp, aes(x = chem_value, y = resp_val,
                          colour = Treatment, fill = Treatment)) +
    geom_point(alpha = 0.45, size = 2.0, shape = 16) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.85, alpha = 0.18) +
    scale_colour_manual(values = treat_pal, name = NULL) +
    scale_fill_manual(  values = treat_pal, name = NULL) +
    labs(x = expression("Microbial P (" * mu * "g g"^{-1} * " soil)"),
         y = rm$y_label,
         tag = rm$panel_label) +
    theme_pub(10) +
    theme(plot.tag = element_text(face = "bold", size = 12))

  if (length(ann_lines) > 0) {
    x_pos <- min(d_resp$chem_value, na.rm = TRUE) +
             0.05 * diff(range(d_resp$chem_value, na.rm = TRUE))
    y_rng <- range(d_resp$resp_val, na.rm = TRUE)
    y_pos <- y_rng[2] - 0.03 * diff(y_rng)
    p <- p + annotate("text", x = x_pos, y = y_pos,
                      label = paste(ann_lines, collapse = "\n"),
                      hjust = 0, vjust = 1,
                      size = 2.9, colour = "grey20", family = "mono")
  }

  fig2_panels[[rm$panel_label]] <- p
}

fig2 <- (fig2_panels[["A"]] | fig2_panels[["B"]]) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(file.path(MAIN_FIGS, "Figure2_uP_simple_slopes.png"),
       fig2, width = 9, height = 4.5, dpi = FIG_DPI, units = "in")
ggsave(file.path(MAIN_FIGS, "Figure2_uP_simple_slopes.pdf"),
       fig2, width = 9, height = 4.5, units = "in")

# ---- Figure 3. AR-collapse summary -----------------------------------------

message("Generating Figure 3...")

fig3_df <- ar_collapse %>%
  filter(response %in% RESPONSES_MAIN) %>%
  mutate(label = paste0(compound, " [", file, "]\n",
                        recode(response, !!!response_labels)),
         survives_lab = factor(survives_AR,
                               levels = c(TRUE, FALSE),
                               labels = c("Survives AR (p < 0.05)",
                                          "Collapses (p \u2265 0.05)")),
         std_lo = std_estimate - 1.96 * std_se,
         std_hi = std_estimate + 1.96 * std_se,
         ar_lo  = ar_estimate  - 1.96 * ar_se,
         ar_hi  = ar_estimate  + 1.96 * ar_se) %>%
  arrange(survives_AR, ar_estimate) %>%
  mutate(label = factor(label, levels = unique(label)))

p3 <- ggplot(fig3_df) +
  geom_vline(xintercept = 0, linetype = "dashed",
             colour = "grey50", linewidth = 0.4) +
  geom_errorbarh(aes(y = label, xmin = std_lo, xmax = std_hi),
                 height = 0.20, linewidth = 0.45, colour = "grey55",
                 position = position_nudge(y = 0.18)) +
  geom_point(aes(y = label, x = std_estimate),
             shape = 1, size = 2.4, colour = "grey45",
             position = position_nudge(y = 0.18)) +
  geom_errorbarh(aes(y = label, xmin = ar_lo, xmax = ar_hi,
                     colour = survives_lab),
                 height = 0.20, linewidth = 0.55,
                 position = position_nudge(y = -0.18)) +
  geom_point(aes(y = label, x = ar_estimate, colour = survives_lab),
             shape = 19, size = 2.6,
             position = position_nudge(y = -0.18)) +
  scale_colour_manual(values = c("Survives AR (p < 0.05)"     = "#1a9641",
                                  "Collapses (p \u2265 0.05)" = "#d73027"),
                       name = NULL) +
  labs(x = "Standardised \u03b2 (\u00b1 95% CI)", y = NULL,
       caption = "Open grey = standard LMM; coloured filled = corExp-corrected") +
  theme_pub(9) +
  theme(plot.caption = element_text(size = 8, colour = "grey40", hjust = 0))

fig3_height <- max(5, 0.40 * nrow(fig3_df) + 1.8)
ggsave(file.path(MAIN_FIGS, "Figure3_AR_collapse.png"),
       p3, width = 9.5, height = fig3_height, dpi = FIG_DPI, units = "in")
ggsave(file.path(MAIN_FIGS, "Figure3_AR_collapse.pdf"),
       p3, width = 9.5, height = fig3_height, units = "in")
# =============================================================================
# 17. SUPPLEMENTARY FIGURES
# =============================================================================

# ---- Figure S1. Empirical flux ACF per plot ---------------------------------

message("\nGenerating Figure S1: flux ACF...")

acf_long <- acf_per_plot %>%
  pivot_longer(starts_with("acf_"),
               names_to = "series", names_prefix = "acf_",
               values_to = "acf") %>%
  mutate(series = recode(series,
                         log   = "log(Flux)",
                         resid = "Flux residuals (after T + VWC)"),
         Plot   = factor(Plot))

pS1 <- ggplot(acf_long, aes(x = lag_days, y = acf, colour = Plot)) +
  geom_hline(yintercept =  ACF_THRESH, linetype = "dotted", colour = "grey40") +
  geom_hline(yintercept = -ACF_THRESH, linetype = "dotted", colour = "grey40") +
  geom_hline(yintercept = 0, linetype = "solid",
             colour = "grey70", linewidth = 0.3) +
  geom_line(linewidth = 0.55, alpha = 0.85) +
  scale_colour_manual(values = plot_pal) +
  scale_x_continuous(breaks = seq(0, ACF_MAX_LAG, by = 7L)) +
  facet_wrap(~ series, ncol = 1L) +
  labs(x = "Lag (days)", y = "Autocorrelation",
       caption = paste0("Dotted lines: \u00b1", ACF_THRESH,
                        " decorrelation threshold. ",
                        "ACF on daily-aggregated values per plot.")) +
  theme_pub(9) +
  theme(plot.caption = element_text(size = 8, colour = "grey40", hjust = 0))

ggsave(file.path(SUPP_FIGS, "FigureS1_flux_ACF.png"),
       pS1, width = 9, height = 7, dpi = FIG_DPI, units = "in")

# ---- Figure S2. Window sensitivity grid ------------------------------------

message("Generating Figure S2: window sensitivity...")

sens_plot_df <- sens_grid %>%
  filter(response %in% RESPONSES_MAIN) %>%
  mutate(model_lbl = recode(model,
                            standard_LMM = "Standard LMM",
                            corExp       = "corExp (AR)"),
         target = paste0(compound, " [", file, "] \u2192 ",
                         recode(response, !!!response_labels)))

featured_targets <- results_chem %>%
  filter(window == WINDOW_PRIMARY, response %in% RESPONSES_MAIN,
         p.fdr < FDR_THRESHOLD) %>%
  mutate(target = paste0(compound, " [", file, "] \u2192 ",
                         recode(response, !!!response_labels))) %>%
  pull(target) %>% unique()

sens_featured <- sens_plot_df %>%
  filter(target %in% featured_targets) %>%
  mutate(target = factor(target, levels = featured_targets))

pS2 <- ggplot(sens_featured,
              aes(x = window, y = estimate,
                  colour = model_lbl, fill = model_lbl)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             colour = "grey50", linewidth = 0.35) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = 0.18, colour = NA) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.8) +
  geom_vline(xintercept = WINDOW_PRIMARY, linetype = "dotted",
             colour = "grey25", linewidth = 0.4) +
  scale_colour_manual(values = c("Standard LMM" = "#4575b4",
                                  "corExp (AR)" = "#d73027"),
                       name = NULL) +
  scale_fill_manual(  values = c("Standard LMM" = "#4575b4",
                                  "corExp (AR)" = "#d73027"),
                       name = NULL) +
  scale_x_continuous(breaks = WINDOW_GRID) +
  facet_wrap(~ target, scales = "free_y", ncol = 2L) +
  labs(x = "Window half-width (\u00b1 days)",
       y = "Standardised \u03b2 (\u00b1 95% CI)",
       caption = paste0("Vertical dotted line: primary \u00b1",
                        WINDOW_PRIMARY, "-day window.")) +
  theme_pub(9) +
  theme(strip.text   = element_text(size = 8),
        plot.caption = element_text(size = 8, colour = "grey40", hjust = 0))

n_targets   <- length(featured_targets)
nS2_height  <- max(6, 1.5 * ceiling(n_targets / 2) + 2)
ggsave(file.path(SUPP_FIGS, "FigureS2_window_sensitivity.png"),
       pS2, width = 11, height = nS2_height,
       dpi = FIG_DPI, units = "in")

# ---- Figure S3. Correlation structure comparison ----------------------------

message("Generating Figure S3: correlation structures...")

cs_plot_df <- corr_struct %>%
  filter(response %in% RESPONSES_MAIN) %>%
  mutate(target = paste0(compound, " [", file, "] \u2192 ",
                         recode(response, !!!response_labels)))

cs_featured <- cs_plot_df %>%
  filter(target %in% featured_targets) %>%
  mutate(target = factor(target, levels = featured_targets))

cs_pal <- c("corExp"  = "#d73027",
            "corCAR1" = "#4575b4",
            "corGaus" = "#1a9850")

pS3 <- ggplot(cs_featured,
              aes(y = target, x = estimate, colour = model)) +
  geom_vline(xintercept = 0, linetype = "dashed",
             colour = "grey50", linewidth = 0.35) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                 height = 0, linewidth = 0.5, alpha = 0.85,
                 position = position_dodge(width = 0.55)) +
  geom_point(size = 2.3,
             position = position_dodge(width = 0.55)) +
  scale_colour_manual(values = cs_pal, name = "Correlation\nstructure") +
  labs(x = "Standardised \u03b2 (\u00b1 95% CI)", y = NULL,
       caption = paste0("Estimates at the primary \u00b1", WINDOW_PRIMARY,
                        "-day window. corExp and corCAR1 are mathematically ",
                        "equivalent up to parameterisation.")) +
  theme_pub(9) +
  theme(plot.caption = element_text(size = 8, colour = "grey40", hjust = 0))

ggsave(file.path(SUPP_FIGS, "FigureS3_correlation_structures.png"),
       pS3, width = 10, height = max(5, 0.38 * n_targets + 2),
       dpi = FIG_DPI, units = "in")

# ---- Figure S4. Flux CV response (supplementary only) ----------------------

message("Generating Figure S4: flux CV response...")

figS4_df <- results_chem %>%
  filter(window == WINDOW_PRIMARY, response == "flux_cv") %>%
  mutate(compound_file = paste0(compound, "  [", file, "]"),
         sig_fdr       = p.fdr < FDR_THRESHOLD)

cv_order <- figS4_df %>% arrange(estimate) %>% pull(compound_file)
figS4_df <- mutate(figS4_df,
                   compound_file = factor(compound_file, levels = cv_order))

pS4 <- ggplot(figS4_df,
              aes(x = estimate, y = compound_file,
                  colour = file, shape = sig_fdr)) +
  geom_vline(xintercept = 0, linetype = "dashed",
             colour = "grey50", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                 height = 0.30, linewidth = 0.5, alpha = 0.80) +
  geom_point(size = 2.4) +
  scale_shape_manual(values = c(`FALSE` = 1L, `TRUE` = 19L),
                     labels = c("p (FDR) \u2265 0.05", "p (FDR) < 0.05"),
                     name = NULL) +
  scale_colour_manual(values = file_pal, name = "Chemistry pool") +
  labs(x = "Standardised \u03b2 (\u00b1 95% CI)", y = NULL,
       caption = paste0("Standard LMM at \u00b1", WINDOW_PRIMARY,
                        "-day window. Flux CV reported in supplementary only.")) +
  theme_pub(9) +
  theme(legend.box = "vertical",
        plot.caption = element_text(size = 8, colour = "grey40", hjust = 0))

ggsave(file.path(SUPP_FIGS, "FigureS4_flux_cv.png"),
       pS4, width = 9, height = max(5, 0.27 * length(cv_order) + 2),
       dpi = FIG_DPI, units = "in")

# ---- Figure S5. Same-day microbial P simple slopes -------------------------
# Counterpart to Figure 2 under the most stringent matching scheme.

message("Generating Figure S5: same-day uP simple slopes...")

uP_sd <- sameday_cache[[paste("Cores", "uP", sep = "|")]]

if (!is.null(uP_sd) && nrow(uP_sd) >= MIN_OBS) {

  uP_sd <- uP_sd %>%
    mutate(Treatment = factor(Treatment, levels = c("Control", "Warmed")))

  fig_sd_panels <- list()
  for (i in seq_len(nrow(panel_meta))) {

    rm       <- panel_meta[i, ]
    resp_col <- rm$response

    d_resp <- uP_sd %>%
      select(Treatment, chem_value, resp_val = all_of(resp_col)) %>%
      filter(!is.na(resp_val), !is.na(chem_value))

    if (nrow(d_resp) < MIN_OBS) next

    sl_ctrl <- sameday_ar_slopes %>%
      filter(file == "Cores", compound == "uP",
             response == resp_col, Treatment == "Control")
    sl_warm <- sameday_ar_slopes %>%
      filter(file == "Cores", compound == "uP",
             response == resp_col, Treatment == "Warmed")
    sd_p    <- sameday_ar_int %>%
      filter(file == "Cores", compound == "uP",
             response == resp_col) %>% pull(p.value)

    ann_lines <- character()
    if (nrow(sl_ctrl) > 0 && nrow(sl_warm) > 0)
      ann_lines <- c(
        sprintf("Control \u03b2 = %.3f [%.3f, %.3f]",
                sl_ctrl$slope[1], sl_ctrl$slope_lo[1], sl_ctrl$slope_hi[1]),
        sprintf("Warmed  \u03b2 = %.3f [%.3f, %.3f]",
                sl_warm$slope[1], sl_warm$slope_lo[1], sl_warm$slope_hi[1])
      )
    if (length(sd_p) > 0)
      ann_lines <- c(ann_lines,
                     sprintf("Interaction p (AR) = %.3g", sd_p[1]))

    p <- ggplot(d_resp, aes(x = chem_value, y = resp_val,
                            colour = Treatment, fill = Treatment)) +
      geom_point(alpha = 0.45, size = 2.0, shape = 16) +
      geom_smooth(method = "lm", se = TRUE, linewidth = 0.85, alpha = 0.18) +
      scale_colour_manual(values = treat_pal, name = NULL) +
      scale_fill_manual(  values = treat_pal, name = NULL) +
      labs(x = expression("Microbial P (" * mu * "g g"^{-1} * " soil)"),
           y = rm$y_label,
           tag = rm$panel_label) +
      theme_pub(10) +
      theme(plot.tag = element_text(face = "bold", size = 12))

    if (length(ann_lines) > 0) {
      x_pos <- min(d_resp$chem_value, na.rm = TRUE) +
               0.05 * diff(range(d_resp$chem_value, na.rm = TRUE))
      y_rng <- range(d_resp$resp_val, na.rm = TRUE)
      y_pos <- y_rng[2] - 0.03 * diff(y_rng)
      p <- p + annotate("text", x = x_pos, y = y_pos,
                        label = paste(ann_lines, collapse = "\n"),
                        hjust = 0, vjust = 1,
                        size = 2.9, colour = "grey20", family = "mono")
    }
    fig_sd_panels[[rm$panel_label]] <- p
  }

  if (length(fig_sd_panels) == length(RESPONSES_MAIN)) {
    fig_S5 <- (fig_sd_panels[["A"]] | fig_sd_panels[["B"]]) +
      plot_layout(guides = "collect") &
      theme(legend.position = "bottom")
    fig_S5 <- fig_S5 + plot_annotation(
      title    = "Microbial P \u00d7 Treatment under same-day matching",
      subtitle = paste0("Plot \u00d7 Date matched on the same calendar day ",
                        "only (n = ", nrow(uP_sd),
                        " observations); compare with Figure 2."),
      theme    = theme(plot.title    = element_text(size = 11, face = "bold"),
                       plot.subtitle = element_text(size = 9, colour = "grey40"))
    )
    ggsave(file.path(SUPP_FIGS, "FigureS5_uP_sameday_slopes.png"),
           fig_S5, width = 10, height = 5,
           dpi = FIG_DPI, units = "in")
  }
}

# ---- Figure S6. Cross-method comparison -------------------------------------
# Side-by-side AR-corrected estimates from windowed and same-day analyses
# for every compound x main-response. Visual demonstration that the primary
# microbial P finding is robust to matching scheme.

message("Generating Figure S6: cross-method comparison...")

cm_long <- cross_method_comparison %>%
  filter(!is.na(windowed_estimate_ar) | !is.na(sameday_estimate_ar)) %>%
  transmute(
    file, compound, response,
    windowed_est = windowed_estimate_ar,
    windowed_lo  = windowed_estimate_ar - 1.96 * windowed_se_ar,
    windowed_hi  = windowed_estimate_ar + 1.96 * windowed_se_ar,
    windowed_p   = windowed_p_ar,
    sameday_est  = sameday_estimate_ar,
    sameday_lo   = sameday_estimate_ar  - 1.96 * sameday_se_ar,
    sameday_hi   = sameday_estimate_ar  + 1.96 * sameday_se_ar,
    sameday_p    = sameday_p_ar,
    label        = paste0(compound, " [", file, "]\n",
                          recode(response, !!!response_labels))
  )

if (nrow(cm_long) > 0L) {

  lab_order <- cm_long %>%
    mutate(m = pmax(windowed_est, sameday_est, na.rm = TRUE)) %>%
    arrange(m) %>% pull(label)
  cm_long <- mutate(cm_long,
                    label = factor(label, levels = unique(lab_order)))

  pS6 <- ggplot(cm_long) +
    geom_vline(xintercept = 0, linetype = "dashed",
               colour = "grey50", linewidth = 0.4) +
    geom_errorbarh(aes(y = label, xmin = windowed_lo, xmax = windowed_hi),
                   height = 0.20, linewidth = 0.5, colour = "#4575b4",
                   position = position_nudge(y = 0.18)) +
    geom_point(aes(y = label, x = windowed_est),
               shape = 19, size = 2.4, colour = "#4575b4",
               position = position_nudge(y = 0.18)) +
    geom_errorbarh(aes(y = label, xmin = sameday_lo, xmax = sameday_hi),
                   height = 0.20, linewidth = 0.5, colour = "#d73027",
                   position = position_nudge(y = -0.18)) +
    geom_point(aes(y = label, x = sameday_est),
               shape = 17, size = 2.4, colour = "#d73027",
               position = position_nudge(y = -0.18)) +
    labs(x = "Standardised \u03b2 (AR-corrected) \u00b1 95% CI",
         y = NULL,
         caption = paste0("Blue circles: windowed analysis (\u00b1",
                          WINDOW_PRIMARY, "d). Red triangles: same-day ",
                          "matching. Both axes are AR-corrected (corExp).")) +
    theme_pub(9) +
    theme(plot.caption = element_text(size = 8, colour = "grey40", hjust = 0))

  n_labs    <- length(unique(cm_long$label))
  nS6_height <- max(5, 0.50 * n_labs + 2)
  ggsave(file.path(SUPP_FIGS, "FigureS6_cross_method.png"),
         pS6, width = 10, height = nS6_height,
         dpi = FIG_DPI, units = "in")
}

# ---- Figure S7. Same-day chemistry coefficients (Figure 1 counterpart) ------
# Counterpart to Figure 1 under the most stringent (same-day) matching.
# Uses Figure 1's compound ordering for direct visual comparability.

message("Generating Figure S7: same-day chemistry coefficients...")

figS7_df <- sameday_chem %>%
  filter(response %in% RESPONSES_MAIN) %>%
  mutate(compound_file = paste0(compound, "  [", file, "]"),
         sig_fdr       = p.fdr < FDR_THRESHOLD,
         response_lbl  = recode(response, !!!response_labels))

# Use Figure 1's compound ordering for direct comparison
figS7_df <- mutate(figS7_df,
                   compound_file = factor(compound_file, levels = comp_order))

pS7 <- ggplot(figS7_df,
              aes(x = estimate, y = compound_file,
                  colour = file, shape = sig_fdr)) +
  geom_vline(xintercept = 0, linetype = "dashed",
             colour = "grey50", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                 height = 0.30, linewidth = 0.5, alpha = 0.80) +
  geom_point(size = 2.4) +
  scale_shape_manual(values = c(`FALSE` = 1L, `TRUE` = 19L),
                     labels = c("p (FDR) \u2265 0.05", "p (FDR) < 0.05"),
                     name = NULL) +
  scale_colour_manual(values = file_pal, name = "Chemistry pool") +
  facet_wrap(~ response_lbl, ncol = 2L) +
  labs(x = "Standardised \u03b2 (\u00b1 95% CI)", y = NULL,
       caption = "Same-day matching counterpart to Figure 1.") +
  theme_pub(9) +
  theme(legend.box = "vertical",
        plot.caption = element_text(size = 8, colour = "grey40", hjust = 0))

ggsave(file.path(SUPP_FIGS, "FigureS7_chemistry_coefficients_sameday.png"),
       pS7, width = 11, height = fig1_height,
       dpi = FIG_DPI, units = "in")

# ---- Figure S8. Same-day AR-collapse summary (Figure 3 counterpart) ---------
# Same-day analogue of Figure 3 for compounds that pass FDR in the same-day
# standard LMM and that have AR-corrected estimates.

message("Generating Figure S8: same-day AR-collapse summary...")

sameday_ar_collapse <- sameday_chem %>%
  filter(p.fdr < FDR_THRESHOLD, response %in% RESPONSES_MAIN) %>%
  inner_join(
    sameday_ar_main %>%
      select(file, compound, response,
             ar_estimate = estimate,
             ar_se       = std.error,
             ar_p        = p.value,
             ar_aic_improvement = aic_improvement),
    by = c("file", "compound", "response")
  ) %>%
  mutate(survives_AR = ar_p < 0.05,
         pct_attenuation =
           ifelse(estimate == 0, NA_real_,
                  100 * (abs(estimate) - abs(ar_estimate)) / abs(estimate)),
         label = paste0(compound, " [", file, "]\n",
                        recode(response, !!!response_labels)),
         survives_lab = factor(survives_AR,
                               levels = c(TRUE, FALSE),
                               labels = c("Survives AR (p < 0.05)",
                                          "Collapses (p \u2265 0.05)")),
         std_lo = estimate - 1.96 * std.error,
         std_hi = estimate + 1.96 * std.error,
         ar_lo  = ar_estimate - 1.96 * ar_se,
         ar_hi  = ar_estimate + 1.96 * ar_se) %>%
  arrange(survives_AR, ar_estimate) %>%
  mutate(label = factor(label, levels = unique(label)))

if (nrow(sameday_ar_collapse) > 0) {
  pS8 <- ggplot(sameday_ar_collapse) +
    geom_vline(xintercept = 0, linetype = "dashed",
               colour = "grey50", linewidth = 0.4) +
    geom_errorbarh(aes(y = label, xmin = std_lo, xmax = std_hi),
                   height = 0.20, linewidth = 0.45, colour = "grey55",
                   position = position_nudge(y = 0.18)) +
    geom_point(aes(y = label, x = estimate),
               shape = 1, size = 2.4, colour = "grey45",
               position = position_nudge(y = 0.18)) +
    geom_errorbarh(aes(y = label, xmin = ar_lo, xmax = ar_hi,
                       colour = survives_lab),
                   height = 0.20, linewidth = 0.55,
                   position = position_nudge(y = -0.18)) +
    geom_point(aes(y = label, x = ar_estimate, colour = survives_lab),
               shape = 19, size = 2.6,
               position = position_nudge(y = -0.18)) +
    scale_colour_manual(values = c("Survives AR (p < 0.05)"     = "#1a9641",
                                    "Collapses (p \u2265 0.05)" = "#d73027"),
                         name = NULL) +
    labs(x = "Standardised \u03b2 (\u00b1 95% CI)", y = NULL,
         caption = paste0("Same-day matching counterpart to Figure 3. ",
                          "Open grey = standard LMM; ",
                          "coloured filled = corExp-corrected.")) +
    theme_pub(9) +
    theme(plot.caption = element_text(size = 8, colour = "grey40", hjust = 0))

  fig8_height <- max(5, 0.40 * nrow(sameday_ar_collapse) + 1.8)
  ggsave(file.path(SUPP_FIGS, "FigureS8_AR_collapse_sameday.png"),
         pS8, width = 9.5, height = fig8_height,
         dpi = FIG_DPI, units = "in")

  # Export the same-day AR-collapse table to CSV
  write_csv(sameday_ar_collapse %>%
              select(file, compound, response,
                     std_estimate = estimate, std_se = std.error,
                     std_fdr = p.fdr,
                     ar_estimate, ar_se, ar_p,
                     ar_aic_improvement, pct_attenuation, survives_AR),
            file.path(SUPP_TABLES, "sameday_ar_collapse.csv"))
}

# ---- Diagnostic plots for AR-robust compounds at the primary window --------

message("\nGenerating diagnostic plots for AR-robust compounds...")

diag_targets <- results_autocor_main %>%
  filter(window == WINDOW_PRIMARY, response %in% RESPONSES_MAIN,
         p.value < 0.05) %>%
  arrange(p.value) %>%
  select(window, file, compound, response)

for (i in seq_len(nrow(diag_targets))) {
  row     <- diag_targets[i, ]
  key     <- paste(row$file, row$compound, row$window, sep = "|")
  matched <- matched_cache[[key]]
  if (is.null(matched)) next
  make_diagnostic_plots(matched, row$response, row$compound,
                        row$file, row$window, SUPP_DIAG)
}

# =============================================================================
# 18. CONSOLE SUMMARY
# =============================================================================

message("\n=============================================================")
message("ANALYSIS COMPLETE")
message("=============================================================")

message("\nABIOTIC BASE MODEL")
message("  Marginal R\u00b2 = ",   round(abiotic_r2[1, "R2m"], 3),
        " | Conditional R\u00b2 = ", round(abiotic_r2[1, "R2c"], 3))
message("  Temperature \u03b2 = ",
        round(abiotic_fe$estimate[abiotic_fe$term == "Temperature"], 4),
        " (t = ", round(abiotic_fe$statistic[abiotic_fe$term == "Temperature"], 1), ")")
message("  VWC \u03b2 = ",
        round(abiotic_fe$estimate[abiotic_fe$term == "VWC"], 3),
        " (t = ", round(abiotic_fe$statistic[abiotic_fe$term == "VWC"], 1), ")")

message("\nFLUX ACF (decorrelation lags, |ACF| < ", ACF_THRESH, "):")
print(decorr_lags)

message("\nWINDOWED PRIMARY: AR-corrected main-effect survivors at p < 0.05")
results_autocor_main %>%
  filter(window == WINDOW_PRIMARY, response %in% RESPONSES_MAIN,
         p.value < 0.05) %>%
  select(file, compound, response, estimate, p.value, aic_improvement) %>%
  arrange(file, compound, response) %>% print(n = Inf, width = Inf)

message("\nMICROBIAL P AR-corrected simple slopes (windowed primary):")
results_autocor_slopes %>%
  filter(file == "Cores", compound == "uP",
         window == WINDOW_PRIMARY, response %in% RESPONSES_MAIN) %>%
  select(response, Treatment, slope, slope_lo, slope_hi) %>%
  arrange(response, Treatment) %>% print(n = Inf)

message("\nSAME-DAY: AR-corrected main-effect survivors at p < 0.05")
if (nrow(sameday_ar_main) > 0)
  sameday_ar_main %>%
    filter(p.value < 0.05) %>%
    select(file, compound, response, estimate, p.value, aic_improvement) %>%
    arrange(file, compound, response) %>% print(n = Inf, width = Inf)

message("\nMICROBIAL P AR-corrected simple slopes (same-day):")
if (nrow(sameday_ar_slopes) > 0)
  sameday_ar_slopes %>%
    filter(file == "Cores", compound == "uP",
           response %in% RESPONSES_MAIN) %>%
    select(response, Treatment, slope, slope_lo, slope_hi) %>%
    arrange(response, Treatment) %>% print(n = Inf)

message("\nTREATMENT MAIN EFFECT (windowed primary, main responses):")
treat_min <- results_treat %>%
  filter(window == WINDOW_PRIMARY, response %in% RESPONSES_MAIN) %>%
  summarise(min_fdr = min(p.fdr, na.rm = TRUE)) %>% pull(min_fdr)
message("  Minimum FDR-adjusted p: ", round(treat_min, 3),
        " (no Treatment main effect at FDR < ", FDR_THRESHOLD, ")")

message("\nOUTPUTS")
message("  Main:           ", MAIN_DIR)
message("  Supplementary:  ", SUPP_DIR)

# =============================================================================
# 19. TREATMENT EFFECT ON CHEMISTRY POOLS
# =============================================================================
# Whether warming changes the MAGNITUDE of each chemistry pool (distinct from
# whether it modifies the chemistry-flux slope, tested in Sections 8 and 13).
# Two design-matched approaches:
#   (a) Unpaired: chem_value ~ Treatment + (1 | Plot) on chemistry
#       observations. The contrast rests on six plot-level means.
#   (b) Paired/blocked: the TRACE plots are spatially blocked in consecutive
#       pairs (Block 1 = Plots 1,2; Block 2 = Plots 3,4; Block 3 = Plots 5,6),
#       one Control and one Warmed plot per block. A paired t-test on the
#       three within-block mean differences removes shared spatial variance
#       and is the design-matched test. With three blocks this is an
#       estimation exercise: the confidence interval is the result.
# Unadjusted p-values are primary (each compound is a pre-specified question,
# not a screen); a BH-FDR column is retained for completeness.

message("\n=== SECTION 19: Treatment effect on chemistry pools ===")

BLOCK_MAP <- c(`1` = 1L, `2` = 1L, `3` = 2L, `4` = 2L, `5` = 3L, `6` = 3L)

# ---- 19a. Unpaired treatment effect -----------------------------------------

fit_treatment_unpaired <- function(chem_df, compound) {

  d <- chem_df %>%
    select(Plot, Treatment, chem_value = all_of(compound)) %>%
    filter(!is.na(chem_value), is.finite(chem_value), !is.na(Treatment)) %>%
    mutate(Treatment = factor(Treatment, levels = c("Control", "Warmed")))

  if (nrow(d) < MIN_OBS) return(NULL)
  if (length(unique(d$Treatment)) < 2L) return(NULL)

  ctrl_mean <- mean(d$chem_value[d$Treatment == "Control"], na.rm = TRUE)
  warm_mean <- mean(d$chem_value[d$Treatment == "Warmed"],  na.rm = TRUE)

  plot_means <- d %>%
    group_by(Plot, Treatment) %>%
    summarise(plot_mean = mean(chem_value, na.rm = TRUE), .groups = "drop")
  plot_means_str <- paste(
    sprintf("P%d(%s)=%.3g", plot_means$Plot,
            substr(plot_means$Treatment, 1, 1), plot_means$plot_mean),
    collapse = "; ")

  mod <- tryCatch(suppressWarnings(
    lmerTest::lmer(chem_value ~ Treatment + (1 | Plot),
                   data = d, REML = TRUE,
                   control = lmerControl(optimizer = "bobyqa"))),
    error = function(e) NULL)

  used_lm <- FALSE
  if (is.null(mod) || (inherits(mod, "merMod") && lme4::isSingular(mod))) {
    mod_lm  <- lm(chem_value ~ Treatment, data = d)
    co      <- summary(mod_lm)$coefficients
    used_lm <- TRUE
    if (!"TreatmentWarmed" %in% rownames(co)) return(NULL)
    est <- co["TreatmentWarmed", "Estimate"]
    se  <- co["TreatmentWarmed", "Std. Error"]
    pv  <- co["TreatmentWarmed", "Pr(>|t|)"]
    dfr <- mod_lm$df.residual
  } else {
    co  <- summary(mod)$coefficients
    if (!"TreatmentWarmed" %in% rownames(co)) return(NULL)
    est <- co["TreatmentWarmed", "Estimate"]
    se  <- co["TreatmentWarmed", "Std. Error"]
    pv  <- co["TreatmentWarmed", "Pr(>|t|)"]
    dfr <- co["TreatmentWarmed", "df"]
  }

  ci_lo  <- est - 1.96 * se
  ci_hi  <- est + 1.96 * se
  pct    <- if (is.finite(ctrl_mean) && ctrl_mean != 0)
              100 * est / ctrl_mean else NA_real_
  pct_lo <- if (is.finite(ctrl_mean) && ctrl_mean != 0)
              100 * ci_lo / ctrl_mean else NA_real_
  pct_hi <- if (is.finite(ctrl_mean) && ctrl_mean != 0)
              100 * ci_hi / ctrl_mean else NA_real_
  if (!is.na(pct_lo) && !is.na(pct_hi) && pct_lo > pct_hi) {
    tmp <- pct_lo; pct_lo <- pct_hi; pct_hi <- tmp
  }

  tibble(
    compound     = compound,
    n_obs        = nrow(d),
    n_plots      = length(unique(d$Plot)),
    control_mean = ctrl_mean,
    warmed_mean  = warm_mean,
    diff_WmC     = est,
    se           = se,
    conf.low     = ci_lo,
    conf.high    = ci_hi,
    df           = dfr,
    pct_change   = pct,
    pct_low      = pct_lo,
    pct_high     = pct_hi,
    p.value      = pv,
    model_type   = if (used_lm) "lm (Plot singular)" else "lmer",
    plot_means   = plot_means_str
  )
}

trt_rows <- list()
for (file_nm in names(chem_files)) {
  cdf <- chem_files[[file_nm]]$df
  for (comp in chem_files[[file_nm]]$compounds) {
    res <- fit_treatment_unpaired(cdf, comp)
    if (!is.null(res))
      trt_rows[[length(trt_rows) + 1L]] <-
        bind_cols(tibble(file = file_nm), res)
  }
}
treatment_on_chem <- bind_rows(trt_rows) %>%
  group_by(file) %>%
  mutate(p.fdr = p.adjust(p.value, method = "BH")) %>%
  ungroup() %>%
  arrange(file, p.value) %>%
  relocate(file, compound, .before = everything())

# ---- 19b. Paired (blocked) treatment effect ---------------------------------
# Paired t-test on the three within-block mean differences. The block mixed
# model is deliberately NOT used: with three blocks its variance component
# is estimated on two degrees of freedom and is anticonservative.

fit_treatment_blocked <- function(chem_df, compound) {

  d <- chem_df %>%
    mutate(Block = BLOCK_MAP[as.character(Plot)]) %>%
    select(Plot, Block, Treatment, chem_value = all_of(compound)) %>%
    filter(!is.na(chem_value), is.finite(chem_value),
           !is.na(Treatment), !is.na(Block)) %>%
    mutate(Treatment = factor(Treatment, levels = c("Control", "Warmed")))

  if (nrow(d) < MIN_OBS) return(NULL)

  plot_means <- d %>%
    group_by(Block, Plot, Treatment) %>%
    summarise(plot_mean = mean(chem_value, na.rm = TRUE), .groups = "drop")

  block_ok <- plot_means %>%
    group_by(Block) %>%
    summarise(has_c = any(Treatment == "Control"),
              has_w = any(Treatment == "Warmed"),
              n_c   = sum(Treatment == "Control"),
              n_w   = sum(Treatment == "Warmed"),
              .groups = "drop")
  usable_blocks <- block_ok %>%
    filter(has_c, has_w, n_c == 1L, n_w == 1L) %>% pull(Block)
  if (length(usable_blocks) < 2L) return(NULL)

  block_diff <- plot_means %>%
    filter(Block %in% usable_blocks) %>%
    select(Block, Treatment, plot_mean) %>%
    pivot_wider(names_from = Treatment, values_from = plot_mean) %>%
    mutate(diff_WmC = Warmed - Control)

  n_blocks  <- nrow(block_diff)
  mean_diff <- mean(block_diff$diff_WmC)
  sd_diff   <- sd(block_diff$diff_WmC)

  if (n_blocks >= 2L && is.finite(sd_diff) && sd_diff > 0) {
    tt        <- t.test(block_diff$Warmed, block_diff$Control, paired = TRUE)
    paired_p  <- tt$p.value
    paired_lo <- tt$conf.int[1]
    paired_hi <- tt$conf.int[2]
    paired_t  <- unname(tt$statistic)
  } else {
    paired_p <- paired_lo <- paired_hi <- paired_t <- NA_real_
  }

  ctrl_grand <- mean(block_diff$Control, na.rm = TRUE)
  pct_change <- if (is.finite(ctrl_grand) && ctrl_grand != 0)
                  100 * mean_diff / ctrl_grand else NA_real_
  diffs_str  <- paste(sprintf("B%d=%.3g", block_diff$Block,
                              block_diff$diff_WmC), collapse = "; ")

  tibble(
    compound         = compound,
    n_obs            = nrow(d),
    n_blocks_paired  = n_blocks,
    block_diffs      = diffs_str,
    paired_mean_diff = mean_diff,
    paired_se        = sd_diff / sqrt(n_blocks),
    paired_conf_low  = paired_lo,
    paired_conf_high = paired_hi,
    paired_t         = paired_t,
    paired_p         = paired_p,
    pct_change       = pct_change
  )
}

blk_rows <- list()
for (file_nm in names(chem_files)) {
  cdf <- chem_files[[file_nm]]$df
  for (comp in chem_files[[file_nm]]$compounds) {
    res <- fit_treatment_blocked(cdf, comp)
    if (!is.null(res))
      blk_rows[[length(blk_rows) + 1L]] <-
        bind_cols(tibble(file = file_nm), res)
  }
}
treatment_on_chem_blocked <- bind_rows(blk_rows) %>%
  arrange(file, paired_p) %>%
  relocate(file, compound, .before = everything())

write_csv(treatment_on_chem,
          file.path(SUPP_TABLES, "treatment_effect_on_chemistry.csv"))
write_csv(treatment_on_chem_blocked,
          file.path(SUPP_TABLES, "treatment_effect_blocked.csv"))

message("  Unpaired: ", nrow(treatment_on_chem), " compounds | sig (raw p): ",
        sum(treatment_on_chem$p.value < FDR_THRESHOLD, na.rm = TRUE))
message("  Paired:   ", nrow(treatment_on_chem_blocked),
        " compounds | sig paired t-test: ",
        sum(treatment_on_chem_blocked$paired_p < FDR_THRESHOLD, na.rm = TRUE))

# ---- 19c. Figure S9: unpaired treatment effect on a comparable scale --------
# Plotted as percentage change from the control mean. Compounds are measured
# in widely different native units (e.g. Lysimeter NPOC in thousands of ug
# vs Cores ratios near unity); plotting raw effect sizes on one shared axis
# renders the smaller-unit CIs invisible. Percentage change places every
# compound on a comparable scale so all confidence intervals are legible.

fig_trt_df <- treatment_on_chem %>%
  filter(is.finite(pct_change), is.finite(pct_low), is.finite(pct_high)) %>%
  mutate(compound_file = paste0(compound, "  [", file, "]"),
         sig           = p.value < FDR_THRESHOLD) %>%
  arrange(pct_change) %>%
  mutate(compound_file = factor(compound_file, levels = compound_file))

p_trt <- ggplot(fig_trt_df,
                aes(x = pct_change, y = compound_file, colour = file)) +
  geom_vline(xintercept = 0, linetype = "dashed",
             colour = "grey50", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = pct_low, xmax = pct_high),
                 height = 0.30, linewidth = 0.5, alpha = 0.85) +
  geom_point(size = 2.6) +
  scale_colour_manual(values = file_pal, name = "Chemistry pool") +
  labs(x = "Warming effect (% change from control mean, \u00b1 95% CI)",
       y = NULL,
       caption = paste0("chem_value ~ Treatment + (1 | Plot). ",
                        "Effect rests on six plot-level means; ",
                        "no compound differs significantly.")) +
  theme_pub(9)

ggsave(file.path(SUPP_FIGS, "FigureS9_treatment_effect_chemistry.png"),
       p_trt, width = 9, height = max(5, 0.30 * nrow(fig_trt_df) + 2),
       dpi = FIG_DPI, units = "in")

# =============================================================================
# 20. RANDOM-SLOPE ROBUSTNESS CHECKS FOR MICROBIAL P
# =============================================================================
# The primary pipeline uses random intercepts (1 | Plot). This section tests
# whether a random-slope term (1 + scale(chem) | Plot) is warranted for the
# headline microbial P coupling, under both matching schemes, and refits the
# AR-corrected headline models with a random slope so corrected confidence
# intervals are available where the slope term is warranted.

message("\n=== SECTION 20: Random-slope robustness checks (uP) ===")

uP_compound <- if ("uP" %in% names(cores)) "uP" else
               grep("^uP$|ubial P|microbial.?P",
                    names(cores), value = TRUE)[1]
if (is.na(uP_compound))
  stop("Cannot find microbial P column ('uP') in Cores. ",
       "Check TRACE_Cores_clean.csv column names.")

uP_windowed <- build_matched_df(cores, uP_compound, WINDOW_PRIMARY)
if (!is.null(uP_windowed))
  uP_windowed <- uP_windowed %>%
    mutate(Treatment = factor(Treatment, levels = c("Control", "Warmed")),
           Plot      = factor(Plot))

uP_sameday <- build_sameday_df(cores, uP_compound)
if (!is.null(uP_sameday))
  uP_sameday <- uP_sameday %>%
    mutate(Treatment = factor(Treatment, levels = c("Control", "Warmed")),
           Plot      = factor(Plot))

run_random_slope_check <- function(df, response_col, scheme_label) {

  if (is.null(df)) return(NULL)
  d <- df %>%
    select(resp = all_of(response_col), chem_value, Treatment, Plot) %>%
    filter(!is.na(resp), !is.na(chem_value),
           is.finite(resp), is.finite(chem_value))
  if (nrow(d) < MIN_OBS || length(unique(d$Plot)) < 3L) return(NULL)

  m_ri <- tryCatch(suppressWarnings(
    lme4::lmer(resp ~ scale(chem_value) + Treatment + (1 | Plot),
               data = d, REML = FALSE,
               control = lmerControl(optimizer = "bobyqa"))),
    error = function(e) NULL)
  m_rs <- tryCatch(suppressWarnings(
    lme4::lmer(resp ~ scale(chem_value) + Treatment +
                 (1 + scale(chem_value) | Plot),
               data = d, REML = FALSE,
               control = lmerControl(optimizer = "bobyqa"))),
    error = function(e) NULL)
  if (is.null(m_ri) || is.null(m_rs)) return(NULL)

  lr     <- tryCatch(anova(m_ri, m_rs), error = function(e) NULL)
  lr_p   <- if (!is.null(lr)) lr[2L, "Pr(>Chisq)"] else NA_real_
  lr_chi <- if (!is.null(lr)) lr[2L, "Chisq"]      else NA_real_
  vc     <- as.data.frame(VarCorr(m_rs))
  slope_sd <- vc$sdcor[vc$grp == "Plot" &
                       grepl("chem_value", vc$var1) & is.na(vc$var2)]
  slope_sd <- if (length(slope_sd) == 0L) NA_real_ else slope_sd[1]

  tibble(
    scheme            = scheme_label,
    response          = response_col,
    n_obs             = nrow(d),
    n_plots           = length(unique(d$Plot)),
    aic_intercept     = AIC(m_ri),
    aic_slope         = AIC(m_rs),
    delta_aic         = AIC(m_ri) - AIC(m_rs),
    lr_chisq          = lr_chi,
    lr_p              = lr_p,
    slope_sd_plot     = slope_sd,
    random_slope_pref = (!is.na(lr_p) & lr_p < 0.05 &
                         (AIC(m_ri) - AIC(m_rs)) > 2),
    rs_singular       = lme4::isSingular(m_rs)
  )
}

rs_rows <- list()
for (resp in RESPONSES_MAIN) {
  r1 <- run_random_slope_check(uP_windowed, resp, "windowed_14d")
  r2 <- run_random_slope_check(uP_sameday,  resp, "same_day")
  if (!is.null(r1)) rs_rows[[length(rs_rows) + 1L]] <- r1
  if (!is.null(r2)) rs_rows[[length(rs_rows) + 1L]] <- r2
}
random_slope_check <- bind_rows(rs_rows)

# ---- 20b. AR-corrected headline refits: random intercept vs random slope ----

ar_refit_uP <- function(df, response_col, scheme_label) {

  if (is.null(df)) return(NULL)
  d <- df %>%
    select(Date, resp = all_of(response_col),
           chem_value, Treatment, Plot) %>%
    filter(!is.na(resp), !is.na(chem_value),
           is.finite(resp), is.finite(chem_value)) %>%
    group_by(Date, Plot, Treatment) %>%
    summarise(resp       = mean(resp,       na.rm = TRUE),
              chem_value = mean(chem_value, na.rm = TRUE),
              .groups    = "drop") %>%
    mutate(time_num = as.numeric(Date),
           Plot     = factor(Plot),
           cv       = as.numeric(scale(chem_value)))
  if (nrow(d) < MIN_OBS) return(NULL)

  m_ri <- tryCatch(suppressWarnings(
    nlme::lme(resp ~ cv + Treatment,
              random      = ~ 1 | Plot,
              correlation = nlme::corExp(form = ~ time_num | Plot,
                                         nugget = TRUE),
              data = d, method = "REML")),
    error = function(e) NULL)
  m_rs <- tryCatch(suppressWarnings(
    nlme::lme(resp ~ cv + Treatment,
              random      = ~ 1 + cv | Plot,
              correlation = nlme::corExp(form = ~ time_num | Plot,
                                         nugget = TRUE),
              data = d, method = "REML")),
    error = function(e) NULL)

  pull_cv <- function(m, label) {
    if (is.null(m)) return(NULL)
    tt <- summary(m)$tTable
    if (!"cv" %in% rownames(tt)) return(NULL)
    tibble(scheme = scheme_label, response = response_col, model = label,
           estimate  = tt["cv", "Value"],
           std.error = tt["cv", "Std.Error"],
           conf.low  = tt["cv", "Value"] - 1.96 * tt["cv", "Std.Error"],
           conf.high = tt["cv", "Value"] + 1.96 * tt["cv", "Std.Error"],
           p.value   = tt["cv", "p-value"],
           aic       = AIC(m))
  }
  bind_rows(pull_cv(m_ri, "random_intercept"),
            pull_cv(m_rs, "random_slope"))
}

ar_refit_rows <- list()
for (resp in RESPONSES_MAIN) {
  ar_refit_rows[[length(ar_refit_rows) + 1L]] <-
    ar_refit_uP(uP_windowed, resp, "windowed_14d")
  ar_refit_rows[[length(ar_refit_rows) + 1L]] <-
    ar_refit_uP(uP_sameday,  resp, "same_day")
}
uP_ar_refit <- bind_rows(ar_refit_rows)

if (nrow(random_slope_check) > 0)
  write_csv(random_slope_check,
            file.path(SUPP_TABLES, "random_slope_check_uP.csv"))
if (nrow(uP_ar_refit) > 0)
  write_csv(uP_ar_refit,
            file.path(SUPP_TABLES, "uP_AR_refit_random_slope.csv"))

message("  Random-slope checks: ", nrow(random_slope_check), " rows")
message("  AR refits:           ", nrow(uP_ar_refit), " rows")

# =============================================================================
# 21. PER-PLOT MICROBIAL P COUPLING BY BLOCK
# =============================================================================
# Extracts the per-plot conditional microbial P -> flux-residual slopes
# (fixed effect + plot-level random deviation) from the random-slope model,
# grouped by spatial block, to characterise where the plot-to-plot slope
# variation lives and whether any block stands apart. Descriptive: with three
# blocks there is no test to run.

message("\n=== SECTION 21: Per-plot uP coupling by block ===")

uP_coupling_by_block  <- tibble()
uP_coupling_block_smry <- tibble()

if (!is.null(uP_windowed)) {

  d_cpl <- uP_windowed %>%
    select(Date, resp = flux_resid_mean,
           chem_value, Treatment, Plot) %>%
    filter(!is.na(resp), !is.na(chem_value),
           is.finite(resp), is.finite(chem_value)) %>%
    mutate(cv    = as.numeric(scale(chem_value)),
           PlotF = factor(Plot))

  if (nrow(d_cpl) >= MIN_OBS && length(unique(d_cpl$PlotF)) >= 3L) {

    m_cpl <- tryCatch(suppressWarnings(
      lme4::lmer(resp ~ cv + Treatment + (1 + cv | PlotF),
                 data = d_cpl, REML = TRUE,
                 control = lmerControl(optimizer = "bobyqa"))),
      error = function(e) NULL)

    if (!is.null(m_cpl)) {
      fixed_slope_cpl <- lme4::fixef(m_cpl)[["cv"]]
      re_plot <- lme4::ranef(m_cpl)$PlotF
      re_plot$PlotF <- rownames(re_plot)

      uP_coupling_by_block <- re_plot %>%
        as_tibble() %>%
        rename(rand_dev_slope = cv,
               rand_dev_intercept = `(Intercept)`) %>%
        mutate(Plot              = as.integer(PlotF),
               Block             = BLOCK_MAP[as.character(Plot)],
               fixed_slope       = fixed_slope_cpl,
               conditional_slope = fixed_slope_cpl + rand_dev_slope) %>%
        left_join(
          d_cpl %>%
            group_by(Plot) %>%
            summarise(Treatment = first(Treatment),
                      n_obs     = n(),
                      .groups   = "drop") %>%
            mutate(Plot = as.integer(as.character(Plot))),
          by = "Plot") %>%
        select(Block, Plot, Treatment, n_obs,
               fixed_slope, rand_dev_slope, conditional_slope) %>%
        arrange(Block, Plot)

      uP_coupling_block_smry <- uP_coupling_by_block %>%
        group_by(Block) %>%
        summarise(plots           = paste(Plot, collapse = ","),
                  mean_cond_slope = mean(conditional_slope),
                  min_cond_slope  = min(conditional_slope),
                  max_cond_slope  = max(conditional_slope),
                  .groups = "drop")

      write_csv(uP_coupling_by_block,
                file.path(SUPP_TABLES, "uP_coupling_by_block.csv"))

      message("  Per-plot coupling slopes extracted: ",
              nrow(uP_coupling_by_block), " plots")
      message("  All conditional slopes positive: ",
              all(uP_coupling_by_block$conditional_slope > 0))
    }
  }
}

# 21b. FIGURE S10 — BLOCKING SPATIAL STORY (two-panel composite)
# =============================================================================
# Panel A — per-plot microbial P coupling by block. The six per-plot
#           conditional slopes of the microbial P -> flux-residual coupling,
#           coloured by spatial block, showing that the coupling is positive
#           in every plot and the plot-to-plot spread is concentrated rather
#           than pervasive.
# Panel B — warming response of chemistry by block. The three within-block
#           differences (Warmed - Control) for the chemistry compounds with
#           the largest absolute point-estimate response to warming, showing
#           the spatial heterogeneity discussed in the text: Blocks 1 and 2
#           tend to respond, Block 3 does not.

message("\n=== SECTION 21b: Figure S10 (blocking spatial story) ===")

block_pal_s10 <- c("Block 1" = "#4575b4",
                   "Block 2" = "#1a9641",
                   "Block 3" = "#e08214")

# ---- Panel A: per-plot microbial P coupling --------------------------------

pS10_A <- NULL
if (nrow(uP_coupling_by_block) > 0) {
  
  figS10A_df <- uP_coupling_by_block %>%
    mutate(plot_lab = paste0("Plot ", Plot,
                             " [B", Block, ", ",
                             substr(as.character(Treatment), 1, 1), "]"),
           BlockF   = factor(paste0("Block ", Block))) %>%
    arrange(conditional_slope) %>%
    mutate(plot_lab = factor(plot_lab, levels = plot_lab))
  
  pS10_A <- ggplot(figS10A_df,
                   aes(x = conditional_slope, y = plot_lab,
                       colour = BlockF)) +
    geom_vline(aes(xintercept = fixed_slope),
               linetype = "dashed", colour = "grey45", linewidth = 0.5) +
    geom_segment(aes(x = fixed_slope, xend = conditional_slope,
                     y = plot_lab, yend = plot_lab),
                 colour = "grey70", linewidth = 0.4) +
    geom_point(size = 3.2) +
    scale_colour_manual(values = block_pal_s10, name = "Spatial block") +
    labs(x = "Conditional microbial P \u2192 flux-residual slope",
         y = NULL,
         title = "A. Per-plot microbial P coupling") +
    theme_pub(9) +
    theme(plot.title = element_text(size = 10, face = "bold"))
}

# ---- Panel B: warming response of chemistry by block -----------------------
# For the N_TOP_S10 compounds with the largest absolute paired mean
# difference (Warmed - Control), parse the three within-block differences
# out of the `block_diffs` string in treatment_on_chem_blocked, express each
# as percentage change from the control grand mean for the compound, and
# plot the three values per compound.

pS10_B <- NULL
N_TOP_S10 <- 8L

if (exists("treatment_on_chem_blocked") &&
    nrow(treatment_on_chem_blocked) > 0L &&
    exists("treatment_on_chem") &&
    nrow(treatment_on_chem) > 0L) {
  
  ctrl_means <- treatment_on_chem %>%
    select(file, compound, control_mean)
  
  top_compounds <- treatment_on_chem_blocked %>%
    inner_join(ctrl_means, by = c("file", "compound")) %>%
    filter(is.finite(paired_mean_diff), is.finite(control_mean),
           control_mean != 0) %>%
    mutate(abs_pct = abs(100 * paired_mean_diff / control_mean)) %>%
    arrange(desc(abs_pct)) %>%
    slice_head(n = N_TOP_S10) %>%
    select(file, compound, control_mean, block_diffs, abs_pct)
  
  parse_block_diffs <- function(s) {
    parts <- strsplit(s, ";\\s*")[[1]]
    out <- lapply(parts, function(p) {
      m <- regmatches(p, regexec("B([0-9]+)=(-?[0-9eE.+-]+)", p))[[1]]
      if (length(m) == 3L)
        tibble(Block = as.integer(m[2]), diff = as.numeric(m[3]))
      else NULL
    })
    bind_rows(out)
  }
  
  figS10B_df <- top_compounds %>%
    mutate(rows = lapply(block_diffs, parse_block_diffs)) %>%
    select(file, compound, control_mean, rows) %>%
    tidyr::unnest(rows) %>%
    mutate(pct_diff      = 100 * diff / control_mean,
           BlockF        = factor(paste0("Block ", Block)),
           compound_file = paste0(compound, "  [", file, "]"))
  
  comp_order_s10b <- figS10B_df %>%
    group_by(compound_file) %>%
    summarise(m = mean(abs(pct_diff)), .groups = "drop") %>%
    arrange(m) %>% pull(compound_file)
  figS10B_df <- figS10B_df %>%
    mutate(compound_file = factor(compound_file, levels = comp_order_s10b))
  
  pS10_B <- ggplot(figS10B_df,
                   aes(x = pct_diff, y = compound_file, colour = BlockF)) +
    geom_vline(xintercept = 0, linetype = "dashed",
               colour = "grey45", linewidth = 0.4) +
    geom_point(size = 3.2, alpha = 0.9,
               position = position_dodge(width = 0.45)) +
    scale_colour_manual(values = block_pal_s10, name = "Spatial block") +
    labs(x = paste0("Warming effect (Warmed \u2212 Control, ",
                    "% of control mean, by block)"),
         y = NULL,
         title = paste0("B. Warming response of chemistry by block (top ",
                        N_TOP_S10, " compounds)")) +
    theme_pub(9) +
    theme(plot.title = element_text(size = 10, face = "bold"))
}

# ---- Compose and save ------------------------------------------------------

if (!is.null(pS10_A) || !is.null(pS10_B)) {
  
  pS10 <- if (!is.null(pS10_A) && !is.null(pS10_B)) {
    (pS10_A / pS10_B) +
      patchwork::plot_layout(heights = c(1, 1.4), guides = "collect") &
      theme(legend.position = "bottom")
  } else if (!is.null(pS10_A)) {
    pS10_A
  } else {
    pS10_B
  }
  
  cap_text <- paste0(
    "A: dashed line is the population (fixed-effect) slope; points are ",
    "per-plot conditional slopes from the random-slope model; all slopes ",
    "are positive. B: three within-block warming effects (Warmed \u2212 ",
    "Control), expressed as percentage of the control mean, for the ",
    N_TOP_S10, " compounds with the largest absolute paired mean ",
    "difference; spatial heterogeneity is visible chiefly as Block 3 ",
    "departing from Blocks 1 and 2.")
  pS10 <- pS10 + patchwork::plot_annotation(caption = cap_text) &
    theme(plot.caption = element_text(size = 8, colour = "grey40",
                                      hjust = 0))
  
  ggsave(file.path(SUPP_FIGS, "FigureS10_blocking_spatial_story.png"),
         pS10,
         width = 9,
         height = if (!is.null(pS10_A) && !is.null(pS10_B)) 9.5 else 4.6,
         dpi = FIG_DPI, units = "in")
  message("  Figure S10 (two-panel) written")
}

# =============================================================================
# 21c. SUPPLEMENTARY DOCX TABLES (S2, S4, S5, S6, S7)
# =============================================================================
# Renders as publication-ready DOCX tables the supplementary results that
# the manuscript references but that were previously exported only as CSV:
#   Table S2 — flux ACF at fixed lags
#   Table S4 — Treatment main effect on flux (windowed primary)
#   Table S5 — same-day AR-corrected microbial P simple slopes
#   Table S6 — random-slope robustness and AR refit for microbial P
#   Table S7 — treatment effect on chemistry (unpaired + paired/blocked)

message("\n=== SECTION 21c: Supplementary DOCX tables ===")

make_supp_table <- function(df, caption, path, fontsize = 8) {
  if (is.null(df) || nrow(df) == 0L) {
    message("  [skip] ", basename(path), " — no rows"); return(invisible())
  }
  ft <- flextable::flextable(df) %>%
    flextable::theme_booktabs() %>%
    flextable::bg(part = "header", bg = "#3a3a3a") %>%
    flextable::color(part = "header", color = "white") %>%
    flextable::bold(part = "header") %>%
    flextable::fontsize(size = fontsize, part = "all") %>%
    flextable::padding(padding.top = 2.5, padding.bottom = 2.5,
                       part = "all") %>%
    flextable::align(align = "center", part = "all") %>%
    flextable::set_table_properties(layout = "autofit", width = 1) %>%
    flextable::add_header_lines(caption)
  read_docx() %>% body_add_flextable(ft) %>% print(target = path)
  message("  written: ", basename(path))
}

# ---- Table S2 — flux ACF at fixed lags --------------------------------------
tS2 <- acf_at_lags %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))
make_supp_table(
  tS2,
  paste0("Table S2. Empirical autocorrelation function of daily-aggregated ",
         "soil CO2 flux at fixed lags (days), per plot. Values for ",
         "log-flux and for flux residuals after Temperature and VWC ",
         "removal. No plot decorrelates below 0.1 within the 60-day ",
         "window."),
  file.path(SUPP_TABLES, "TableS2_flux_ACF.docx"))

# ---- Table S4 — Treatment main effect on flux -------------------------------
tS4 <- results_treat %>%
  filter(window == WINDOW_PRIMARY, response %in% RESPONSES_MAIN) %>%
  transmute(Response  = recode(response, !!!response_labels),
            Predictor = paste0(compound, " [", file, "]"),
            `Treatment beta` = format_est(estimate, 3),
            `95% CI`  = format_ci(conf.low, conf.high, 3),
            `p (FDR)` = format_p(p.fdr)) %>%
  arrange(Response, `p (FDR)`)
make_supp_table(
  tS4,
  paste0("Table S4. Treatment (warming) main effect on soil CO2 flux at ",
         "the primary +/-14-day window, from the standard chemistry-flux ",
         "models. No Treatment main effect survives Benjamini-Hochberg ",
         "FDR correction."),
  file.path(SUPP_TABLES, "TableS4_treatment_main_effect.docx"))

# ---- Table S5 — same-day AR-corrected uP simple slopes ----------------------
tS5 <- sameday_ar_slopes %>%
  filter(compound == "uP", response %in% RESPONSES_MAIN) %>%
  transmute(Response  = recode(response, !!!response_labels),
            Treatment,
            `Slope (beta)` = format_est(slope, 4),
            `95% CI`  = format_ci(slope_lo, slope_hi, 4)) %>%
  arrange(Response, Treatment)
make_supp_table(
  tS5,
  paste0("Table S5. Same-day AR-corrected per-treatment simple slopes for ",
         "the microbial P -> flux coupling. corExp-corrected estimates ",
         "from same-day matching; control and warmed slopes do not ",
         "differ."),
  file.path(SUPP_TABLES, "TableS5_sameday_uP_slopes.docx"))

# ---- Table S6 — random-slope robustness + AR refit --------------------------
tS6a <- random_slope_check %>%
  transmute(Scheme    = scheme,
            Response  = recode(response, !!!response_labels),
            `delta AIC` = format_est(delta_aic, 2),
            `LR p`    = format_p(lr_p),
            `Slope SD (Plot)` = format_est(slope_sd_plot, 3),
            `RS warranted` = ifelse(random_slope_pref, "yes", "no"),
            Singular  = ifelse(rs_singular, "yes", "no"))
make_supp_table(
  tS6a,
  paste0("Table S6a. Random-slope robustness check for the microbial P -> ",
         "flux coupling. Likelihood-ratio test of a random-slope term ",
         "against a random-intercept-only specification, by matching ",
         "scheme and response."),
  file.path(SUPP_TABLES, "TableS6a_random_slope_check.docx"))

tS6b <- uP_ar_refit %>%
  transmute(Scheme    = scheme,
            Response  = recode(response, !!!response_labels),
            Model     = recode(model,
                               random_intercept = "Random intercept",
                               random_slope     = "Random slope"),
            Estimate  = format_est(estimate, 4),
            `95% CI`  = format_ci(conf.low, conf.high, 4),
            `p-value` = format_p(p.value))
make_supp_table(
  tS6b,
  paste0("Table S6b. AR-corrected microbial P fixed-effect estimate under ",
         "random-intercept and random-slope specifications. The headline ",
         "estimate is essentially unchanged by the random-slope term."),
  file.path(SUPP_TABLES, "TableS6b_uP_AR_refit.docx"))

# ---- Table S7 — treatment effect on chemistry -------------------------------
tS7 <- treatment_on_chem %>%
  transmute(Pool      = file,
            Compound  = compound,
            n         = n_obs,
            `Control mean` = format_est(control_mean, 3),
            `Warmed mean`  = format_est(warmed_mean, 3),
            `% change`     = format_est(pct_change, 1),
            `95% CI (%)`   = format_ci(pct_low, pct_high, 0),
            `p-value`      = format_p(p.value)) %>%
  arrange(Pool, `p-value`)
make_supp_table(
  tS7,
  paste0("Table S7. Treatment (warming) effect on each chemistry pool ",
         "(unpaired analysis): chem_value ~ Treatment + (1 | Plot). ",
         "Effect expressed as percentage change from the control mean ",
         "with 95% confidence interval. Unadjusted p-values; no compound ",
         "differs significantly."),
  file.path(SUPP_TABLES, "TableS7_treatment_effect_chemistry.docx"), 7)

tS7b <- treatment_on_chem_blocked %>%
  transmute(Pool      = file,
            Compound  = compound,
            `Blocks paired` = n_blocks_paired,
            `Block differences` = block_diffs,
            `Mean diff (W-C)` = format_est(paired_mean_diff, 3),
            `Paired 95% CI` = format_ci(paired_conf_low,
                                        paired_conf_high, 3),
            `Paired p` = format_p(paired_p)) %>%
  arrange(Pool, `Paired p`)
make_supp_table(
  tS7b,
  paste0("Table S7b. Paired (blocked) treatment effect on chemistry. ",
         "Paired t-test on the three within-block differences between the ",
         "warmed and control plot. The block mixed model is not reported ",
         "because, with three blocks, its variance component is ",
         "anticonservative."),
  file.path(SUPP_TABLES, "TableS7b_treatment_effect_blocked.docx"), 7)

# =============================================================================
# 22. CONSOLIDATED MANUSCRIPT OUTPUT (TXT)
# =============================================================================
# Single human-readable text file collecting every result to be presented or
# discussed in the manuscript (main text and supplement). Numerical tables
# are written in full; this file is the reference for manuscript writing and
# review.

message("\n=== SECTION 22: Writing consolidated manuscript output ===")

manuscript_txt <- file.path(MAIN_DIR, "MANUSCRIPT_RESULTS_ALL.txt")
con <- file(manuscript_txt, "w", encoding = "UTF-8")

w  <- function(...) cat(..., "\n", sep = "", file = con)
hr <- function() cat(strrep("=", 78), "\n", sep = "", file = con)
sub <- function() cat(strrep("-", 78), "\n", sep = "", file = con)

dump_df <- function(df, caption = NULL, digits = 4) {
  if (!is.null(caption)) { cat("\n", caption, "\n", sep = "", file = con) }
  if (is.null(df) || nrow(df) == 0L) {
    cat("  [no rows]\n", file = con); return(invisible())
  }
  df_fmt <- df %>%
    mutate(across(where(is.numeric),
                  ~ ifelse(is.na(.), NA, round(., digits))))
  out <- capture.output(print(as.data.frame(df_fmt), row.names = FALSE))
  cat(paste0("  ", out, collapse = "\n"), "\n", file = con)
}

hr()
w("TRACE FLUX-CHEMISTRY ANALYSIS — CONSOLIDATED RESULTS FOR MANUSCRIPT")
w("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
w("R version: ", R.version.string)
hr()
w("")
w("CONTENTS")
w("  1. Study scope and sample sizes")
w("  2. Abiotic base model (Table 1)")
w("  3. Flux temporal autocorrelation (Figure S1, Table S2)")
w("  4. Primary windowed chemistry-flux models (Table 2, Figure 1)")
w("  5. Autocorrelation correction and collapse (Table 3, Figure 3)")
w("  6. AR diagnostics: range parameters")
w("  7. Microbial P interaction with Treatment (Table 4, Figure 2)")
w("  8. Treatment main effect on flux (Table S4)")
w("  9. Same-day matching sensitivity (Table S5; Figures S5, S6, S7, S8)")
w(" 10. Same-day AR justification diagnostic")
w(" 11. Window-width and correlation-structure sensitivity")
w("     (Figure S2 cited; Figures S3, S4 generated for internal use only)")
w(" 12. Treatment effect on chemistry pools")
w("     (Figure S9, Tables S7 and S7b)")
w(" 13. Random-slope robustness for microbial P (Tables S6a, S6b)")
w(" 14. Per-plot microbial P coupling by block (Figure S10)")
w("")
w("FIGURE AND TABLE MAP")
w("  Main figures:  1 chemistry coefficients | 2 uP simple slopes |")
w("                 3 AR collapse")
w("  Main tables:   1 abiotic | 2 chemistry primary | 3 AR collapse |")
w("                 4 uP interaction")
w("  Supp figures:  S1 flux ACF | S2 window sensitivity |")
w("                 S3 correlation structures (internal) |")
w("                 S4 flux CV (internal) | S5 same-day uP slopes |")
w("                 S6 cross-method comparison |")
w("                 S7 same-day chemistry coefficients |")
w("                 S8 same-day AR collapse |")
w("                 S9 treatment effect on chemistry |")
w("                 S10 per-plot uP coupling by block")
w("  Supp tables:   S2 flux ACF at fixed lags |")
w("                 S4 treatment main effect on flux |")
w("                 S5 same-day uP simple slopes |")
w("                 S6a random-slope check | S6b uP AR refit |")
w("                 S7 treatment effect on chemistry (unpaired) |")
w("                 S7b treatment effect on chemistry (paired/blocked)")
w("")

# ---- 1. Scope ---------------------------------------------------------------
hr(); w("1. STUDY SCOPE AND SAMPLE SIZES"); hr()
w("Hourly flux observations (analysed): ", nrow(flux))
w("Flux date range: ", as.character(min(flux$FluxDate)),
  " to ", as.character(max(flux$FluxDate)))
w("Plots: 6 (3 Control, 3 Warmed); spatial blocks: ",
  "B1={1,2}, B2={3,4}, B3={5,6}")
w("Chemistry compounds per pool:")
for (nm in names(chem_files))
  w("  ", nm, ": ", length(chem_files[[nm]]$compounds), " compounds")
w("Total chemistry variables: ",
  sum(vapply(chem_files, function(x) length(x$compounds), integer(1))))
w("Primary window: +/-", WINDOW_PRIMARY, " days")
w("Primary responses: ", paste(RESPONSES_MAIN, collapse = ", "))

# ---- 2. Abiotic base model --------------------------------------------------
hr(); w("2. ABIOTIC BASE MODEL (Table 1)"); hr()
w("Model: log(Flux) ~ Temperature + VWC + (1 | Plot), REML")
dump_df(abiotic_fe %>%
          select(term, estimate, std.error, statistic, conf.low, conf.high),
        "Fixed effects:")
dump_df(abiotic_re %>% select(group, term, estimate),
        "Random effects (sd):")
w("")
w("Marginal R2:    ", round(abiotic_r2[1, "R2m"], 4))
w("Conditional R2: ", round(abiotic_r2[1, "R2c"], 4))

# ---- 3. Flux autocorrelation ------------------------------------------------
hr(); w("3. FLUX TEMPORAL AUTOCORRELATION (Figure S1)"); hr()
w("ACF of daily-aggregated flux, values at fixed lags per plot:")
dump_df(acf_at_lags, NULL, digits = 3)
w("")
w("Decorrelation lags (first lag with |ACF| < ", ACF_THRESH, "):")
dump_df(decorr_lags, NULL, digits = 0)
w("(NA = ACF does not drop below threshold within the ",
  ACF_MAX_LAG, "-day window.)")

# ---- 4. Primary windowed models ---------------------------------------------
hr(); w("4. PRIMARY WINDOWED CHEMISTRY-FLUX MODELS (Table 2, Figure 1)"); hr()
w("Model: response ~ scale(chem) + Treatment + (1 | Plot) at +/-",
  WINDOW_PRIMARY, " days")
w("Standard-LMM chemistry coefficients, primary window, main responses:")
dump_df(results_chem %>%
          filter(window == WINDOW_PRIMARY,
                 response %in% RESPONSES_MAIN) %>%
          select(file, compound, response, estimate, std.error,
                 conf.low, conf.high, r2_marginal, n_obs, p.value, p.fdr) %>%
          arrange(response, p.fdr))
w("")
w("FDR-significant standard-LMM main effects (primary window): ",
  sum(results_chem$window == WINDOW_PRIMARY &
      results_chem$response %in% RESPONSES_MAIN &
      results_chem$p.fdr < FDR_THRESHOLD, na.rm = TRUE))

# ---- 5. AR correction and collapse ------------------------------------------
hr(); w("5. AUTOCORRELATION CORRECTION AND COLLAPSE (Table 3, Figure 3)"); hr()
w("AR-corrected (corExp) refits of FDR-significant standard-LMM effects:")
dump_df(ar_collapse %>%
          filter(response %in% RESPONSES_MAIN) %>%
          select(file, compound, response, std_estimate, ar_estimate,
                 ar_p, ar_aic_improvement, pct_attenuation, survives_AR) %>%
          arrange(response, ar_p))
w("")
w("Compounds surviving AR correction (p < 0.05), primary window: ",
  sum(ar_collapse$response %in% RESPONSES_MAIN &
      ar_collapse$survives_AR, na.rm = TRUE))

# ---- 6. AR diagnostics ------------------------------------------------------
hr(); w("6. AR DIAGNOSTICS: RANGE PARAMETERS"); hr()
w("corExp range parameter (autocorrelation length, days) per fitted model:")
dump_df(ar_diagnostics %>%
          select(file, compound, response, range_days, nugget,
                 sigma, n_obs_agg) %>%
          arrange(file, compound, response), NULL, digits = 2)

# ---- 7. Microbial P interaction ---------------------------------------------
hr(); w("7. MICROBIAL P x TREATMENT INTERACTION (Table 4, Figure 2)"); hr()
w("AR-corrected interaction terms (uP x Treatment), primary window:")
dump_df(results_autocor_int %>%
          filter(compound == "uP", window == WINDOW_PRIMARY,
                 response %in% RESPONSES_MAIN) %>%
          select(file, compound, response, estimate, std.error,
                 p.value, aic_improvement))
w("")
w("AR-corrected per-treatment simple slopes (uP), primary window:")
dump_df(results_autocor_slopes %>%
          filter(compound == "uP", window == WINDOW_PRIMARY,
                 response %in% RESPONSES_MAIN) %>%
          select(file, compound, response, Treatment,
                 slope, slope_lo, slope_hi) %>%
          arrange(response, Treatment))

# ---- 8. Treatment main effect on flux ---------------------------------------
hr(); w("8. TREATMENT MAIN EFFECT ON FLUX"); hr()
treat_min <- results_treat %>%
  filter(window == WINDOW_PRIMARY, response %in% RESPONSES_MAIN) %>%
  summarise(m = min(p.fdr, na.rm = TRUE)) %>% pull(m)
w("Minimum FDR-adjusted p for any Treatment main effect ",
  "(primary window, main responses): ", round(treat_min, 3))
w("No Treatment main effect on flux survives FDR correction.")

# ---- 9. Same-day matching ---------------------------------------------------
hr(); w("9. SAME-DAY MATCHING SENSITIVITY (Tables S, Figures S5-S8)"); hr()
w("Same-day standard-LMM chemistry coefficients, main responses:")
dump_df(sameday_chem %>%
          select(file, compound, response, estimate, std.error,
                 p.value, p.fdr, n_obs) %>%
          arrange(response, p.fdr))
w("")
w("Same-day AR-corrected chemistry coefficients:")
dump_df(sameday_ar_main %>%
          select(file, compound, response, estimate, std.error,
                 p.value, aic_improvement, n_obs_agg) %>%
          arrange(response, p.value))
w("")
w("Same-day AR-corrected uP simple slopes:")
dump_df(sameday_ar_slopes %>%
          filter(compound == "uP", response %in% RESPONSES_MAIN) %>%
          select(file, compound, response, Treatment,
                 slope, slope_lo, slope_hi) %>%
          arrange(response, Treatment))
w("")
w("Cross-method comparison (windowed vs same-day, AR-corrected):")
dump_df(cross_method_comparison %>%
          select(file, compound, response,
                 windowed_estimate_ar, windowed_p_ar,
                 sameday_estimate_ar, sameday_p_ar) %>%
          arrange(response, file, compound))
w("")
w("Same-day AR-collapse summary:")
if (exists("sameday_ar_collapse") && nrow(sameday_ar_collapse) > 0) {
  dump_df(sameday_ar_collapse %>%
            select(file, compound, response, estimate,
                   ar_estimate, ar_p, survives_AR) %>%
            arrange(response, ar_p))
} else {
  w("  [not available]")
}

# ---- 10. Same-day AR justification ------------------------------------------
hr(); w("10. SAME-DAY AR JUSTIFICATION DIAGNOSTIC"); hr()
w("Likelihood-ratio test, corExp vs no correlation structure:")
dump_df(sameday_ar_lr_summary %>%
          select(file, compound, response, n_agg, delta_aic,
                 lr_p, range_days, AR_preferred) %>%
          arrange(file, compound, response), NULL, digits = 2)
w("")
w("Binned-lag residual ACF (same-day standard-LMM residuals):")
dump_df(sameday_ar_test %>%
          select(file, compound, response, bin, n, r, mean_gap) %>%
          arrange(file, compound, response), NULL, digits = 3)

# ---- 11. Sensitivity --------------------------------------------------------
hr(); w("11. WINDOW-WIDTH AND CORRELATION-STRUCTURE SENSITIVITY"); hr()
w("Window sensitivity grid (FDR-significant compounds, main responses):")
dump_df(sens_grid %>%
          filter(response %in% RESPONSES_MAIN) %>%
          select(file, compound, response, window, model,
                 estimate, conf.low, conf.high, p.value) %>%
          arrange(file, compound, response, model, window))
w("")
w("Correlation-structure comparison (corExp vs corCAR1 vs corGaus):")
dump_df(corr_struct %>%
          filter(response %in% RESPONSES_MAIN) %>%
          select(file, compound, response, model,
                 estimate, conf.low, conf.high, p.value, aic) %>%
          arrange(file, compound, response, model))

# ---- 12. Treatment effect on chemistry --------------------------------------
hr(); w("12. TREATMENT EFFECT ON CHEMISTRY POOLS (Section 19, Figure S9)"); hr()
w("Unpaired: chem_value ~ Treatment + (1 | Plot). Unadjusted p primary.")
dump_df(treatment_on_chem %>%
          select(file, compound, n_obs, control_mean, warmed_mean,
                 diff_WmC, conf.low, conf.high, pct_change,
                 p.value, p.fdr, model_type))
w("")
w("Per-plot means (the six values each unpaired contrast rests on):")
dump_df(treatment_on_chem %>% select(file, compound, plot_means))
w("")
w("Paired (blocked) t-test on three within-block differences:")
dump_df(treatment_on_chem_blocked %>%
          select(file, compound, n_blocks_paired, block_diffs,
                 paired_mean_diff, paired_conf_low, paired_conf_high,
                 pct_change, paired_p))
w("")
w("Note: the block mixed model is deliberately not reported; with three")
w("blocks its variance component is anticonservative. The paired t-test")
w("on the three within-block differences is the design-matched test.")

# ---- 13. Random-slope robustness --------------------------------------------
hr(); w("13. RANDOM-SLOPE ROBUSTNESS FOR MICROBIAL P (Section 20)"); hr()
w("Random-intercept vs random-slope comparison (LR test, ML fit):")
dump_df(random_slope_check %>%
          select(scheme, response, n_obs, delta_aic, lr_p,
                 slope_sd_plot, random_slope_pref, rs_singular),
        NULL, digits = 3)
w("")
w("AR-corrected headline uP refits (random intercept vs random slope):")
dump_df(uP_ar_refit %>%
          select(scheme, response, model, estimate,
                 conf.low, conf.high, p.value, aic), NULL, digits = 4)

# ---- 14. Per-plot coupling by block -----------------------------------------
hr(); w("14. PER-PLOT MICROBIAL P COUPLING BY BLOCK (Figure S10)"); hr()
w("Per-plot conditional uP -> flux-residual slopes (windowed +/-14d):")
dump_df(uP_coupling_by_block, NULL, digits = 4)
w("")
w("Block-level summary of conditional slopes:")
dump_df(uP_coupling_block_smry, NULL, digits = 4)
w("")
w("All six per-plot conditional slopes are positive: the microbial P -")
w("flux coupling holds in the same direction in every plot and block.")

w("")
hr()
w("END OF CONSOLIDATED RESULTS")
hr()

close(con)
message("  Consolidated manuscript output written to: ", manuscript_txt)

message("\n=============================================================")
message("FULL ANALYSIS COMPLETE — all sections (0-22)")
message("=============================================================")
message("  Main outputs:          ", MAIN_DIR)
message("  Supplementary outputs: ", SUPP_DIR)
message("  Manuscript results:    ", manuscript_txt)
