# ============================================================================
# GLMM_tables.R
# ----------------------------------------------------------------------------
# Publication tables (Times New Roman, black) from the GLMM outputs.
# Run after the primary GLMM script. Reads OUTPUT_GLMM_FINAL/tables/*.csv and
# the RDS of fitted objects.
#
# Model sourcing:
#   Overall warming effects        pooled models   (Treatment + Year_f + season)
#   Per-year / per-season effects  interaction models (Treatment * Year_f
#                                    + Treatment * season)
#   Q10 and VWC optimum            mechanistic modification model (centred
#                                    covariates temp_c, vwc_c)
#   Residual warming effect        mechanistic pathway model
#   Variability                    empirical per-cell log(SD) analysis
#
# Produces (in OUTPUT_GLMM_FINAL/):
#   Manuscript_Tables_Main.docx
#     Table 1  Warming effect on CO2 flux (overall, by season, by year)
#     Table 2  Warming effect on soil temperature and VWC (overall, season, year)
#     Table 3  Apparent Q10 and VWC optimum by treatment
#   Manuscript_Tables_Supplementary.docx
#     Table S1 Mechanistic model coefficients (pathway and modification)
#     Table S2 Model convergence and AR(1) summary (eight models)
#     Table S3 Per-plot flux summary
#     Table S4 Data coverage by pipeline, plot and year
#     Table S5 Treatment effect on day-to-day variability (empirical)
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(tibble)
  library(purrr); library(stringr); library(officer); library(flextable)
  library(glmmTMB)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

setwd("D:/USDA/TRACE_DM_APR_26/CO2")
out_dir <- "OUTPUT_GLMM_FINAL"
tab_dir <- file.path(out_dir, "tables")

obj <- readRDS(file.path(out_dir, "objects", "glmm_v3_objects.rds"))

read_safe <- function(p) if (file.exists(p)) read_csv(p, show_col_types = FALSE) else NULL

fmt_p <- function(p) ifelse(is.na(p), "-",
                     ifelse(p < 0.001, "<0.001", formatC(p, digits = 3, format = "fg")))
fmt_n <- function(x, d = 2) ifelse(is.na(x), "-", formatC(x, digits = d, format = "f"))
ci    <- function(lo, hi, d = 2) ifelse(is.na(lo) | is.na(hi), "-",
                                 sprintf("[%s, %s]", fmt_n(lo, d), fmt_n(hi, d)))

style_ft <- function(ft) {
  ft |>
    theme_booktabs() |>
    font(fontname = "Times New Roman", part = "all") |>
    fontsize(size = 9, part = "all") |>
    bold(part = "header") |>
    color(color = "black", part = "all") |>
    align(align = "left", part = "all") |>
    padding(padding = 3, part = "all") |>
    autofit()
}

# ============================================================================
# TABLE 1: Warming effect on CO2 flux (log estimates back-transformed)
#   Overall row from the pooled flux model; per-season and per-year rows from
#   the flux interaction model.
# ============================================================================
flux_marg   <- read_safe(file.path(tab_dir, "contrast_flux_marginal.csv"))
flux_season <- read_safe(file.path(tab_dir, "contrast_flux_by_season.csv"))
flux_year   <- read_safe(file.path(tab_dir, "contrast_flux_by_year.csv"))

add_ratio <- function(df) {
  if (is.null(df) || "status" %in% names(df)) return(NULL)
  df %>% mutate(
    ratio_warmed_control = exp(estimate),
    ratio_lower          = exp(lower.CL),
    ratio_upper          = exp(upper.CL),
    pct_change           = (ratio_warmed_control - 1) * 100
  )
}
flux_marg   <- add_ratio(flux_marg)
flux_season <- add_ratio(flux_season)
flux_year   <- add_ratio(flux_year)

row_overall <- if (!is.null(flux_marg)) tibble(
  Scope   = "Overall",
  Ratio   = fmt_n(flux_marg$ratio_warmed_control[1]),
  CI      = ci(flux_marg$ratio_lower[1], flux_marg$ratio_upper[1]),
  Pct     = fmt_n(flux_marg$pct_change[1], 1),
  Model_p = fmt_p(flux_marg$p.value[1])
) else tibble()
rows_season <- if (!is.null(flux_season)) flux_season |>
  transmute(Scope   = paste0("Season: ", season),
            Ratio   = fmt_n(ratio_warmed_control),
            CI      = ci(ratio_lower, ratio_upper),
            Pct     = fmt_n(pct_change, 1),
            Model_p = fmt_p(p.value)) else tibble()
rows_year <- if (!is.null(flux_year)) flux_year |>
  transmute(Scope   = paste0("Year: ", Year_f),
            Ratio   = fmt_n(ratio_warmed_control),
            CI      = ci(ratio_lower, ratio_upper),
            Pct     = fmt_n(pct_change, 1),
            Model_p = fmt_p(p.value)) else tibble()

tab1 <- bind_rows(row_overall, rows_season, rows_year) |>
  rename(`Warmed:control ratio` = Ratio, `95% CI` = CI, `% change` = Pct,
         `p (unadjusted)` = Model_p)
ft1 <- style_ft(flextable(tab1))

# ============================================================================
# TABLE 2: Warming effect on temperature and VWC (overall, per season, per year)
#   Overall from pooled models; per-year and per-season from interaction models.
#   VWC rescaled vwc_pct -> m^3 m^-3 (divide by 100).
# ============================================================================
temp_marg   <- read_safe(file.path(tab_dir, "contrast_temp_marginal.csv"))
temp_year   <- read_safe(file.path(tab_dir, "contrast_temp_by_year.csv"))
temp_season <- read_safe(file.path(tab_dir, "contrast_temp_by_season.csv"))
vwc_marg    <- read_safe(file.path(tab_dir, "contrast_vwc_marginal.csv"))
vwc_year    <- read_safe(file.path(tab_dir, "contrast_vwc_by_year.csv"))
vwc_season  <- read_safe(file.path(tab_dir, "contrast_vwc_by_season.csv"))

drop_status <- function(df) if (is.null(df) || "status" %in% names(df)) NULL else df
temp_marg <- drop_status(temp_marg); temp_year <- drop_status(temp_year)
temp_season <- drop_status(temp_season)
vwc_marg <- drop_status(vwc_marg); vwc_year <- drop_status(vwc_year)
vwc_season <- drop_status(vwc_season)

rescale_vwc <- function(df) {
  if (is.null(df)) return(NULL)
  df %>% mutate(across(any_of(c("estimate", "lower.CL", "upper.CL", "SE")),
                       ~ .x / 100))
}
vwc_marg   <- rescale_vwc(vwc_marg)
vwc_year   <- rescale_vwc(vwc_year)
vwc_season <- rescale_vwc(vwc_season)

build_rows <- function(marg, by_season, by_year, response, unit_d) {
  out <- tibble()
  if (!is.null(marg)) out <- bind_rows(out, tibble(
    Response = response, Scope = "Overall",
    Estimate = fmt_n(marg$estimate[1], unit_d),
    CI       = ci(marg$lower.CL[1], marg$upper.CL[1], unit_d),
    p        = fmt_p(marg$p.value[1])))
  if (!is.null(by_season)) out <- bind_rows(out, by_season |>
    transmute(Response = response,
              Scope    = paste0("Season: ", season),
              Estimate = fmt_n(estimate, unit_d),
              CI       = ci(lower.CL, upper.CL, unit_d),
              p        = fmt_p(p.value)))
  if (!is.null(by_year)) out <- bind_rows(out, by_year |>
    transmute(Response = response,
              Scope    = paste0("Year: ", Year_f),
              Estimate = fmt_n(estimate, unit_d),
              CI       = ci(lower.CL, upper.CL, unit_d),
              p        = fmt_p(p.value)))
  out
}
tab2 <- bind_rows(
  build_rows(temp_marg, temp_season, temp_year, "Soil temperature (\u00B0C)",      2),
  build_rows(vwc_marg,  vwc_season,  vwc_year,  "VWC (m\u00B3 m\u207B\u00B3)",      3)
) |>
  rename(`Warming effect` = Estimate, `95% CI` = CI, `p (unadjusted)` = p)
ft2 <- style_ft(flextable(tab2)) |> merge_v(j = "Response")

# ============================================================================
# TABLE 3: Apparent Q10 and VWC optimum by treatment (modification model)
# ============================================================================
q10_tbl <- read_safe(file.path(tab_dir, "q10_by_treatment.csv"))
vwc_opt <- read_safe(file.path(tab_dir, "vwc_optimum_by_treatment.csv"))

q10_wood <- tibble(
  treatment   = c("control", "warmed"),
  Q10_wood    = c(2.51, 0.71),
  Q10_wood_sd = c(1.23, 1.30)
)

tab3 <- q10_tbl |>
  left_join(q10_wood, by = "treatment") |>
  left_join(vwc_opt |> select(treatment, optimum), by = "treatment") |>
  transmute(
    Treatment              = treatment,
    `Q10 (this study)`     = fmt_n(Q10),
    `Q10 95% CI`           = ci(Q10_lo, Q10_hi),
    `Q10 (Wood 2025)`      = fmt_n(Q10_wood),
    `VWC optimum (m³ m⁻³)` = fmt_n(optimum, 3)
  )
ft3 <- style_ft(flextable(tab3))

q10_red_pct <- if (!is.null(q10_tbl) && all(c("control","warmed") %in% q10_tbl$treatment))
  100 * (q10_tbl$Q10[q10_tbl$treatment == "control"] -
         q10_tbl$Q10[q10_tbl$treatment == "warmed"]) /
        q10_tbl$Q10[q10_tbl$treatment == "control"] else NA_real_
vwc_opt_shift <- if (!is.null(vwc_opt) && all(c("control","warmed") %in% vwc_opt$treatment))
  vwc_opt$optimum[vwc_opt$treatment == "control"] -
  vwc_opt$optimum[vwc_opt$treatment == "warmed"] else NA_real_

tab3_note <- sprintf(
  "Apparent Q10 changed by %s%% under warming (control %s, warmed %s). VWC optimum shifted toward drier soil by %s m\u00B3 m\u207B\u00B3 (control %s, warmed %s). Wood et al. (2025) reported a first-experimental-year Q10 of 2.51 (control) and 0.71 (warmed); SD 1.23 and 1.30.",
  fmt_n(q10_red_pct, 1),
  fmt_n(q10_tbl$Q10[q10_tbl$treatment == "control"], 2),
  fmt_n(q10_tbl$Q10[q10_tbl$treatment == "warmed"],  2),
  fmt_n(vwc_opt_shift, 3),
  fmt_n(vwc_opt$optimum[vwc_opt$treatment == "control"], 3),
  fmt_n(vwc_opt$optimum[vwc_opt$treatment == "warmed"],  3)
)

q10_c  <- read_safe(file.path(tab_dir, "q10_contrast.csv"))
vopt_c <- read_safe(file.path(tab_dir, "vwc_optimum_contrast.csv"))
contrast_note <- if (!is.null(q10_c) && !is.null(vopt_c) &&
                     !"status" %in% names(q10_c) && !"status" %in% names(vopt_c)) {
  sprintf(
    "Temperature-sensitivity contrast (warmed minus control slope of log flux on temperature): delta = %s (SE %s, 95%% CI [%s, %s]), p = %s; equivalent Q10 ratio = %s (95%% CI [%s, %s]). Moisture-optimum contrast (warmed minus control quadratic peak): delta = %s m\u00B3 m\u207B\u00B3 (SE %s, 95%% CI [%s, %s]), p = %s.",
    fmt_n(q10_c$delta_slope, 5), fmt_n(q10_c$SE, 5),
    fmt_n(q10_c$delta_lo,    5), fmt_n(q10_c$delta_hi, 5),
    fmt_p(q10_c$p_value),
    fmt_n(q10_c$Q10_ratio,    3), fmt_n(q10_c$Q10_ratio_lo, 3),
    fmt_n(q10_c$Q10_ratio_hi, 3),
    fmt_n(vopt_c$delta_opt,   3), fmt_n(vopt_c$SE,         3),
    fmt_n(vopt_c$delta_lo,    3), fmt_n(vopt_c$delta_hi,   3),
    fmt_p(vopt_c$p_value))
} else "(contrast tests not available)"

# ============================================================================
# TABLE S1: Mechanistic coefficients (pathway and modification)
# ============================================================================
mech_path_coef <- read_safe(file.path(tab_dir, "coef_mech_pathway.csv"))
mech_mod_coef  <- read_safe(file.path(tab_dir, "coef_mech_modification.csv"))
build_coef <- function(df, model_label) {
  if (is.null(df)) return(NULL)
  df |>
    filter(component == "conditional") |>
    transmute(Model    = model_label,
              Term     = term,
              Estimate = fmt_n(Estimate, 3),
              SE       = fmt_n(`Std. Error`, 3),
              p        = fmt_p(`Pr(>|z|)`))
}
tabS1 <- bind_rows(
  build_coef(mech_path_coef, "Pathway"),
  build_coef(mech_mod_coef,  "Modification")
)
if (nrow(tabS1) == 0) tabS1 <- tibble(Term = "(coefficient files missing)")
ftS1 <- style_ft(flextable(tabS1)) |> merge_v(j = "Model")

# ============================================================================
# TABLE S2: Convergence and AR(1) summary (tau in STEPS; eight models)
# ============================================================================
summarize_fit_row <- function(label, res) {
  if (is.null(res) || is.null(res$fit)) {
    return(tibble(Model = label, Status = (res$tactic %||% "FAILED"),
                  pdHess = "-", AIC = "-",
                  `sigma(Plot)` = "-", `sigma(Plot_Year)` = "-",
                  phi = "-", `tau (steps)` = "-"))
  }
  fit <- res$fit
  vc <- tryCatch(VarCorr(fit)$cond, error = function(e) NULL)
  sigma_plot <- NA_real_; sigma_py <- NA_real_; phi <- NA_real_; tau <- NA_real_
  if (!is.null(vc)) {
    for (i in seq_along(vc)) {
      sd_i <- as.numeric(attr(vc[[i]], "stddev"))[1]
      cm <- attr(vc[[i]], "correlation")
      if (!is.null(cm) && is.matrix(cm) && nrow(cm) >= 2) {
        phi      <- as.numeric(cm[1, 2])
        sigma_py <- sd_i
        if (!is.na(phi) && phi > 0 && phi < 1) tau <- -1 / log(phi)
      } else {
        sigma_plot <- sd_i
      }
    }
  }
  pd  <- tryCatch(isTRUE(fit$sdr$pdHess), error = function(e) FALSE)
  aic <- tryCatch(AIC(fit),               error = function(e) NA_real_)
  tibble(Model = label,
         Status            = (res$tactic %||% "standard"),
         pdHess            = if (pd) "TRUE" else "FALSE",
         AIC               = fmt_n(aic,        1),
         `sigma(Plot)`     = fmt_n(sigma_plot, 4),
         `sigma(Plot_Year)`= fmt_n(sigma_py,   4),
         phi               = fmt_n(phi,        4),
         `tau (steps)`     = fmt_n(tau,        2))
}
tabS2 <- bind_rows(
  summarize_fit_row("Temperature (pooled)",      obj$fits$temp_pooled),
  summarize_fit_row("Temperature (interaction)", obj$fits$temp_inter),
  summarize_fit_row("VWC (pooled)",              obj$fits$vwc_pooled),
  summarize_fit_row("VWC (interaction)",         obj$fits$vwc_inter),
  summarize_fit_row("Flux (pooled)",             obj$fits$flux_pooled),
  summarize_fit_row("Flux (interaction)",        obj$fits$flux_inter),
  summarize_fit_row("Mechanistic (pathway)",     obj$fits$mech_path),
  summarize_fit_row("Mechanistic (modification)",obj$fits$mech_mod)
)
ftS2 <- style_ft(flextable(tabS2))

# ============================================================================
# TABLE S3: Per-plot flux summary
# ============================================================================
d_flux       <- obj$d_flux
ctrl_geomean <- d_flux %>% filter(Treatment == "control") %>%
  summarise(g = exp(mean(log_flux, na.rm = TRUE))) %>% pull(g)
per_plot <- d_flux %>%
  group_by(Plot, Treatment) %>%
  summarise(geo_flux = exp(mean(log_flux, na.rm = TRUE)),
            n_days   = dplyr::n(), .groups = "drop") %>%
  mutate(ratio_to_control = if_else(Treatment == "warmed",
                                    geo_flux / ctrl_geomean, NA_real_))
tabS3 <- per_plot %>%
  transmute(Plot                    = as.character(Plot),
            Treatment,
            `Geometric-mean flux`   = fmt_n(geo_flux),
            `Ratio to control mean` = ifelse(is.na(ratio_to_control), "-",
                                             fmt_n(ratio_to_control)),
            `n plot-days`           = n_days)
ftS3 <- style_ft(flextable(tabS3))

# ============================================================================
# TABLE S4: Data coverage
# ============================================================================
cov_tv <- read_safe(file.path(tab_dir, "coverage_tv_year.csv")) %>%
  { if (is.null(.)) NULL else mutate(., Pipeline = "Temp/VWC") }
cov_fl <- read_safe(file.path(tab_dir, "coverage_flux_year.csv")) %>%
  { if (is.null(.)) NULL else mutate(., Pipeline = "Flux") }
cov_combined <- bind_rows(cov_tv, cov_fl)
tabS4 <- if (nrow(cov_combined) > 0) {
  cov_combined %>%
    select(Pipeline, Plot, Treatment, Year, n_days) %>%
    pivot_wider(names_from = Year, values_from = n_days, values_fill = 0L) %>%
    mutate(Plot = as.character(Plot)) %>%
    arrange(Pipeline, Treatment, Plot)
} else tibble(Pipeline = "(no coverage data)")
ftS4 <- style_ft(flextable(tabS4)) |> merge_v(j = "Pipeline")

# ============================================================================
# TABLE S5: Treatment effect on day-to-day variability (empirical)
# ============================================================================
var_csv <- read_safe(file.path(tab_dir, "variability_summary.csv"))
scale_lab <- c(temperature   = "SD, \u00B0C",
               vwc           = "SD, vwc_pct",
               `flux (log)`  = "SD, log units")
tabS5 <- if (!is.null(var_csv) && !"status" %in% names(var_csv)) {
  var_csv |>
    transmute(
      Response            = model,
      Scale               = unname(scale_lab[model]),
      `n cells`           = n_cells,
      `SD ratio (W:C)`    = fmt_n(SD_ratio, 3),
      `95% CI`            = ci(SD_ratio_lo, SD_ratio_hi, 3),
      p                   = fmt_p(p_value)
    )
} else tibble(Response = "(variability summary missing)")
ftS5 <- style_ft(flextable(tabS5))

# ============================================================================
# MAIN TABLES DOCX
# ============================================================================
doc_main <- read_docx() |>
  body_add_par("Table 1. Warming effect on soil CO2 flux.", style = "heading 2") |>
  body_add_par("Geometric-mean warmed:control flux ratio with 95% confidence intervals and percent change. The overall row is the balanced long-term marginal effect from the pooled flux model (log_flux ~ Treatment + Year_f + season); per-season and per-year rows are from the interaction model (log_flux ~ Treatment * Year_f + Treatment * season). All models: (1|Plot) + ar1(time_ou + 0 | Plot_Year), dispformula ~ 1. p-values are unadjusted for pre-specified comparisons.", style = "Normal") |>
  body_add_flextable(ft1) |>
  body_add_par("", style = "Normal") |>
  body_add_par("Table 2. Warming effect on soil temperature and soil moisture.", style = "heading 2") |>
  body_add_par("Overall rows are marginal effects from the pooled models (response ~ Treatment + Year_f + season); per-season and per-year rows are from the interaction models (response ~ Treatment * Year_f + Treatment * season). VWC fitted on the vwc_pct scale (vwc_mean * 100) and reported here divided by 100 (m3 m-3). All models: (1|Plot) + ar1(time_ou + 0 | Plot_Year), dispformula ~ 1. p-values are unadjusted.", style = "Normal") |>
  body_add_flextable(ft2) |>
  body_add_par("", style = "Normal") |>
  body_add_par("Table 3. Apparent temperature sensitivity (Q10) and VWC optimum of soil respiration.", style = "heading 2") |>
  body_add_par(tab3_note, style = "Normal") |>
  body_add_par("Q10 = exp(10 * slope) of log flux on soil temperature per treatment, from the mechanistic modification model (log_flux ~ Treatment + Year_f + season + temp_c + vwc_c + I(vwc_c^2) + Treatment:(temp_c + vwc_c + I(vwc_c^2)) + (1|Plot) + ar1(time_ou + 0 | Plot_Year), dispformula ~ 1; temp_c and vwc_c grand-mean centred). Control slope is temp_c; warmed slope is temp_c + Treatmentwarmed:temp_c. VWC optimum is the per-treatment quadratic peak -b1/(2 b2) on the centred scale plus the centring constant. 95% CIs by delta method.", style = "Normal") |>
  body_add_flextable(ft3) |>
  body_add_par("Contrast tests (pre-specified):", style = "Normal") |>
  body_add_par(contrast_note, style = "Normal")

print(doc_main, target = file.path(out_dir, "Manuscript_Tables_Main.docx"))

# ============================================================================
# SUPPLEMENTARY TABLES DOCX
# ============================================================================
doc_supp <- read_docx() |>
  body_add_par("Table S1. Mechanistic model coefficients (conditional component).",
               style = "heading 2") |>
  body_add_par("Fixed-effect coefficients on the log-flux scale for the mechanistic pathway model (common environmental slopes) and modification model (treatment-specific slopes). Covariates temp_c and vwc_c are grand-mean centred. The modification model's Treatment:temp_c gives the temperature-sensitivity difference (Table 3); the Treatment:vwc_c and Treatment:I(vwc_c^2) terms give the moisture-optimum difference.", style = "Normal") |>
  body_add_flextable(ftS1) |>
  body_add_par("", style = "Normal") |>
  body_add_par("Table S2. Model convergence and AR(1) random-structure summary.",
               style = "heading 2") |>
  body_add_par("Random structure for every model: (1|Plot) + ar1(time_ou + 0 | Plot_Year), dispformula ~ 1. phi is the AR(1) lag-1 autocorrelation; tau = -1/log(phi) is the decorrelation time in observation steps (consecutive retained plot-days), not calendar days, because ar1() treats consecutive numFactor levels as equally spaced. Status reports the fitting tactic; a dispersed-start retry is used where the standard fit does not yield a positive-definite Hessian.", style = "Normal") |>
  body_add_flextable(ftS2) |>
  body_add_par("", style = "Normal") |>
  body_add_par("Table S3. Per-plot flux summary.", style = "heading 2") |>
  body_add_par("Geometric-mean daily flux per plot from the flux pipeline (Sep 2018 onward); warmed plots also expressed as a ratio to the pooled control geometric mean, to show plot-to-plot (slope-position) heterogeneity.", style = "Normal") |>
  body_add_flextable(ftS3) |>
  body_add_par("", style = "Normal") |>
  body_add_par("Table S4. Data coverage (plot-days per plot per year).",
               style = "heading 2") |>
  body_add_par("Pipeline-wise coverage; temp/VWC and flux datasets are filtered identically by date (>= 2018-09-01) and daily QC (>= 12 hourly records per plot-day) but differ in that the flux pipeline additionally requires Flux > 0.", style = "Normal") |>
  body_add_flextable(ftS4) |>
  body_add_par("", style = "Normal") |>
  body_add_par("Table S5. Treatment effect on day-to-day variability (exploratory).",
               style = "heading 2") |>
  body_add_par("Warmed:control ratio of the within-plot-year-season standard deviation of the daily series, from a mixed model log(SD) ~ Treatment + season + (1|Plot) on the per-cell dispersions with a finite-sample bias correction. SD ratio = exp(Treatment coefficient). Dispersion is in each response's modelling scale (degrees Celsius, vwc_pct, natural-log units for flux; the log-scale SD approximates the coefficient of variation of flux). Exploratory: cell-level SD uncertainty and residual temporal dependence are not propagated.", style = "Normal") |>
  body_add_flextable(ftS5)

print(doc_supp, target = file.path(out_dir, "Manuscript_Tables_Supplementary.docx"))

cat("\nTables written:\n")
cat(" -", file.path(out_dir, "Manuscript_Tables_Main.docx"),          "\n")
cat(" -", file.path(out_dir, "Manuscript_Tables_Supplementary.docx"), "\n")
