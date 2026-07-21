# ============================================================================
# GLMM_VIS.R
# ----------------------------------------------------------------------------
# Publication figures from the GLMM outputs. Reads the RDS of fitted objects
# and the contrast CSVs; fits no new inferential models.
#
# Model sourcing:
#   Treatment x season EMMs        interaction models
#   Year-resolved warming effects  interaction-model contrast CSVs
#   Q10 and response surface        mechanistic modification model (centred
#                                    covariates temp_c, vwc_c; axes shown on the
#                                    real temperature and moisture scale)
#   Variability                    variability_summary.csv (empirical per-cell SD)
#
# Time-series figures: segments connecting measurement points more than one day
# apart are drawn dashed and in a lighter tint; area fills use a light tint.
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(tibble)
  library(stringr); library(lubridate); library(purrr)
  library(ggplot2); library(patchwork); library(viridisLite)
  library(emmeans); library(glmmTMB)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

setwd("D:/USDA/TRACE_DM_APR_26/CO2")
out_dir <- "OUTPUT_GLMM_FINAL"
fig_dir <- file.path(out_dir, "figures")
tab_dir <- file.path(out_dir, "tables")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

trt_cols <- c(control = "#0571b0", warmed = "#ca0020")

# Lighten a colour toward white by a fraction (for gap segments and fills).
lighten <- function(col, amount = 0.55) {
  rgb_v <- col2rgb(col) / 255
  out   <- rgb_v + (1 - rgb_v) * amount
  rgb(out[1], out[2], out[3])
}

theme_pub <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title       = element_text(face = "bold", size = base_size + 1),
      plot.subtitle    = element_text(size = base_size - 1, colour = "grey30"),
      axis.title       = element_text(face = "bold"),
      axis.text        = element_text(colour = "black"),
      strip.background = element_rect(fill = "grey90"),
      strip.text       = element_text(face = "bold"),
      legend.position  = "bottom"
    )
}

save_fig <- function(plot, name, width, height, dpi = 600) {
  ggsave(file.path(fig_dir, paste0(name, ".png")), plot,
         width = width, height = height, dpi = dpi)
  ggsave(file.path(fig_dir, paste0(name, ".pdf")), plot,
         width = width, height = height)
}

read_safe <- function(path) {
  if (!file.exists(path)) { warning("missing: ", path); return(NULL) }
  read_csv(path, show_col_types = FALSE)
}

# ---------------------------------------------------------------------------
# Gap-aware line segmentation. Builds an explicit segment table between
# consecutive points; segments spanning more than gap_days are flagged is_gap.
# A gap-styled line is then two layers: solid segments (is_gap FALSE) and
# dashed lightened segments (is_gap TRUE).
# ---------------------------------------------------------------------------
segment_series <- function(df, date_col = "Date", y_col, group_col = NULL,
                           gap_days = 1) {
  build_one <- function(d) {
    d <- d[order(d[[date_col]]), ]
    d <- d[is.finite(d[[y_col]]), ]
    n <- nrow(d)
    if (n < 2) return(NULL)
    gap <- as.numeric(d[[date_col]][-1] - d[[date_col]][-n])
    tibble(x    = d[[date_col]][-n], xend = d[[date_col]][-1],
           y    = d[[y_col]][-n],    yend = d[[y_col]][-1],
           is_gap = gap > gap_days)
  }
  if (is.null(group_col)) return(build_one(df))
  df %>% group_split(.data[[group_col]]) %>%
    map_dfr(function(g) {
      s <- build_one(g)
      if (is.null(s)) return(NULL)
      s[[group_col]] <- g[[group_col]][1]
      s
    })
}

# Add a gap-styled line (solid + dashed-light) for one coloured series.
add_gap_line <- function(seg, colour, linewidth = 0.6) {
  list(
    geom_segment(data = dplyr::filter(seg, !is_gap),
                 aes(x = x, y = y, xend = xend, yend = yend),
                 colour = colour, linewidth = linewidth),
    geom_segment(data = dplyr::filter(seg, is_gap),
                 aes(x = x, y = y, xend = xend, yend = yend),
                 colour = lighten(colour, 0.55), linewidth = linewidth,
                 linetype = "dashed")
  )
}

# ---------------------------- load outputs ----------------------------------
obj <- readRDS(file.path(out_dir, "objects", "glmm_v3_objects.rds"))

d_analysis <- obj$d_analysis
d_flux     <- obj$d_flux
m_temp_i   <- obj$fits$temp_inter$fit
m_vwc_i    <- obj$fits$vwc_inter$fit
m_flux_i   <- obj$fits$flux_inter$fit
m_mech     <- obj$fits$mech_mod$fit
TEMP_CENTRE <- obj$config$TEMP_CENTRE
VWC_CENTRE  <- obj$config$VWC_CENTRE

# ============================================================================
# MAIN FIGURE 1: Treatment x season estimated marginal means
#   From the interaction models (all three responses now carry Treatment:season).
#   VWC EMMs rescaled vwc_pct -> m^3 m^-3.
# ============================================================================
emm_flux_s <- as.data.frame(emmeans(m_flux_i, ~ Treatment | season, data = d_flux)) |>
  mutate(est = exp(emmean), lo = exp(lower.CL), hi = exp(upper.CL))
emm_temp_s <- as.data.frame(emmeans(m_temp_i, ~ Treatment | season, data = d_analysis)) |>
  mutate(est = emmean, lo = lower.CL, hi = upper.CL)
emm_vwc_s  <- as.data.frame(emmeans(m_vwc_i, ~ Treatment | season, data = d_analysis)) |>
  mutate(est = emmean / 100, lo = lower.CL / 100, hi = upper.CL / 100)

panel_season <- function(df, ylab, title) {
  ggplot(df, aes(x = season, y = est, colour = Treatment, group = Treatment)) +
    geom_point(position = position_dodge(0.3), size = 3) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.18,
                  position = position_dodge(0.3), linewidth = 0.8) +
    geom_line(position = position_dodge(0.3), linewidth = 0.9) +
    scale_colour_manual(values = trt_cols) +
    labs(x = NULL, y = ylab, title = title) +
    theme_pub()
}

fig1 <- (panel_season(emm_flux_s, expression(CO[2]~flux~(mu*mol~m^-2~s^-1)), "(A) CO2 flux") /
         panel_season(emm_temp_s, "Soil temperature (\u00B0C)",              "(B) Soil temperature") /
         panel_season(emm_vwc_s,  expression(VWC~(m^3~m^-3)),                 "(C) Soil moisture")) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")
save_fig(fig1, "Figure1_TreatmentSeason", width = 6, height = 12)

# ============================================================================
# MAIN FIGURE 2: Year-resolved warming effect (from interaction-model contrasts)
# ============================================================================
valid_contrast <- function(x, required = c("Year_f", "estimate", "lower.CL", "upper.CL")) {
  !is.null(x) && !"status" %in% names(x) && all(required %in% names(x))
}
cy_flux <- read_safe(file.path(tab_dir, "contrast_flux_by_year.csv"))
cy_temp <- read_safe(file.path(tab_dir, "contrast_temp_by_year.csv"))
cy_vwc  <- read_safe(file.path(tab_dir, "contrast_vwc_by_year.csv"))
if (valid_contrast(cy_flux)) cy_flux <- cy_flux %>%
  mutate(Year_f = factor(Year_f), ratio_warmed_control = exp(estimate),
         ratio_lower = exp(lower.CL), ratio_upper = exp(upper.CL))
if (valid_contrast(cy_temp)) cy_temp <- cy_temp %>% mutate(Year_f = factor(Year_f))
if (valid_contrast(cy_vwc)) cy_vwc <- cy_vwc %>%
  mutate(Year_f = factor(Year_f), estimate = estimate / 100,
         lower.CL = lower.CL / 100, upper.CL = upper.CL / 100)

forest_panel <- function(df, xval, xlo, xhi, ref, xlab, title, col) {
  ggplot(df, aes(x = .data[[xval]], y = Year_f)) +
    geom_vline(xintercept = ref, linetype = 2, colour = "grey40") +
    geom_errorbarh(aes(xmin = .data[[xlo]], xmax = .data[[xhi]]),
                   height = 0.25, linewidth = 0.7, colour = col) +
    geom_point(size = 2.6, colour = col) +
    scale_y_discrete(limits = rev) +
    labs(x = xlab, y = NULL, title = title) +
    theme_pub() + theme(legend.position = "none")
}
if (valid_contrast(cy_flux) && valid_contrast(cy_temp) && valid_contrast(cy_vwc)) {
p2a <- forest_panel(cy_flux, "ratio_warmed_control", "ratio_lower", "ratio_upper",
                    ref = 1, xlab = "Warmed : control flux ratio",
                    title = "(A) CO2 flux", col = "#ca0020")
p2b <- forest_panel(cy_temp, "estimate", "lower.CL", "upper.CL",
                    ref = 0, xlab = "Warming effect on soil temperature (\u00B0C)",
                    title = "(B) Soil temperature", col = "#b2182b")
p2c <- forest_panel(cy_vwc, "estimate", "lower.CL", "upper.CL",
                    ref = 0, xlab = expression(Warming~effect~on~VWC~(m^3~m^-3)),
                    title = "(C) Soil moisture", col = "#2166ac")
fig2 <- (p2a | p2b | p2c) +
  plot_annotation(title = "Year-resolved warming effects",
                  theme = theme(plot.title = element_text(face = "bold", size = 13)))
save_fig(fig2, "Figure2_YearForest", width = 13, height = 5)
} else {
  warning("Year-resolved contrast output unavailable; Figure 2 was not created.")
}

# ============================================================================
# MAIN FIGURE 2b: Season-resolved warming effect (from interaction contrasts)
#   New: all three responses carry Treatment:season, so per-season warmed-
#   control effects are shown together.
# ============================================================================
cs_flux <- read_safe(file.path(tab_dir, "contrast_flux_by_season.csv"))
cs_temp <- read_safe(file.path(tab_dir, "contrast_temp_by_season.csv"))
cs_vwc  <- read_safe(file.path(tab_dir, "contrast_vwc_by_season.csv"))
if (!is.null(cs_flux) && !is.null(cs_temp) && !is.null(cs_vwc) &&
    !"status" %in% names(cs_flux) && !"status" %in% names(cs_temp) &&
    !"status" %in% names(cs_vwc)) {
  cs_flux <- cs_flux %>% mutate(ratio = exp(estimate),
                                lo = exp(lower.CL), hi = exp(upper.CL))
  cs_vwc  <- cs_vwc  %>% mutate(estimate = estimate / 100,
                                lower.CL = lower.CL / 100, upper.CL = upper.CL / 100)
  s2a <- ggplot(cs_flux, aes(x = ratio, y = season)) +
    geom_vline(xintercept = 1, linetype = 2, colour = "grey40") +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.2, colour = "#ca0020") +
    geom_point(size = 3, colour = "#ca0020") +
    labs(x = "Warmed : control flux ratio", y = NULL, title = "(A) CO2 flux") +
    theme_pub() + theme(legend.position = "none")
  s2b <- ggplot(cs_temp, aes(x = estimate, y = season)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
    geom_errorbarh(aes(xmin = lower.CL, xmax = upper.CL), height = 0.2, colour = "#b2182b") +
    geom_point(size = 3, colour = "#b2182b") +
    labs(x = "Warming effect (\u00B0C)", y = NULL, title = "(B) Soil temperature") +
    theme_pub() + theme(legend.position = "none")
  s2c <- ggplot(cs_vwc, aes(x = estimate, y = season)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
    geom_errorbarh(aes(xmin = lower.CL, xmax = upper.CL), height = 0.2, colour = "#2166ac") +
    geom_point(size = 3, colour = "#2166ac") +
    labs(x = expression(Warming~effect~on~VWC~(m^3~m^-3)), y = NULL,
         title = "(C) Soil moisture") +
    theme_pub() + theme(legend.position = "none")
  fig2b <- (s2a | s2b | s2c) +
    plot_annotation(title = "Season-resolved warming effects",
                    theme = theme(plot.title = element_text(face = "bold", size = 13)))
  save_fig(fig2b, "Figure2b_SeasonForest", width = 13, height = 4)
}

# ============================================================================
# MAIN FIGURE 3: Apparent Q10 by treatment
# ============================================================================
q10_trt <- read_safe(file.path(tab_dir, "q10_by_treatment.csv"))
if (!is.null(q10_trt) && !"status" %in% names(q10_trt) &&
    all(c("treatment", "Q10", "Q10_lo", "Q10_hi") %in% names(q10_trt))) {
  q10_trt <- q10_trt %>% rename(Treatment = treatment)
  p3 <- ggplot(q10_trt, aes(x = Treatment, y = Q10, colour = Treatment)) +
    geom_hline(yintercept = 1, linetype = 3, colour = "grey50") +
    geom_errorbar(aes(ymin = Q10_lo, ymax = Q10_hi), width = 0.12, linewidth = 0.9) +
    geom_point(size = 4) +
    scale_colour_manual(values = trt_cols, guide = "none") +
    labs(x = NULL, y = expression(Apparent~Q[10]),
         title    = "Apparent temperature sensitivity of soil respiration",
         subtitle = "Multi-year estimates from the mechanistic modification model.") +
    theme_pub()
  save_fig(p3, "Figure3_Q10_by_treatment", width = 6, height = 6)
}

# ============================================================================
# MAIN FIGURE 4: Percent-change response surface (modification model)
#   Built on centred covariates temp_c, vwc_c; axes shown on the real scale
#   (temp_c + TEMP_CENTRE, vwc_c + VWC_CENTRE).
# ============================================================================
vwc_seq_c  <- seq(min(d_flux$vwc_c,  na.rm = TRUE),
                  max(d_flux$vwc_c,  na.rm = TRUE), length.out = 40)
temp_seq_c <- seq(min(d_flux$temp_c, na.rm = TRUE),
                  max(d_flux$temp_c, na.rm = TRUE), length.out = 40)

emm_surf <- emmeans(m_mech, ~ Treatment * vwc_c * temp_c,
                    at       = list(vwc_c = vwc_seq_c, temp_c = temp_seq_c),
                    data     = d_flux, rg.limit = 200000)
surf <- as.data.frame(emm_surf) %>%
  mutate(flux = exp(emmean),
         vwc_mean  = vwc_c  + VWC_CENTRE,
         temp_mean = temp_c + TEMP_CENTRE)
surf_wide <- surf %>%
  select(vwc_mean, temp_mean, Treatment, flux) %>%
  pivot_wider(names_from = Treatment, values_from = flux) %>%
  mutate(diff = warmed - control,
         pct  = (warmed - control) / control * 100)
write_csv(surf_wide, file.path(tab_dir, "mech_surface_difference.csv"))

p4 <- ggplot(surf_wide, aes(x = vwc_mean, y = temp_mean, fill = pct)) +
  geom_raster(interpolate = TRUE) +
  geom_contour(aes(z = pct), colour = "white", alpha = 0.5, linewidth = 0.3) +
  scale_fill_gradient2(low = "#0571b0", mid = "white", high = "#ca0020",
                       midpoint = 0, name = "% change") +
  labs(x = expression(Soil~VWC~(m^3~m^-3)), y = "Soil temperature (\u00B0C)",
       title    = "Warming effect on CO2 flux across the moisture-temperature domain",
       subtitle = "Percent change in model-estimated geometric-mean flux (warmed vs control), from the mechanistic modification model") +
  theme_pub() + theme(legend.position = "right")
save_fig(p4, "Figure4_PercentChange_Surface", width = 8, height = 6)

# ============================================================================
# MAIN FIGURE 5: Daily delta time series with gap-aware line styling
#   Gap segments (> 1 day between measured plot-days) are dashed and lightened;
#   the delta ribbon uses a light tint.
# ============================================================================
if (file.exists("HURCN.csv")) {
  hurcn <- read.csv("HURCN.csv", stringsAsFactors = FALSE)
  if ("Start.Date" %in% names(hurcn))
    hurcn <- hurcn |> rename(Start = Start.Date, End = Finish.Date)
  hurcn <- hurcn |>
    mutate(Start = as.Date(parse_date_time(Start, orders = c("dmy","mdy","ymd"))),
           End   = as.Date(parse_date_time(End,   orders = c("dmy","mdy","ymd")))) |>
    filter(!is.na(Start), !is.na(End))
  data_start <- min(c(d_flux$Date, d_analysis$Date), na.rm = TRUE)
  data_end   <- max(c(d_flux$Date, d_analysis$Date), na.rm = TRUE)
  # Retain only events overlapping the observed data window and clip ranges so
  # annotations cannot extend the time-series axis beyond the measurements.
  hurcn <- hurcn |> filter(End >= data_start, Start <= data_end) |>
    mutate(Start = pmax(Start, data_start), End = pmin(End, data_end))
  disturbances_ranges <- hurcn |> filter(Start != End)
  disturbances_points <- hurcn |> filter(Start == End)
} else {
  warning("HURCN.csv not found; delta plots drawn without disturbance overlays.")
  hurcn               <- data.frame()
  disturbances_ranges <- data.frame(Start = as.Date(character()), End = as.Date(character()))
  disturbances_points <- data.frame(Start = as.Date(character()))
}

data_years    <- sort(unique(c(d_flux$Year, d_analysis$Year)))
custom_breaks <- sort(as.Date(paste0(data_years, "-01-01")))

create_delta_plot <- function(daily_data, var_name, y_lab, title_txt,
                              line_col, show_x = FALSE) {
  # Daily treatment differences are shown descriptively. The ribbon is the
  # same-day between-plot SE, not a model-based confidence interval.
  agg <- daily_data |>
    group_by(Date, Treatment) |>
    summarise(mean_val = mean(.data[[var_name]], na.rm = TRUE),
              se_val   = sd(.data[[var_name]], na.rm = TRUE) /
                         sqrt(sum(!is.na(.data[[var_name]]))),
              .groups  = "drop") |>
    pivot_wider(names_from = Treatment, values_from = c(mean_val, se_val)) |>
    mutate(Delta     = mean_val_warmed - mean_val_control,
           Delta_SE  = sqrt(se_val_warmed^2 + se_val_control^2),
           Delta_Min = Delta - Delta_SE,
           Delta_Max = Delta + Delta_SE) |>
    filter(is.finite(Delta), is.finite(Delta_SE)) |>
    arrange(Date)

  seg <- segment_series(agg, date_col = "Date", y_col = "Delta", gap_days = 1)

  p <- ggplot() +
    geom_rect(data = disturbances_ranges,
              aes(xmin = Start, xmax = End, ymin = -Inf, ymax = Inf,
                  fill = "Storm/Hurricane"), alpha = 0.4, inherit.aes = FALSE) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "black", linewidth = 0.5) +
    geom_ribbon(data = agg, aes(x = Date, ymin = Delta_Min, ymax = Delta_Max),
                fill = lighten(line_col, 0.6), alpha = 0.5) +
    add_gap_line(seg, line_col, linewidth = 0.5) +
    {if (nrow(disturbances_points) > 0)
      geom_vline(data = disturbances_points, aes(xintercept = Start),
                 colour = "grey30", linewidth = 0.4, linetype = "solid")} +
    scale_x_date(breaks = custom_breaks, date_labels = "%Y") +
    scale_fill_manual(name = NULL, values = c("Storm/Hurricane" = "wheat")) +
    labs(y = y_lab, x = NULL, title = title_txt) +
    theme_pub() +
    theme(legend.position = if (show_x) "bottom" else "none")
  if (!show_x) p <- p + theme(axis.text.x = element_blank())
  p
}

f5a <- create_delta_plot(d_flux,     "flux_mean",
        expression(Delta~flux~(mu*mol~m^-2~s^-1)),      "(A) CO2 flux",          "#006400")
f5b <- create_delta_plot(d_analysis, "temp_mean",
        expression(Delta~temperature~"("*degree*C*")"), "(B) Soil temperature",  "#330066")
f5c <- create_delta_plot(d_analysis, "vwc_mean",
        expression(Delta~VWC~(m^3~m^-3)),               "(C) Soil moisture",     "#E69F00",
        show_x = TRUE)
fig5 <- f5a / f5b / f5c +
  plot_annotation(subtitle = "Descriptive daily warmed-minus-control differences; ribbons are same-day between-plot SEs, not model-based confidence intervals.")
save_fig(fig5, "Figure5_DailyDelta", width = 12, height = 10)

# ============================================================================
# MAIN FIGURE 6: Treatment effect on day-to-day variability
#   Point estimate and CI of the warmed:control SD ratio from the pipeline
#   variability_summary.csv (empirical per-cell log-SD analysis).
# ============================================================================
var_sum <- read_safe(file.path(tab_dir, "variability_summary.csv"))
if (!is.null(var_sum) && !"status" %in% names(var_sum)) {
  lab_map <- c(temperature = "Soil temperature", vwc = "Soil moisture",
               `flux (log)` = "CO2 flux")
  var_sum <- var_sum %>%
    mutate(response = recode(model, !!!lab_map))
  fig6 <- ggplot(var_sum, aes(x = SD_ratio, y = response)) +
    geom_vline(xintercept = 1, linetype = 2, colour = "grey40") +
    geom_errorbarh(aes(xmin = SD_ratio_lo, xmax = SD_ratio_hi),
                   height = 0.18, linewidth = 0.8, colour = "grey25") +
    geom_point(size = 3.2, colour = "#ca0020") +
    labs(x = "Warmed : control SD ratio", y = NULL,
         title    = "Treatment effect on day-to-day variability",
         subtitle = "Warmed:control ratio of within-plot-year-season daily SD (exploratory; ratio > 1 = warmed more variable)") +
    theme_pub()
  save_fig(fig6, "Figure6_Variability", width = 8, height = 4)
}

# ============================================================================
# SUPPLEMENTARY S1: Mechanistic response curves (centred covariates; real axes)
# ============================================================================
make_response_curve <- function(model, focal_c, focal_centre, focal_lab, title) {
  other_c <- setdiff(c("temp_c", "vwc_c"), focal_c)
  focal_seq_c <- seq(min(d_flux[[focal_c]], na.rm = TRUE),
                     max(d_flux[[focal_c]], na.rm = TRUE), length.out = 60)
  at_list <- setNames(list(focal_seq_c, 0), c(focal_c, other_c))
  emm <- emmeans(model, as.formula(paste("~ Treatment *", focal_c)),
                 at = at_list, data = d_flux)
  df <- as.data.frame(emm) |>
    mutate(flux = exp(emmean), lo = exp(lower.CL), hi = exp(upper.CL),
           focal_real = .data[[focal_c]] + focal_centre)
  ggplot(df, aes(x = focal_real, y = flux, colour = Treatment, fill = Treatment)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.22, colour = NA) +
    geom_line(linewidth = 1.1) +
    scale_colour_manual(values = trt_cols) + scale_fill_manual(values = trt_cols) +
    labs(x = focal_lab, y = expression(CO[2]~flux~(mu*mol~m^-2~s^-1)), title = title,
         subtitle = "Model-estimated geometric-mean flux; other driver fixed at its mean") +
    theme_pub()
}
s1 <- (make_response_curve(m_mech, "vwc_c",  VWC_CENTRE, expression(Soil~VWC~(m^3~m^-3)),
                           "(A) Response to soil moisture") |
       make_response_curve(m_mech, "temp_c", TEMP_CENTRE, "Soil temperature (\u00B0C)",
                           "(B) Response to temperature")) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")
save_fig(s1, "FigureS1_ResponseCurves", width = 12, height = 6)

# ============================================================================
# SUPPLEMENTARY S2: Control/warmed surfaces + absolute difference
# ============================================================================
p_s2a <- ggplot(surf, aes(x = vwc_mean, y = temp_mean, z = flux)) +
  geom_contour_filled(bins = 12) +
  facet_wrap(~ Treatment, labeller = as_labeller(c(control = "Control", warmed = "Warmed"))) +
  scale_fill_viridis_d(option = "plasma", name = expression(CO[2]~flux)) +
  labs(x = expression(Soil~VWC~(m^3~m^-3)), y = "Soil temperature (\u00B0C)",
       title    = "(A) Flux response surfaces",
       subtitle = "From the mechanistic modification model") +
  theme_pub() + theme(legend.position = "right")
p_s2b <- ggplot(surf_wide, aes(x = vwc_mean, y = temp_mean, fill = diff)) +
  geom_raster(interpolate = TRUE) +
  geom_contour(aes(z = diff), colour = "white", alpha = 0.5, linewidth = 0.3) +
  scale_fill_gradient2(low = "#0571b0", mid = "white", high = "#ca0020", midpoint = 0,
                       name = expression(Delta~flux)) +
  labs(x = expression(Soil~VWC~(m^3~m^-3)), y = "Soil temperature (\u00B0C)",
       title    = "(B) Absolute warming effect",
       subtitle = "Warmed minus control") +
  theme_pub() + theme(legend.position = "right")
save_fig(p_s2a / p_s2b, "FigureS2_Surfaces", width = 10, height = 10)

# ============================================================================
# SUPPLEMENTARY S3: Hourly descriptive correlations (NOT the inferential model)
# ============================================================================
flux_raw <- read_csv("FLUX.csv", show_col_types = FALSE) |>
  mutate(DayHour     = ymd_hms(DayHour, quiet = TRUE, tz = "UTC"),
         Treatment   = factor(Treatment, levels = c("control", "warmed")),
         Flux        = as.numeric(Flux),
         Temperature = as.numeric(Temperature),
         VWC         = as.numeric(VWC))
d_hourly <- flux_raw |>
  filter(Flux > 0, is.finite(Flux), is.finite(VWC), is.finite(Temperature)) |>
  mutate(log_flux = log(Flux))

descriptive_panel <- function(data, x_col, x_lab, smooth_rhs, title) {
  ggplot(data, aes(x = .data[[x_col]], y = log_flux, colour = Treatment)) +
    geom_point(alpha = 0.04, size = 0.7) +
    geom_smooth(method = "lm", formula = smooth_rhs, colour = "black",
                linewidth = 1.0, se = TRUE) +
    facet_wrap(~ Treatment) +
    scale_colour_manual(values = trt_cols, guide = "none") +
    labs(x = x_lab, y = expression(log(CO[2]~flux)), title = title,
         subtitle = "Descriptive hourly relationship; not the inferential model") +
    theme_pub()
}
s3a <- descriptive_panel(d_hourly, "VWC",         expression(Soil~VWC~(m^3~m^-3)),
                         y ~ poly(x, 2, raw = TRUE), "(A) log-flux vs soil moisture")
s3b <- descriptive_panel(d_hourly, "Temperature", "Soil temperature (\u00B0C)",
                         y ~ x,                       "(B) log-flux vs soil temperature")
save_fig(s3a / s3b, "FigureS3_HourlyCorrelations", width = 12, height = 10)

# ============================================================================
# SUPPLEMENTARY S4: Data coverage / gap structure over time (both pipelines)
# ============================================================================
cov_dat <- bind_rows(
  d_analysis |> distinct(Plot, Date, Treatment) |> mutate(Pipeline = "Temp/VWC"),
  d_flux     |> distinct(Plot, Date, Treatment) |> mutate(Pipeline = "Flux")
) |> mutate(Plot = factor(Plot))
p_s4 <- ggplot(cov_dat, aes(x = Date, y = Plot, colour = Treatment)) +
  geom_point(shape = 15, size = 0.5) +
  facet_wrap(~ Pipeline, ncol = 1) +
  scale_colour_manual(values = trt_cols) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(x = NULL, y = "Plot",
       title    = "Daily measurement coverage by plot",
       subtitle = "Each mark is a retained plot-day (>= 12 valid hours); blank spans are gaps") +
  theme_pub()
save_fig(p_s4, "FigureS4_DataCoverage", width = 12, height = 6)

# ============================================================================
# SUPPLEMENTARY S5: Mediation path schematic (descriptive)
# ============================================================================
temp_marg  <- read_safe(file.path(tab_dir, "contrast_temp_marginal.csv"))
vwc_marg   <- read_safe(file.path(tab_dir, "contrast_vwc_marginal.csv"))
flux_marg  <- read_safe(file.path(tab_dir, "contrast_flux_marginal.csv"))
mech_marg  <- read_safe(file.path(tab_dir, "mech_residual_marginal.csv"))
q10_tbl    <- read_safe(file.path(tab_dir, "q10_by_treatment.csv"))
vwc_opt    <- read_safe(file.path(tab_dir, "vwc_optimum_by_treatment.csv"))

fmt <- function(x, d = 2) ifelse(is.na(x), "n/a", formatC(x, digits = d, format = "f"))

lab_w_temp   <- if (!is.null(temp_marg))
  sprintf("+%s \u00B0C", fmt(temp_marg$estimate[1])) else "n/a"
lab_w_vwc    <- if (!is.null(vwc_marg))
  sprintf("%s m\u00B3 m\u207B\u00B3", fmt(vwc_marg$estimate[1] / 100, 3)) else "n/a"
lab_temp_flux <- if (!is.null(q10_tbl))
  sprintf("Q10 %s / %s",
          fmt(q10_tbl$Q10[q10_tbl$treatment == "control"]),
          fmt(q10_tbl$Q10[q10_tbl$treatment == "warmed"])) else "n/a"
lab_vwc_flux  <- if (!is.null(vwc_opt))
  sprintf("opt %s -> %s",
          fmt(vwc_opt$optimum[vwc_opt$treatment == "control"], 3),
          fmt(vwc_opt$optimum[vwc_opt$treatment == "warmed"],  3)) else "n/a"
lab_direct    <- if (!is.null(mech_marg))
  sprintf("direct: x%s", fmt(exp(mech_marg$estimate[1]))) else "n/a"
lab_total     <- if (!is.null(flux_marg))
  sprintf("total: x%s", fmt(exp(flux_marg$estimate[1]))) else "n/a"

nodes <- tibble(
  x = c(0, 1, 1, 2), y = c(0.5, 0.82, 0.18, 0.5),
  label = c("Warming", "Soil\ntemperature", "Soil\nmoisture", "CO2 flux"))
edges <- tibble(
  x      = c(0.20, 0.20, 1.12, 1.12, 0.20),
  y      = c(0.56, 0.44, 0.80, 0.20, 0.50),
  xend   = c(0.86, 0.86, 1.80, 1.80, 1.80),
  yend   = c(0.80, 0.20, 0.56, 0.44, 0.50),
  lab    = c(lab_w_temp, lab_w_vwc, lab_temp_flux, lab_vwc_flux, lab_direct),
  lx     = c(0.5, 0.5, 1.5, 1.5, 1.0),
  ly     = c(0.74, 0.26, 0.74, 0.26, 0.53),
  dashed = c(FALSE, FALSE, FALSE, FALSE, TRUE))

p_s5 <- ggplot() +
  geom_segment(data = edges,
               aes(x = x, y = y, xend = xend, yend = yend, linetype = dashed),
               arrow = arrow(length = grid::unit(0.18, "cm"), type = "closed"),
               linewidth = 0.6, colour = "grey25") +
  geom_label(data = edges, aes(x = lx, y = ly, label = lab),
             size = 3.2, label.size = 0, fill = "white") +
  geom_label(data = nodes, aes(x = x, y = y, label = label),
             size = 4, fontface = "bold", label.padding = grid::unit(0.4, "lines")) +
  scale_linetype_manual(values = c("FALSE" = "solid", "TRUE" = "22"), guide = "none") +
  scale_x_continuous(limits = c(-0.25, 2.25)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(title    = "Descriptive pathways of the warming effect on CO2 flux",
       subtitle = paste0(lab_total,
                         "; conditional treatment association after T and VWC (not a causal direct effect). A causal mediation analysis would require randomized treatment, pre-treatment mediator--outcome confounder control, longitudinal mediator/outcome models, and sensitivity analysis for mediator--outcome confounding.")) +
  theme_void() +
  theme(plot.title    = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9, colour = "grey30"),
        plot.margin   = margin(12, 12, 12, 12))
save_fig(p_s5, "FigureS5_Mediation_Schematic", width = 9, height = 5)

# ============================================================================
# SUPPLEMENTARY S6: Per-plot Q10 (descriptive; computed on the fly)
# ============================================================================
# A descriptive partial-pooling model with plot-specific temperature slopes and
# the same within-Plot-Year AR(1) structure used by the primary models.
m_plot_q10 <- tryCatch(
  glmmTMB(log_flux ~ Treatment + Year_f + season + temp_c + vwc_c + I(vwc_c^2) +
            Treatment:temp_c + Treatment:vwc_c + Treatment:I(vwc_c^2) +
            (1 + temp_c | Plot) + ar1(time_ou + 0 | Plot_Year),
          dispformula = ~ 1, data = d_flux, REML = TRUE),
  error = function(e) e)
if (inherits(m_plot_q10, "glmmTMB") && isTRUE(m_plot_q10$sdr$pdHess)) {
  fixed <- fixef(m_plot_q10)$cond
  trt_temp_term <- grep("(^Treatmentwarmed:temp_c$|^temp_c:Treatmentwarmed$)",
                        names(fixed), value = TRUE)
  if (length(trt_temp_term) != 1L) {
    stop("Could not identify the warmed-by-temperature fixed-effect term.")
  }
  plot_re <- ranef(m_plot_q10)$cond$Plot
  slopes <- setNames(plot_re[, "temp_c"], rownames(plot_re))
  plot_trt <- d_flux |> distinct(Plot, Treatment)
  q10_plot <- plot_trt |>
    mutate(fixed_slope = if_else(Treatment == "warmed",
                                 fixed[["temp_c"]] + fixed[[trt_temp_term]],
                                 fixed[["temp_c"]]),
           slope = fixed_slope + slopes[as.character(Plot)],
           Q10 = exp(10 * slope)) |>
    filter(is.finite(Q10))
  write_csv(q10_plot, file.path(tab_dir, "q10_per_plot_supplementary.csv"))
  p_s6 <- ggplot(q10_plot, aes(x = factor(Plot), y = Q10, colour = Treatment)) +
    geom_hline(yintercept = 1, linetype = 3, colour = "grey50") +
    geom_point(size = 3) +
    scale_colour_manual(values = trt_cols) +
    labs(x = "Plot", y = expression(Apparent~Q[10]),
         title = "Per-plot apparent temperature sensitivity",
         subtitle = "Partial-pooling plot slopes from an AR(1)-aware hierarchical model; descriptive estimates, not plot-level hypothesis tests.") +
    theme_pub()
  save_fig(p_s6, "FigureS6_Q10_perplot", width = 7, height = 5)
} else {
  warning("The AR(1)-aware hierarchical per-plot Q10 model did not certify; Figure S6 was not created.")
}

# ============================================================================
# SUPPLEMENTARY: Climate incident table figure (HURCN)
# ============================================================================
if (nrow(hurcn) > 0 && "Event.Description" %in% names(hurcn)) {
  plot_data <- hurcn |>
    arrange(Start) |>
    mutate(Date_Label  = format(Start, "%b %d, %Y"),
           Final_Label = str_trim(Event.Description),
           Row_ID      = row_number())
  p_table <- ggplot(plot_data) +
    geom_text(aes(x = 0.8, y = Row_ID, label = Date_Label),
              hjust = 1, size = 3.5, colour = "grey30") +
    geom_text(aes(x = 0.9, y = Row_ID, label = Final_Label),
              hjust = 0, size = 3.5, fontface = "bold", colour = "black") +
    scale_y_reverse() + scale_x_continuous(limits = c(0, 5)) +
    coord_cartesian(clip = "off") + theme_void() +
    theme(plot.margin = margin(15, 15, 15, 15))
  ggsave(file.path(fig_dir, "FigureS_ClimateTable.png"), p_table,
         width = 8, height = 1 + nrow(plot_data) * 0.35, dpi = 300)
}

cat("\nVisualisation complete. Figures in:", normalizePath(fig_dir), "\n")
