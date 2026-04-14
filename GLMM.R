# ============================================================================
# COMPREHENSIVE PUBLICATION-READY ANALYSIS SCRIPT
# Soil CO2 Flux Warming Experiment
# ============================================================================

# REQUIRED PACKAGES
library(tidyverse)
library(glmmTMB)
library(DHARMa)
library(emmeans)
library(multcomp)
library(MuMIn)
library(ggeffects)
library(patchwork)
library(flextable)
library(officer)
library(performance)
library(viridis)

# Setup directories
setwd("D:/USDA/TRACE_DM_DEC/STAT")
out_dir <- "OUTPUT_GLMM"
dir.create(out_dir, showWarnings = FALSE)
dir.create(file.path(out_dir, "figures"), showWarnings = FALSE)
dir.create(file.path(out_dir, "tables"), showWarnings = FALSE)
dir.create(file.path(out_dir, "diagnostics"), showWarnings = FALSE)

# ============================================================================
# PART 1: DATA LOADING & PROCESSING
# ============================================================================

# 1) CLIMATE DATA & SEASON FACTOR CALCULATION
cat("\n>>> PROCESSING CLIMATE DATA (CWD ONLY)...\n")

# Settings
baseline_start <- 2000
baseline_end   <- 2024

# Load & Process Raw Climate Data
monthly_clim_raw <- read.csv("D:/USDA/Data/CLIMATE_DATA/terraclimate_loq.csv") %>%
  mutate(
    Year = as.integer(Year),
    Month = as.integer(Month),
    precip_mm = as.numeric(ppt.mm.),
    cwd_mm    = as.numeric(def.mm.) # Climatic Water Deficit
  ) %>%
  arrange(Year, Month)

# --- CALCULATE THRESHOLD (Type 8 Median) ---
baseline_data <- monthly_clim_raw %>% 
  filter(Year >= baseline_start, Year <= baseline_end)

thresh_CWD_med <- as.numeric(quantile(baseline_data$cwd_mm, 0.5, na.rm=TRUE, type=8, names=FALSE))

cat(sprintf("CWD Median Threshold (Type 8): %.4f mm\n", thresh_CWD_med))

# --- CONSTRUCT 'Season' FACTOR ---
climate_key <- monthly_clim_raw %>%
  mutate(
    # The Consensus Choice: CWD Median Split
    Season = factor(ifelse(cwd_mm > thresh_CWD_med, "Dry", "Wet"), 
                    levels = c("Dry", "Wet"))
  ) %>%
  dplyr::select(Year, Month, precip_monthly = precip_mm, cwd_monthly = cwd_mm, Season)

# 2) FLUX DATA PREPARATION
cat("\n>>> LOADING & AGGREGATING FLUX DATA...\n")

raw_data <- read.csv("D:/USDA/TRACE_DM_DEC/STAT/FLUX.csv", stringsAsFactors = FALSE)

d_daily <- raw_data %>%
  mutate(
    DayHour = lubridate::ymd_hms(DayHour, tz = "UTC"), 
    Date = as.Date(DayHour),
    Year = year(DayHour), 
    Month = month(DayHour), 
    DOY = yday(DayHour),
    Plot = factor(Plot), 
    Treatment = factor(Treatment, levels = c("control", "warmed"))
  ) %>%
  filter(Flux >= 1e-6, DayHour >= ymd_hms("2018-09-27T00:00:00Z")) %>%
  # Join Climate Key
  left_join(climate_key, by = c("Year", "Month")) %>%
  group_by(Plot, Treatment, Year, Month, DOY, Date) %>%
  summarise(
    Flux = mean(Flux, na.rm = TRUE), 
    Temperature = mean(Temperature, na.rm = TRUE),
    VWC = mean(VWC, na.rm = TRUE), 
    precip_monthly = first(precip_monthly),
    cwd_monthly = first(cwd_monthly),
    Season = first(Season), # Carry the Season factor forward
    n_hours = n(), 
    .groups = "drop"
  ) %>%
  filter(n_hours >= 12)

# Stats summary
stats_clim <- d_daily %>%
  summarise(
    temp_mean = mean(Temperature, na.rm = TRUE), temp_sd = sd(Temperature, na.rm = TRUE),
    vwc_mean = mean(VWC, na.rm = TRUE), vwc_sd = sd(VWC, na.rm = TRUE),
    precip_mean = mean(precip_monthly, na.rm = TRUE), precip_sd = sd(precip_monthly, na.rm = TRUE),
    flux_mean = mean(Flux, na.rm = TRUE), flux_sd = sd(Flux, na.rm = TRUE)
  )

cat("\n>>> GENERATING LAGS, SCALING, AND GAP_LOG...\n")


d_complete <- d_daily %>%
  # 1. Base transformations
  mutate(
    log_flux = log(Flux),
    sin_doy  = sin(2 * pi * DOY / 365), 
    cos_doy  = cos(2 * pi * DOY / 365),
    Year_c   = Year - 2021
  ) %>%
  arrange(Plot, Date) %>%
  
  # 2. Generate RAW Lags (Per Plot, Unscaled)
  group_by(Plot) %>%
  mutate(
    Gap_Days = as.numeric(Date - lag(Date, default = first(Date))),
    Gap_log  = log1p(Gap_Days),
    
    # Raw Lags - Flux
    lag_flux_1  = lag(log_flux, 1),
    lag_flux_3  = lag(log_flux, 3),
    lag_flux_7  = lag(log_flux, 7),
    lag_flux_14 = lag(log_flux, 14),
    lag_flux_30 = lag(log_flux, 30),
    
    # Raw Lags - Temp
    lag_temp_1  = lag(Temperature, 1),
    lag_temp_3  = lag(Temperature, 3),
    lag_temp_7  = lag(Temperature, 7),
    lag_temp_14 = lag(Temperature, 14),
    lag_temp_30 = lag(Temperature, 30),
    
    # Raw Lags - VWC
    lag_vwc_1   = lag(VWC, 1),
    lag_vwc_3   = lag(VWC, 3),
    lag_vwc_7   = lag(VWC, 7),
    lag_vwc_14  = lag(VWC, 14),
    lag_vwc_30  = lag(VWC, 30)
  ) %>%
  ungroup() %>% # <--- CRITICAL STEP: Ungroup before scaling
  
  # 3. Global Scaling (Compare all plots on the same ruler)
  mutate(
    # Scale Main Drivers
    temp_c   = as.numeric(scale(Temperature)),
    vwc_c    = as.numeric(scale(VWC)),
    precip_c = as.numeric(scale(precip_monthly)),
    
    # Scale Flux Lags
    Flux_lag1_c  = as.numeric(scale(lag_flux_1)),
    Flux_lag3_c  = as.numeric(scale(lag_flux_3)),
    Flux_lag7_c  = as.numeric(scale(lag_flux_7)),
    Flux_lag14_c = as.numeric(scale(lag_flux_14)),
    Flux_lag30_c = as.numeric(scale(lag_flux_30)),
    
    # Scale Temp Lags
    Temp_lag1_c  = as.numeric(scale(lag_temp_1)),
    Temp_lag3_c  = as.numeric(scale(lag_temp_3)),
    Temp_lag7_c  = as.numeric(scale(lag_temp_7)),
    Temp_lag14_c = as.numeric(scale(lag_temp_14)),
    Temp_lag30_c = as.numeric(scale(lag_temp_30)),
    
    # Scale VWC Lags
    VWC_lag1_c   = as.numeric(scale(lag_vwc_1)),
    VWC_lag3_c   = as.numeric(scale(lag_vwc_3)),
    VWC_lag7_c   = as.numeric(scale(lag_vwc_7)),
    VWC_lag14_c  = as.numeric(scale(lag_vwc_14)),
    VWC_lag30_c  = as.numeric(scale(lag_vwc_30))
  ) %>%
  # Remove temporary raw lag columns to keep dataframe clean
  dplyr::select(-starts_with("lag_flux_"), -starts_with("lag_temp_"), -starts_with("lag_vwc_"))

# FINAL OUTPUT SPLITS
d_flux_final <- d_complete %>% 
  filter(!is.na(Flux_lag30_c))
d_flux_final <- d_flux_final %>% 
  drop_na(vwc_c, temp_c)
d_temp_final <- d_complete %>% 
  filter(!is.na(Temp_lag7_c))%>% 
  drop_na(vwc_c, temp_c)
d_vwc_final <- d_complete %>% 
  filter(!is.na(VWC_lag30_c))%>% 
  drop_na(vwc_c, temp_c)

# ============================================================================
# PART 2: FIT MODELS
# ============================================================================

cat("\n>>> FITTING MODELS...\n")
ctrl <- glmmTMBControl(optCtrl = list(iter.max = 20000, eval.max = 20000), 
                       optimizer = "nlminb")

glmm_TOTAL <- glmmTMB(
  log_flux ~ Treatment * Season + Treatment * Year_c +
    Flux_lag1_c * Gap_log + Flux_lag3_c + Flux_lag7_c + Flux_lag14_c + Flux_lag30_c +
    (1 | Plot) + (1 | Plot:Year),
  dispformula = ~ Year_c + Treatment + temp_c + vwc_c + sin_doy,
  data = d_flux_final,
  family = gaussian(link = "identity"),
  control = ctrl)
summary(glmm_TOTAL)

glmm_MECH <- glmmTMB(
  log_flux ~ Treatment * Season + Treatment * Year_c + 
    Treatment * poly(vwc_c, 2) * temp_c + 
    Flux_lag1_c * Gap_log + Flux_lag3_c + Flux_lag7_c + Flux_lag14_c + Flux_lag30_c + 
    (1 | Plot) + (1 | Plot:Year),
  dispformula = ~ Year_c + Treatment + temp_c + vwc_c + sin_doy,
  data = d_flux_final,
  family = gaussian(link = "identity"),
  control = ctrl)
summary(glmm_MECH)

glmm_TEMP <- glmmTMB(
  Temperature ~ Treatment*Year_c + precip_c + Season + vwc_c +
    Temp_lag1_c*Gap_log + Temp_lag3_c + Temp_lag7_c +
    (1|Plot) + (1|Plot:Year),
  dispformula = ~ vwc_c + Year_c + Treatment + precip_c, 
  data = d_temp_final
)
summary(glmm_TEMP)

glmm_VWC <- glmmTMB(
  VWC ~ Treatment*Season + Year_c + Treatment*poly(temp_c, 2) + 
    VWC_lag1_c*Gap_log + VWC_lag3_c + VWC_lag30_c +
    (1|Plot) + (1|Plot:Year),
  dispformula = ~ vwc_c + Year_c + Season + Treatment + cos_doy + precip_c, 
  data = d_vwc_final
)
summary(glmm_VWC)

# ============================================================================
# PART 3: HELPER FUNCTIONS
# ============================================================================

theme_publication <- function(base_size = 20) {
  theme_bw(base_size = base_size) +
    theme(axis.text = element_text(size = 20, color = "black"),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(color = "gray90", size = 0.3),
          axis.title = element_text(face = "bold"),
          legend.position = "right",
          legend.title = element_text(face = "bold"),
          strip.background = element_rect(fill = "gray95", color = "black"),
          strip.text = element_text(face = "bold", size = base_size - 1)
    )
}

create_model_table <- function(model, title_text) {
  sum_tbl <- as.data.frame(summary(model)$coefficients$cond)
  sum_tbl$Parameter <- rownames(sum_tbl)
  sum_tbl <- sum_tbl %>%
    dplyr::select(Parameter, Estimate, `Std. Error`, `z value`, `Pr(>|z|)`) %>%
    mutate(
      Estimate = round(Estimate, 4),
      `Std. Error` = round(`Std. Error`, 4),
      `z value` = round(`z value`, 3),
      `p-value` = format.pval(`Pr(>|z|)`, digits = 3, eps = 0.001),
      Significance = case_when(
        `Pr(>|z|)` < 0.001 ~ "***",
        `Pr(>|z|)` < 0.01 ~ "**",
        `Pr(>|z|)` < 0.05 ~ "*",
        `Pr(>|z|)` < 0.1 ~ ".",
        TRUE ~ ""
      )
    ) %>%
    dplyr::select(-`Pr(>|z|)`)
  
  ft <- flextable(sum_tbl) %>%
    set_caption(caption = title_text) %>%
    theme_booktabs() %>%
    autofit() %>%
    bold(part = "header") %>%
    align(align = "center", part = "all") %>%
    align(j = "Parameter", align = "left", part = "body")
  
  return(ft)
}

create_contrast_table <- function(contrast_df, title_text) {
  contrast_df <- contrast_df %>%
    mutate(across(where(is.numeric), ~round(., 4)))
  
  ft <- flextable(contrast_df) %>%
    set_caption(caption = title_text) %>%
    theme_booktabs() %>%
    autofit() %>%
    bold(part = "header") %>%
    align(align = "center", part = "all")
  
  return(ft)
}


# ============================================================================
# PART 4: COMPREHENSIVE DIAGNOSTICS
# ============================================================================

cat("\n>>> RUNNING DIAGNOSTICS...\n")

run_diagnostics <- function(model, model_name, data_df, time_var) {
  
  sim_resid <- simulateResiduals(model, n = 1000, plot = FALSE)
  
  # === uniformity ===
  uni <- testUniformity(sim_resid)
  
  # === dispersion ===
  disp <- testDispersion(sim_resid)
  
  # === temporal autocorrelation: per-plot loop ===
  temp_ac_list <- lapply(split(data_df, data_df$Plot), function(df_sub) {
    
    idx <- which(data_df$Plot == unique(df_sub$Plot))
    sub_resid <- recalculateResiduals(sim_resid, sel = idx)
    
    tvals <- df_sub[[time_var]]
    
    # skip if non-unique
    if (anyDuplicated(tvals)) {
      return(NA_real_)
    }
    
    tac <- testTemporalAutocorrelation(sub_resid, time = tvals)
    tac$p.value
  })
  
  # summarize across plots
  temp_ac_vals <- unlist(temp_ac_list)
  temp_ac_min <- suppressWarnings(min(temp_ac_vals, na.rm = TRUE))
  
  # === R2 ===
  r2_vals <- suppressWarnings(r.squaredGLMM(model))
  
  # === ICC ===
  icc_val <- performance::icc(model)
  
  # === diagnostics PDF ===
  pdf(file.path(out_dir, "diagnostics", paste0(model_name, "_dharma.pdf")),
      width = 12, height = 6)
  plot(sim_resid, rank = TRUE)
  dev.off()
  
  # === return ===
  dplyr::tibble(
    model = model_name,
    ks_p = uni$p.value,
    disp_p = disp$p.value,
    temp_ac_min_p = temp_ac_min,
    r2_marg = r2_vals[1, "R2m"],
    r2_cond = r2_vals[1, "R2c"],
    icc = suppressWarnings({
      icc_tmp <- performance::icc(model)
      if (inherits(icc_tmp, "data.frame")) {
        icc_tmp$ICC_conditional
      } else {
        NA_real_
      }
    })
  )
}


diag_list <- list(
  run_diagnostics(glmm_TOTAL, "TOTAL", d_flux_final, "DOY"),
  run_diagnostics(glmm_MECH, "MECH", d_flux_final, "DOY"),
  run_diagnostics(glmm_TEMP, "TEMP", d_temp_final, "DOY"),
  run_diagnostics(glmm_VWC, "VWC", d_vwc_final, "DOY")
)

diag_df <- dplyr::bind_rows(diag_list)
write.csv(diag_df, file.path(out_dir, "tables", "diagnostics.csv"), row.names = FALSE)

# ============================================================================
# PART 5: COMPLETE EMMEANS & CONTRASTS
# ============================================================================

# ============================================================================
# COMMON CONDITIONING VALUES (EMPIRICAL MEDIANS / TYPICAL CONDITIONS)
# ============================================================================
typical <- list(
  Season = "Dry",
  Year_c = 0,
  # centered/scaled continuous: use empirical medians from the analysis dataset
  precip_c = median(d_flux_final$precip_c, na.rm = TRUE),
  temp_c   = median(d_flux_final$temp_c,   na.rm = TRUE),
  vwc_c    = median(d_flux_final$vwc_c,    na.rm = TRUE),
  Gap_log  = median(d_flux_final$Gap_log,  na.rm = TRUE),
  # lag terms (important: use medians, not 0)
  Flux_lag1_c  = median(d_flux_final$Flux_lag1_c,  na.rm = TRUE),
  Flux_lag3_c  = median(d_flux_final$Flux_lag3_c,  na.rm = TRUE),
  Flux_lag7_c  = median(d_flux_final$Flux_lag7_c,  na.rm = TRUE),
  Flux_lag14_c = median(d_flux_final$Flux_lag14_c, na.rm = TRUE),
  Flux_lag30_c = median(d_flux_final$Flux_lag30_c, na.rm = TRUE)
)
cond_without <- function(cond, vars) {
  cond[!names(cond) %in% vars]
}

cat("\n>>> RUNNING EMMEANS & CONTRASTS...\n")

# --- TOTAL FLUX MODEL ---
cat("  Processing TOTAL model...\n")

# Treatment x Season
emm_trt_season_total <- emmeans(
  glmm_TOTAL,
  ~ Treatment | Season,
  at = list(
    Year_c     = typical$Year_c,
    precip_c   = typical$precip_c,
    temp_c     = typical$temp_c,
    vwc_c      = typical$vwc_c,
    Gap_log    = typical$Gap_log,
    Flux_lag1_c  = typical$Flux_lag1_c,
    Flux_lag3_c  = typical$Flux_lag3_c,
    Flux_lag7_c  = typical$Flux_lag7_c,
    Flux_lag14_c = typical$Flux_lag14_c,
    Flux_lag30_c = typical$Flux_lag30_c
  )
)
pairs_trt_season_total <- pairs(emm_trt_season_total, adjust = "tukey")

# Treatment marginal
emm_trt_total <- emmeans(glmm_TOTAL, ~ Treatment)
pairs_trt_total <- pairs(emm_trt_total, adjust = "tukey")

# Save all
write.csv(as.data.frame(summary(emm_trt_season_total)), 
          file.path(out_dir, "tables", "TOTAL_emmeans_TrtBySeason.csv"), row.names = FALSE)
write.csv(as.data.frame(summary(pairs_trt_season_total)), 
          file.path(out_dir, "tables", "TOTAL_contrasts_TrtWithinSeason.csv"), row.names = FALSE)
write.csv(as.data.frame(summary(pairs_trt_total)), 
          file.path(out_dir, "tables", "TOTAL_contrasts_TrtMarginal.csv"), row.names = FALSE)

# --- MECH FLUX MODEL ---
cat("  Processing MECH model...\n")

emm_trt_season_mech <- emmeans(
  glmm_MECH,
  ~ Treatment | Season,
  at = list(
    Year_c     = typical$Year_c,
    precip_c   = typical$precip_c,
    temp_c     = typical$temp_c,
    vwc_c      = typical$vwc_c,
    Gap_log    = typical$Gap_log,
    Flux_lag1_c  = typical$Flux_lag1_c,
    Flux_lag3_c  = typical$Flux_lag3_c,
    Flux_lag7_c  = typical$Flux_lag7_c,
    Flux_lag14_c = typical$Flux_lag14_c,
    Flux_lag30_c = typical$Flux_lag30_c
  )
)
pairs_trt_season_mech <- pairs(emm_trt_season_mech, adjust = "tukey")

# Emtrends at multiple focal moisture points (captures nonlinearity locally)
vwc_focal <- as.numeric(quantile(d_flux_final$vwc_c, probs = c(0.1, 0.5, 0.9), na.rm = TRUE))

emtrends_vwc_mech_focal <- emtrends(
  glmm_MECH,
  ~ Treatment | Season * vwc_c,
  var = "vwc_c",
  at = list(vwc_c = vwc_focal, temp_c = typical$temp_c)
)

pairs_trends_vwc_focal <- pairs(emtrends_vwc_mech_focal, by = c("Season", "vwc_c"), adjust = "tukey")

write.csv(as.data.frame(summary(emtrends_vwc_mech_focal)),
          file.path(out_dir, "tables", "MECH_emtrends_VWC_focal.csv"),
          row.names = FALSE)
write.csv(as.data.frame(summary(pairs_trends_vwc_focal)),
          file.path(out_dir, "tables", "MECH_contrasts_VWC_slopes_focal.csv"),
          row.names = FALSE)

write.csv(as.data.frame(summary(emm_trt_season_mech)), 
          file.path(out_dir, "tables", "MECH_emmeans_TrtBySeason.csv"), row.names = FALSE)
write.csv(as.data.frame(summary(pairs_trt_season_mech)), 
          file.path(out_dir, "tables", "MECH_contrasts_TrtWithinSeason.csv"), row.names = FALSE)

# --- VWC MODEL ---
cat("  Processing VWC model...\n")

emm_trt_season_vwc <- emmeans(glmm_VWC, ~ Treatment | Season)
emm_trt_vwc <- emmeans(glmm_VWC, ~ Treatment)

write.csv(as.data.frame(summary(emm_trt_season_vwc)), 
          file.path(out_dir, "tables", "VWC_emmeans_TrtBySeason.csv"), row.names = FALSE)

