# ============================================================================
# Primary GLMM comprehensive diagnostics -- full model set (definitive)
# ----------------------------------------------------------------------------
# dispformula = ~Treatment uniformly (decided). Grand-mean centring for the
# mechanistic covariates. Scope per model class (not every section on every
# model -- only meaningful combinations):
#
#  MODEL SET (8):
#   pooled:      temp/vwc/log_flux ~ Treatment + Year_f + season
#   interaction: temp/vwc/log_flux ~ Treatment*Year_f + Treatment*season
#   mech pathway:    log_flux ~ Treatment + Year_f + season + temp_c+vwc_c+I(vwc_c^2)
#   mech modification: + Treatment:temp_c + Treatment:vwc_c + Treatment:I(vwc_c^2)
#
#  S1 CONVERGENCE CONFIRMATION (all 8, dispformula ~Treatment): pdHess, phi, key effect.
#  S2 MULTI-START REPRODUCIBILITY (all 8): dispersed nlminb starts,
#     using each model's retained dispersion formula.
#  S3 ENCODING + GAP HANDLING (3 interaction models): calendar-day ar1 vs
#     within-year full-grid ar1 vs daily OU vs observed-rank ar1; and
#     numFactor vs rank-order equivalence -- completes VWC & flux.
#  S4 PER-YEAR MODEL-VS-RAW (3 interaction models): observation-weighted and
#     plot-balanced raw contrasts.
#  S5 RESIDUAL AUTOCORR lag 1/7/14 (all final fits).
#  S6 MECHANISTIC QUANTITIES under new spec: Q10 by treatment, VWC vertex by
#     treatment, residual treatment effect under + Year_f + season.
#
# Marginals/effects: EQUAL weights (emmeans default) for marginal treatment contrasts.
# Inputs : FLUX.csv, T_VWC.csv, CLIMATE.csv
# Outputs: OUTPUT_DIAG/GLMM_diagnostics_REPORT.txt
# ============================================================================
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(lubridate)
  library(glmmTMB); library(emmeans)
})
setwd("D:/USDA/TRACE_DM_APR_26/CO2")
out_dir<-"OUTPUT_DIAG"; if(!dir.exists(out_dir)) dir.create(out_dir)
ctrl<-glmmTMBControl(optCtrl=list(iter.max=3000,eval.max=3000))
FIT_TIME_LIMIT_SECONDS<-900  # soft R time limit; compiled optimizer code may finish its current call before interruption
DISP<-~Treatment
REPORT<-file.path(out_dir,"GLMM_diagnostics_REPORT.txt")
wr<-function(...) { cat(..., file=REPORT, append=TRUE); flush.console() }
cat("PRIMARY GLMM DIAGNOSTICS -- FULL MODEL SET\nGenerated: ",format(Sys.time()),
    "\ndispformula = ~Treatment uniformly\n",strrep("=",68),"\n", file=REPORT)

clim<-read_csv("CLIMATE.csv",show_col_types=FALSE) %>%
  mutate(Year=as.integer(Year),Month=as.integer(Month),
         season=factor(season,levels=c("Dry","Wet"))) %>% select(Year,Month,season)
season_of<-function(d) left_join(tibble(Year=year(d),Month=month(d)),clim,by=c("Year","Month"))$season
mean_or_na<-function(x) if(all(is.na(x))) NA_real_ else mean(x,na.rm=TRUE)
add_time<-function(df){ gmin<-min(df$Date)
  df %>% mutate(day_num=as.integer(Date-gmin)+1L,time_num=glmmTMB::numFactor(day_num)) %>%
    group_by(Plot_Year) %>% arrange(Date,.by_group=TRUE) %>%
    mutate(rank_in_py=row_number(),day_in_py=as.integer(Date-min(Date))+1L,
           time_rank=glmmTMB::numFactor(rank_in_py),time_grid=glmmTMB::numFactor(day_in_py)) %>% ungroup() }

tvwc<-read_csv("T_VWC.csv",show_col_types=FALSE) %>%
  mutate(Date=as.Date(DayHour),Treatment=factor(Treatment,levels=c("control","warmed")),Plot=factor(Plot)) %>%
  filter(Date>=as.Date("2018-09-01")) %>%
  group_by(Date,Plot,Treatment) %>%
  summarise(temp_mean=mean_or_na(Temperature),vwc_mean=mean_or_na(VWC),.groups="drop") %>%
  mutate(vwc_pct=vwc_mean*100,Year=year(Date),Year_f=factor(Year),season=season_of(Date),
         Plot_Year=factor(paste(Plot,Year,sep="_"))) %>% filter(!is.na(season),!is.na(temp_mean),!is.na(vwc_mean)) %>% arrange(Date) %>% add_time()
flux<-read_csv("FLUX.csv",show_col_types=FALSE) %>%
  mutate(Date=as.Date(DayHour),Treatment=factor(Treatment,levels=c("control","warmed")),Plot=factor(Plot)) %>%
  filter(Date>=as.Date("2018-09-01"),!is.na(Flux),Flux>0) %>%
  group_by(Date,Plot,Treatment) %>%
  summarise(flux=mean_or_na(Flux),temp_mean=mean_or_na(Temperature),vwc_mean=mean_or_na(VWC),.groups="drop") %>%
  mutate(log_flux=log(flux),Year=year(Date),Year_f=factor(Year),season=season_of(Date),
         Plot_Year=factor(paste(Plot,Year,sep="_"))) %>% filter(!is.na(season),!is.na(log_flux),!is.na(temp_mean),!is.na(vwc_mean)) %>% arrange(Date) %>% add_time()
# grand-mean centred covariates for mechanistic models
tbar<-mean(flux$temp_mean,na.rm=TRUE); vbar<-mean(flux$vwc_mean,na.rm=TRUE)
flux<-flux %>% mutate(temp_c=temp_mean-tbar, vwc_c=vwc_mean-vbar)
wr("Centring: temp_c on ",round(tbar,3),"degC; vwc_c on ",round(vbar,4)," m3/m3\n")

get_phi<-function(m){
  vc<-VarCorr(m)$cond
  if(!"Plot_Year" %in% names(vc)) return(NA_real_)
  cm<-attr(vc[["Plot_Year"]],"correlation")
  if(!is.null(cm)&&nrow(cm)>=2) cm[1,2] else NA_real_}
disp_formula_for<-function(nm) if(identical(disp_used[[nm]],"~1")) ~1 else DISP
get_marg<-function(m){
  em<-emmeans::emmeans(m,~Treatment)
  levs<-as.character(em@grid$Treatment)
  if(!identical(levs,c("control","warmed"))){
    stop("Expected Treatment levels c('control','warmed'); got ",paste(levs,collapse=", "))
  }
  mg<-as.data.frame(summary(emmeans::contrast(em,method="revpairwise"),infer=c(TRUE,TRUE)))
  if(nrow(mg)!=1 || !all(c("estimate","SE") %in% names(mg))){
    stop("Unexpected Treatment contrast table: ",paste(names(mg),collapse=", "))
  }
  if(!is.null(mg$contrast) && !grepl("warmed.*control",mg$contrast[1])){
    stop("Expected warmed-control contrast; got '",mg$contrast[1],"'")
  }
  c(est=mg$estimate[1],se=mg$SE[1])}
fit_capped<-function(expr){
  setTimeLimit(cpu=FIT_TIME_LIMIT_SECONDS,elapsed=FIT_TIME_LIMIT_SECONDS,transient=TRUE)
  on.exit(setTimeLimit(cpu=Inf,elapsed=Inf),add=TRUE)
  tryCatch(eval.parent(substitute(expr)),error=function(e) e)}
fitm<-function(rhs,data,time_term="ar1(time_num+0|Plot_Year)",start=NULL,dispformula=DISP){
  f<-as.formula(paste(rhs,"+ (1|Plot) +",time_term))
  fit_capped(glmmTMB(f,dispformula=dispformula,data=data,REML=TRUE,start=start,control=ctrl))}

# model set
MODELS<-list(
  temp_pooled = list(rhs="temp_mean ~ Treatment + Year_f + season", data=tvwc, cls="pooled"),
  temp_inter  = list(rhs="temp_mean ~ Treatment*Year_f + Treatment*season", data=tvwc, cls="inter"),
  vwc_pooled  = list(rhs="vwc_pct ~ Treatment + Year_f + season", data=tvwc, cls="pooled"),
  vwc_inter   = list(rhs="vwc_pct ~ Treatment*Year_f + Treatment*season", data=tvwc, cls="inter"),
  flux_pooled = list(rhs="log_flux ~ Treatment + Year_f + season", data=flux, cls="pooled"),
  flux_inter  = list(rhs="log_flux ~ Treatment*Year_f + Treatment*season", data=flux, cls="inter"),
  flux_mech_path = list(rhs="log_flux ~ Treatment + Year_f + season + temp_c + vwc_c + I(vwc_c^2)", data=flux, cls="mech"),
  flux_mech_mod  = list(rhs=paste("log_flux ~ Treatment + Year_f + season + temp_c + vwc_c + I(vwc_c^2)",
                                  "+ Treatment:temp_c + Treatment:vwc_c + Treatment:I(vwc_c^2)"), data=flux, cls="mech"))

# ---- S1 convergence confirmation -------------------------------------------
wr("\n",strrep("=",68),"\nS1  CONVERGENCE CONFIRMATION (dispformula ~Treatment)\n",strrep("=",68),"\n")
fits<-list()
disp_used<-setNames(rep("~Treatment", length(MODELS)), names(MODELS))
for(nm in names(MODELS)){
  m<-fitm(MODELS[[nm]]$rhs, MODELS[[nm]]$data)
  if(inherits(m,"error")){wr(sprintf("  %-16s FAILED/TIMEOUT: %s\n",nm,conditionMessage(m)));next}
  fits[[nm]]<-m; mg<-get_marg(m)
  wr(sprintf("  %-16s pdHess=%-5s  phi=%.4f  Treatment-marginal=%.4f (SE %.4f)  AIC=%.1f\n",
             nm,isTRUE(m$sdr$pdHess),get_phi(m),mg["est"],mg["se"],AIC(m)))
  # --- fallback: if a model did not certify under ~Treatment, try disp ~1 ---
  if(!isTRUE(m$sdr$pdHess)){
    m1<-fit_capped(glmmTMB(as.formula(paste(MODELS[[nm]]$rhs,"+ (1|Plot) + ar1(time_num+0|Plot_Year)")),
                           dispformula=~1, data=MODELS[[nm]]$data, REML=TRUE, control=ctrl))
    if(!inherits(m1,"error") && isTRUE(m1$sdr$pdHess)){
      mg1<-get_marg(m1)
      wr(sprintf("  %-16s [fallback disp~1] pdHess=TRUE  phi=%.4f  Treatment-marginal=%.4f (SE %.4f)  AIC=%.1f\n",
                 nm,get_phi(m1),mg1["est"],mg1["se"],AIC(m1)))
      fits[[nm]]<-m1; disp_used[[nm]]<-"~1"
    } else {
      wr(sprintf("  %-16s [FAILED] neither ~Treatment nor ~1 certified with pdHess=TRUE; excluding from downstream diagnostics\n",nm))
      fits[[nm]]<-NULL; disp_used[[nm]]<-"FAILED"
      next
    }
  }
}

# ---- S2 multi-start ---------------------------------------------------------
wr("\n",strrep("=",68),"\nS2  MULTI-START REPRODUCIBILITY (all models)\n",strrep("=",68),"\n")
for(nm in names(MODELS)){
  base<-fits[[nm]]; if(is.null(base)){wr("  ",nm," no base fit\n");next}
  th0<-getME(base,"theta");be0<-fixef(base)$cond;bd0<-fixef(base)$disp;np<-length(th0)
  set.seed(1); res<-list()
  addfit<-function(tag,m){ if(inherits(m,"error")) return(NULL)
    if(!isTRUE(m$sdr$pdHess)) return(NULL)
    c(ll=as.numeric(logLik(m)),phi=get_phi(m),est=as.numeric(get_marg(m)["est"])) }
  b<-addfit("base",base); if(!is.null(b)) res[["base"]]<-b
  for(i in 1:5){ st<-list(beta=be0,betadisp=bd0,theta=th0+rnorm(np,0,0.5))
    m<-fitm(MODELS[[nm]]$rhs,MODELS[[nm]]$data,start=st,dispformula=disp_formula_for(nm)); a<-addfit(i,m); if(!is.null(a)) res[[paste0("s",i)]]<-a }
  if(length(res)>=2){M<-do.call(rbind,res)
    wr(sprintf("  %-16s %d pdHess fits: logLik spread=%.4f phi spread=%.5f marginal spread=%.5f  %s\n",
        nm,nrow(M),diff(range(M[,"ll"])),diff(range(M[,"phi"])),diff(range(M[,"est"])),
        ifelse(diff(range(M[,"ll"]))<0.01 && diff(range(M[,"est"]))<5e-3,"REPRODUCIBLE","CHECK")))
  } else wr(sprintf("  %-16s <2 pdHess fits; cannot judge\n",nm))
}

# ---- S3 encoding + gap handling (interaction models) -----------------------
wr("\n",strrep("=",68),"\nS3  ENCODING + GAP HANDLING (interaction models: temp, vwc, flux)\n",strrep("=",68),"\n")
encs<-c(calendar_day_ar1="ar1(time_num+0|Plot_Year)",within_year_grid_ar1="ar1(time_grid+0|Plot_Year)",
        daily_ou="ou(time_num+0|Plot_Year)",observed_rank_ar1="ar1(time_rank+0|Plot_Year)")
for(nm in c("temp_inter","vwc_inter","flux_inter")){
  wr("\n-- ",nm," --\n")
  for(en in names(encs)){
    m<-fitm(MODELS[[nm]]$rhs,MODELS[[nm]]$data,encs[en],dispformula=disp_formula_for(nm))
    if(inherits(m,"error")){wr(sprintf("  %-14s FAIL/TIMEOUT\n",en));next}
    mg<-get_marg(m)
    wr(sprintf("  %-14s pdHess=%-5s  phi=%.4f  marginal=%.4f (SE %.4f)  AIC=%.1f\n",
               en,isTRUE(m$sdr$pdHess),get_phi(m),mg["est"],mg["se"],AIC(m)))
  }
}
# structural ordering (encoding-independent) for the three data frames
wr("\n  structural ordering check (numFactor vs rank), per response frame:\n")
for(lab in c("tvwc","flux")){ d<-get(lab)
  s<-d %>% group_by(Plot_Year) %>% arrange(Date,.by_group=TRUE) %>%
    summarise(a=paste(order(as.integer(as.character(day_num))),collapse=","),
              b=paste(order(rank_in_py),collapse=","),.groups="drop") %>% mutate(id=a==b)
  wr(sprintf("    %-6s %d/%d Plot_Year groups identical ordering\n",lab,sum(s$id),nrow(s))) }

# ---- S4 per-year validation (interaction models) ---------------------------
wr("\n",strrep("=",68),"\nS4  PER-YEAR MODEL-VS-RAW (interaction models)\n",strrep("=",68),"\n")
for(nm in c("temp_inter","vwc_inter","flux_inter")){
  m<-fits[[nm]]; if(is.null(m)) next
  resp<-all.vars(as.formula(MODELS[[nm]]$rhs))[1]; data<-MODELS[[nm]]$data
  mc<-as.data.frame(summary(emmeans::contrast(emmeans::emmeans(m,~Treatment|Year_f),method="revpairwise")))
  raw_obs<-data %>% group_by(Year_f,Treatment) %>% summarise(mu=mean_or_na(.data[[resp]]),.groups="drop") %>%
    pivot_wider(names_from=Treatment,values_from=mu) %>% mutate(raw_obs_diff=warmed-control)
  raw_plot<-data %>% group_by(Year_f,Treatment,Plot) %>% summarise(plot_mu=mean_or_na(.data[[resp]]),.groups="drop") %>%
    group_by(Year_f,Treatment) %>% summarise(mu=mean_or_na(plot_mu),.groups="drop") %>%
    pivot_wider(names_from=Treatment,values_from=mu) %>% mutate(raw_plot_diff=warmed-control)
  cmp<-mc %>% select(Year_f,model_diff=estimate) %>%
    left_join(raw_obs %>% select(Year_f,raw_obs_diff),by="Year_f") %>%
    left_join(raw_plot %>% select(Year_f,raw_plot_diff),by="Year_f") %>%
    mutate(discrepancy_obs=model_diff-raw_obs_diff,discrepancy_plot=model_diff-raw_plot_diff)
  wr("\n-- ",nm," (max |disc obs|=",sprintf("%.4f",max(abs(cmp$discrepancy_obs),na.rm=TRUE)),
     "; max |disc plot|=",sprintf("%.4f",max(abs(cmp$discrepancy_plot),na.rm=TRUE)),") --\n")
  wr(paste(capture.output(print(as.data.frame(cmp %>% mutate(across(where(is.numeric),~round(.,4)))),row.names=FALSE)),collapse="\n"),"\n")
}

# ---- S5 residual autocorrelation (all final fits) --------------------------
wr("\n",strrep("=",68),"\nS5  RESIDUAL AUTOCORRELATION (within Plot_Year, calendar lag)\n",strrep("=",68),"\n")
for(nm in names(fits)){
  m<-fits[[nm]]; data<-MODELS[[nm]]$data; data$.res<-residuals(m,type="response")
  acf_at<-function(lag){ d<-data %>% select(Plot_Year,Date,.res)
    j<-d %>% mutate(Date=Date+lag)
    mm<-inner_join(d,j,by=c("Plot_Year","Date"),suffix=c("",".lag"))
    if(nrow(mm)<10) NA_real_ else cor(mm$.res,mm$.res.lag) }
  wr(sprintf("  %-16s lag1=%.3f  lag7=%.3f  lag14=%.3f\n",nm,acf_at(1),acf_at(7),acf_at(14)))
}

# ---- S6 mechanistic quantities under new spec ------------------------------
wr("\n",strrep("=",68),"\nS6  MECHANISTIC QUANTITIES (new spec: + Year_f + season)\n",strrep("=",68),"\n")
mm<-fits[["flux_mech_mod"]]
if(!is.null(mm)){
  if(!isTRUE(mm$sdr$pdHess)){
    wr("  [WARNING] flux_mech_mod did not certify with pdHess=TRUE; S6 mechanistic quantities are not reported.\n")
  } else {
  b<-fixef(mm)$cond; V<-vcov(mm)$cond
  if(!"temp_c" %in% names(b)) stop("Q10 calculation expects temp_c in degrees C")
  wr("  [Q10 assumes temp_c is in degrees C and response is natural-log flux.]\n")
  nm_t<-"temp_c"; nm_tw<-grep("Treatmentwarmed:temp_c",names(b),value=TRUE)
  nm_v<-"vwc_c"; nm_vw<-grep("Treatmentwarmed:vwc_c",names(b),value=TRUE)
  nm_q<-"I(vwc_c^2)"; nm_qw<-grep("Treatmentwarmed:I\\(vwc_c\\^2\\)",names(b),value=TRUE)
  slope_c<-b[nm_t]; slope_w<-b[nm_t]+ (if(length(nm_tw)) b[nm_tw] else 0)
  q10_c<-exp(10*slope_c); q10_w<-exp(10*slope_w)
  # slope difference (warmed-control) = the Treatment:temp_c interaction, with SE
  if(length(nm_tw)){ dslope<-b[nm_tw]; dse<-sqrt(V[nm_tw,nm_tw])
    wr(sprintf("  temp slope diff (warmed-control)=%.5f (SE %.5f, p=%.3f)\n",
               dslope,dse,2*pnorm(-abs(dslope/dse)))) }
  wr(sprintf("  Q10 control=%.3f  warmed=%.3f  (slope_c=%.5f slope_w=%.5f)\n",q10_c,q10_w,slope_c,slope_w))
  # VWC optimum only valid as a MAX when quadratic coef < 0; report sign
  b1c<-b[nm_v]; b2c<-b[nm_q]
  b1w<-b[nm_v]+(if(length(nm_vw)) b[nm_vw] else 0); b2w<-b[nm_q]+(if(length(nm_qw)) b[nm_qw] else 0)
  vertex<-function(b1,b2) (-b1/(2*b2))+vbar
  wr(sprintf("  quadratic coef: control b2=%.4f (%s)  warmed b2=%.4f (%s)\n",
             b2c, ifelse(b2c<0,"concave: vertex=MAX","convex: vertex=MIN"),
             b2w, ifelse(b2w<0,"concave: vertex=MAX","convex: vertex=MIN")))
  wr(sprintf("  VWC vertex control=%.4f  warmed=%.4f  m3/m3 (interpret per sign above)\n",
             vertex(b1c,b2c),vertex(b1w,b2w)))
  wr(sprintf("  residual Treatment main effect at mean T,VWC (log): %.4f\n",
             b[grep("^Treatmentwarmed$",names(b))]))
  }
} else wr("  mechanistic modification model did not converge\n")

wr("\n",strrep("=",68),"\nSUMMARY: full model set, dispformula ~Treatment, encoding+gap for temp/vwc/flux,\n")
wr("  dispformula used per model: ",
   paste(sprintf("%s=%s",names(disp_used),unlist(disp_used)),collapse="; "),"\n")
wr("mechanistic quantities under new spec.\n")
writeLines(capture.output(sessionInfo()),file.path(out_dir,"sessionInfo_diag.txt"))
cat("Done. Report at",REPORT,"\n")
