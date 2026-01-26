############################################################
# thesis_commodity_credit_pipeline.R
# Goal: Produce TODAY (1) evidence commodity shocks affect sovereign credit,
#       (2) evidence is stronger when export exposure is higher,
#       (3) tail/state-dependence -> later motivates jump intensity / jumps.
#
# Inputs:
#   - data_cleaning.R  -> creates final_dt (weekly wide)
#   - data/Oil_copper_exposure.csv  (annual 2013-2024 exposures from UN Comtrade)
#
# Outputs (writes to ./output):
#   - per_country_results.csv                (betas + crash betas + delta test)
#   - tail_prob_results.csv                  (P(tail | crash) ratios)
#   - exposure_beta_oil.csv / exposure_beta_copper.csv
#   - pooled_models.txt                      (key pooled interaction regressions)
#   - per_country_summaries.txt              (country regression coefficient tables)
#   - eventstudy_oil.png / eventstudy_copper.png
#   - robustness_pre2025.txt
############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(lmtest)
  library(sandwich)
  library(fixest)
  library(car)
})

#========================
# 0) USER SETTINGS
#========================
OUT_DIR <- "output"
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR)

# crash threshold (use 0.10 for power; you can also run 0.05)
SHOCK_Q <- 0.10
TAIL_Q  <- 0.05     # tail event in EMBI returns
NW_LAG  <- 4        # weekly -> 4 weeks

EVENT_WINDOW <- -2:2

# Country pool
OIL_SET    <- c("Angola","Nigeria","Ecuador","Colombia")
COPPER_SET <- c("Chile","Peru","Zambia")
CTRL_SET   <- c("Brazil","Egypt","Morocco")
COUNTRIES  <- unique(c(OIL_SET, COPPER_SET, CTRL_SET))

# Column names in final_dt (your Bloomberg export style)
# IMPORTANT: use EXACT names as in final_dt
EMBI_MAP <- data.table(
  country  = c("Angola","Nigeria","Ecuador","Chile","Peru","Morocco","Egypt","Brazil","Colombia","Zambia","Mexico"),
  embi_col = c("JPGCAO_Index","JPGCNG_Index","JPGCEC_Index","JPGCCH_Index","JPGCPR_Index",
               "JPGCMOR_Index","JPGCEG_Index","JPGCBZ_Index","JPGCCO_Index","JPGCZM_Index","JPGCMX_Index")
)

# FX columns (if missing, we auto-create a dummy = 1)
FX_MAP <- data.table(
  country = c("Angola","Nigeria","Ecuador","Chile","Peru","Morocco","Egypt","Brazil","Colombia","Zambia"),
  fx_col  = c("USDAOA_Curncy","USDNGN_Curncy","ECU_FX_DUMMY","USDCLP_Curncy","USDPEN_Curncy",
              "USDMAD_Curncy","USDEGP_Curncy","USDBRL_Curncy","USDCOP_Curncy","USDZMW_Curncy")
)

GLOBAL_COLS <- c("date","USGG10YR_Index","SPX_Index","VIXCLS","CO1_Comdty","LP1_Comdty")

#========================
# 1) LOAD DATA
#========================
source("data_cleaning.R")
dt <- copy(final_dt)
setDT(dt)

# Exposure file (annual)
exp <- fread("data/Oil_copper_exposure.csv")
setDT(exp)
setnames(exp, c("Year","country"), c("year","country"), skip_absent = TRUE)
exp[, country := trimws(country)]
exp[, year := as.integer(year)]

for (cc in intersect(c("Oil exposure","Copper exposure"), names(exp))) {
  exp[, (cc) := as.numeric(gsub(",", "", as.character(get(cc))))]
}
exp_use <- exp[, .(country, year,
                   oil_exp = `Oil exposure`,
                   copper_exp = `Copper exposure`)]

#========================
# 2) CLEAN BLOOMBERG EXPORT ARTIFACTS
#========================
bad_col_pattern <- "Requesting Data|#N/A|N/A"
bad_cols <- names(dt)[grepl(bad_col_pattern, names(dt), ignore.case = TRUE)]
if (length(bad_cols) > 0) dt[, (bad_cols) := NULL]

dt[, date := as.Date(date)]
setorder(dt, date)

# Ecuador FX dummy
dt[, ECU_FX_DUMMY := 1]

#========================
# 3) SUBSET MAPS TO YOUR COUNTRY POOL
#========================
embi_map <- EMBI_MAP[country %in% COUNTRIES]
fx_map   <- FX_MAP[country %in% COUNTRIES]

# create dummy FX if missing
for (fx in fx_map$fx_col) {
  if (!fx %in% names(dt)) {
    message("FX column missing: ", fx, " -> creating dummy=1")
    dt[, (fx) := 1]
  }
}

# checks
miss_global <- setdiff(GLOBAL_COLS, names(dt))
if (length(miss_global) > 0) stop("Missing global columns:\n", paste(miss_global, collapse=", "))

miss_embi <- setdiff(embi_map$embi_col, names(dt))
if (length(miss_embi) > 0) stop("Missing EMBI index columns:\n", paste(miss_embi, collapse=", "))

#========================
# 4) WIDE -> LONG PANEL (index levels -> returns)
#========================
need_cols <- unique(c(GLOBAL_COLS, embi_map$embi_col, fx_map$fx_col))
panel_w <- dt[, ..need_cols]

embi_long <- melt(panel_w,
                  id.vars = GLOBAL_COLS,
                  measure.vars = embi_map$embi_col,
                  variable.name = "embi_col",
                  value.name = "embi_level")
embi_long <- merge(embi_long, embi_map, by="embi_col", all.x=TRUE)

fx_long <- melt(panel_w,
                id.vars = "date",
                measure.vars = fx_map$fx_col,
                variable.name = "fx_col",
                value.name = "fx_usdxxx")
fx_long <- merge(fx_long, fx_map, by="fx_col", all.x=TRUE)

panel <- merge(
  embi_long[, .(date,country,embi_level,USGG10YR_Index,SPX_Index,VIXCLS,CO1_Comdty,LP1_Comdty)],
  fx_long[, .(date,country,fx_usdxxx)],
  by=c("date","country"),
  all.x=TRUE
)

panel <- panel[country %in% COUNTRIES]
setorder(panel, country, date)

panel[, embi_level := as.numeric(embi_level)]
panel <- panel[!is.na(embi_level) & embi_level > 0]

#========================
# 5) RETURNS + CONTROLS
#========================
panel[, r_embi := log(embi_level) - log(shift(embi_level, 1)), by=country]
panel[, r_fx   := log(fx_usdxxx)  - log(shift(fx_usdxxx, 1)), by=country]

panel[, r_spx  := log(SPX_Index)  - log(shift(SPX_Index, 1)), by=country]
panel[, r_vix  := log(VIXCLS)     - log(shift(VIXCLS, 1)), by=country]
panel[, d_us10y := USGG10YR_Index - shift(USGG10YR_Index, 1), by=country]

panel[, r_oil   := log(CO1_Comdty) - log(shift(CO1_Comdty, 1)), by=country]
panel[, r_copper:= log(LP1_Comdty) - log(shift(LP1_Comdty, 1)), by=country]

panel <- panel[!is.na(r_embi) & !is.na(r_oil) & !is.na(r_copper) &
                 !is.na(r_spx) & !is.na(r_vix) & !is.na(d_us10y)]

panel[, year := as.integer(format(date, "%Y"))]

#========================
# 6) MERGE EXPOSURES (2013-2024) + FILL 2025
#========================
panel <- merge(panel, exp_use, by=c("country","year"), all.x=TRUE)

setorder(panel, country, date)
panel[, oil_exp := nafill(oil_exp, type="locf"), by=country]
panel[, oil_exp := nafill(oil_exp, type="nocb"), by=country]
panel[, copper_exp := nafill(copper_exp, type="locf"), by=country]
panel[, copper_exp := nafill(copper_exp, type="nocb"), by=country]
panel[, exp_is_imputed := as.integer(year == 2025)]

#========================
# 7) DEFINE CRASH WEEKS
#========================
u <- unique(panel[, .(date, r_oil, r_copper)])
q_oil <- quantile(u$r_oil, SHOCK_Q, na.rm=TRUE)
q_cu  <- quantile(u$r_copper, SHOCK_Q, na.rm=TRUE)

panel[, oil_crash := as.integer(r_oil < q_oil)]
panel[, copper_crash := as.integer(r_copper < q_cu)]

#========================
# 8) PER-COUNTRY REGRESSIONS (baseline + crash slope) + DELTA TEST
#========================
run_one_country <- function(cty) {
  d <- panel[country == cty]
  
  # choose driver + crash variable
  if (cty %in% COPPER_SET) {
    x <- "r_copper"; crash <- "copper_crash"
  } else if (cty %in% OIL_SET) {
    x <- "r_oil"; crash <- "oil_crash"
  } else {
    x <- "r_oil"; crash <- "oil_crash"
  }
  
  # If crash has (almost) no variation after filtering, skip crash interaction
  n_crash <- sum(d[[crash]] == 1, na.rm = TRUE)
  n_total <- nrow(d)
  crash_ok <- is.finite(n_crash) && n_crash >= 10 && n_crash <= (n_total - 10)
  
  # -------------------------
  # Baseline regression
  # -------------------------
  f1 <- as.formula(paste0("r_embi ~ ", x, " + r_spx + r_vix + d_us10y + r_fx"))
  m1 <- lm(f1, data = d)
  V1 <- NeweyWest(m1, lag = NW_LAG, prewhite = FALSE)
  ct1 <- coeftest(m1, vcov. = V1)
  
  beta <- ct1[x, 1]
  pval <- ct1[x, 4]
  
  # -------------------------
  # Crash interaction regression (safe)
  # r_embi ~ x + x*crash + controls
  # -------------------------
  if (crash_ok) {
    d[, x_crash_int := get(x) * get(crash)]
    
    f2 <- as.formula(paste0("r_embi ~ ", x, " + x_crash_int + r_spx + r_vix + d_us10y + r_fx"))
    m2 <- lm(f2, data = d)
    V2 <- NeweyWest(m2, lag = NW_LAG, prewhite = FALSE)
    ct2 <- coeftest(m2, vcov. = V2)
    
    # normal slope = beta_x
    beta_norm <- ct2[x, 1]
    p_norm    <- ct2[x, 4]
    
    # delta = interaction coefficient
    # crash slope = beta_norm + delta
    if ("x_crash_int" %in% rownames(ct2)) {
      delta <- ct2["x_crash_int", 1]
      p_delta <- ct2["x_crash_int", 4]
      beta_crash <- beta_norm + delta
      
      # p-value for crash slope itself isn't directly ct2; keep NA (or compute later if needed)
      p_crash <- NA_real_
    } else {
      # interaction got dropped => cannot identify delta
      delta <- NA_real_
      p_delta <- NA_real_
      beta_crash <- NA_real_
      p_crash <- NA_real_
    }
  } else {
    beta_norm <- NA_real_
    p_norm <- NA_real_
    beta_crash <- NA_real_
    p_crash <- NA_real_
    delta <- NA_real_
    p_delta <- NA_real_
  }
  
  data.table(
    country = cty,
    driver = x,
    n = n_total,
    n_crash = n_crash,
    
    beta = as.numeric(beta),
    p = as.numeric(pval),
    
    beta_norm = as.numeric(beta_norm),
    p_norm = as.numeric(p_norm),
    
    beta_crash = as.numeric(beta_crash),
    p_crash = as.numeric(p_crash),
    
    delta_crash_minus_norm = as.numeric(delta),
    p_delta = as.numeric(p_delta)
  )
}

res_main <- rbindlist(lapply(COUNTRIES, run_one_country))
res_main[, group := fifelse(country %in% OIL_SET, "Oil exporter",
                            fifelse(country %in% COPPER_SET, "Copper exporter", "Control"))]
setorder(res_main, group, p)

fwrite(res_main, file.path(OUT_DIR, "per_country_results.csv"))

# Controls: run copper too (so you have placebo copper)
run_control_driver <- function(cty, x, crash) {
  d <- panel[country == cty]
  f1 <- as.formula(paste0("r_embi ~ ", x, " + r_spx + r_vix + d_us10y + r_fx"))
  m1 <- lm(f1, data=d)
  ct1 <- coeftest(m1, vcov.=NeweyWest(m1, lag=NW_LAG, prewhite=FALSE))
  
  d[, x_norm := get(x) * (1 - get(crash))]
  d[, x_crash := get(x) * get(crash)]
  m2 <- lm(r_embi ~ x_norm + x_crash + r_spx + r_vix + d_us10y + r_fx, data=d)
  V2 <- NeweyWest(m2, lag=NW_LAG, prewhite=FALSE)
  ct2 <- coeftest(m2, vcov.=V2)
  lh <- linearHypothesis(m2, "x_crash = x_norm", vcov. = V2)
  
  data.table(
    country=cty, driver=x,
    beta = ct1[x,1], p = ct1[x,4],
    beta_norm = ct2["x_norm",1], p_norm = ct2["x_norm",4],
    beta_crash= ct2["x_crash",1], p_crash= ct2["x_crash",4],
    delta = coef(m2)["x_crash"] - coef(m2)["x_norm"],
    p_delta = lh$`Pr(>F)`[2]
  )
}

ctrl_oil <- rbindlist(lapply(CTRL_SET, run_control_driver, x="r_oil", crash="oil_crash"))
ctrl_cu  <- rbindlist(lapply(CTRL_SET, run_control_driver, x="r_copper", crash="copper_crash"))
fwrite(ctrl_oil, file.path(OUT_DIR, "controls_oil_results.csv"))
fwrite(ctrl_cu,  file.path(OUT_DIR, "controls_copper_results.csv"))

# Save coefficient tables in one text file
sink(file.path(OUT_DIR, "per_country_summaries.txt"))
cat("=== PER-COUNTRY RESULTS (sorted) ===\n\n")
print(res_main[order(group, p)])
cat("\n\n=== CONTROLS: Oil ===\n\n")
print(ctrl_oil[order(p)])
cat("\n\n=== CONTROLS: Copper ===\n\n")
print(ctrl_cu[order(p)])
sink()

#========================
# 9) TAIL PROBABILITY TEST (justifies 'jump intensity' later)
#========================
tail_prob <- function(cty) {
  d <- panel[country == cty]
  thr <- quantile(d$r_embi, TAIL_Q, na.rm=TRUE)
  d[, bad := as.integer(r_embi <= thr)]
  
  if (cty %in% COPPER_SET) {
    crash <- "copper_crash"
  } else if (cty %in% OIL_SET) {
    crash <- "oil_crash"
  } else {
    crash <- "oil_crash"
  }
  
  p_crash <- d[get(crash)==1, mean(bad, na.rm=TRUE)]
  p_norm  <- d[get(crash)==0, mean(bad, na.rm=TRUE)]
  data.table(
    country=cty,
    crash_var=crash,
    p_bad_crash=p_crash,
    p_bad_norm=p_norm,
    ratio=p_crash/p_norm
  )
}

tail_res <- rbindlist(lapply(COUNTRIES, tail_prob))
tail_res[, group := fifelse(country %in% OIL_SET, "Oil exporter",
                            fifelse(country %in% COPPER_SET, "Copper exporter", "Control"))]
fwrite(tail_res, file.path(OUT_DIR, "tail_prob_results.csv"))

#========================
# 10) EXPOSURE -> BETA CHECK (cross-sectional sanity)
#========================
exp_avg <- exp_use[country %in% COUNTRIES,
                   .(oil_exp_avg = mean(oil_exp, na.rm=TRUE),
                     copper_exp_avg = mean(copper_exp, na.rm=TRUE)),
                   by=country]

res_xs <- merge(res_main, exp_avg, by="country", all.x=TRUE)

oil_xs <- res_xs[driver=="r_oil"][, .(country, group, beta, beta_crash, delta_crash_minus_norm, p, p_crash, p_delta, oil_exp_avg)]
cu_xs  <- res_xs[driver=="r_copper"][, .(country, group, beta, beta_crash, delta_crash_minus_norm, p, p_crash, p_delta, copper_exp_avg)]

fwrite(oil_xs, file.path(OUT_DIR, "exposure_beta_oil.csv"))
fwrite(cu_xs,  file.path(OUT_DIR, "exposure_beta_copper.csv"))

# plots
p1 <- ggplot(oil_xs, aes(oil_exp_avg, beta, label=country)) +
  geom_point(size=3) + geom_smooth(method="lm", se=FALSE) +
  geom_text(vjust=-0.6) + theme_minimal() +
  labs(title="Oil exposure vs oil beta (per-country)", x="Avg oil exposure share", y="Beta on r_oil")
ggsave(file.path(OUT_DIR, "xs_oil_exposure_beta.png"), p1, width=8, height=5, dpi=150)

p2 <- ggplot(cu_xs, aes(copper_exp_avg, beta, label=country)) +
  geom_point(size=3) + geom_smooth(method="lm", se=FALSE) +
  geom_text(vjust=-0.6) + theme_minimal() +
  labs(title="Copper exposure vs copper beta (per-country)", x="Avg copper exposure share", y="Beta on r_copper")
ggsave(file.path(OUT_DIR, "xs_copper_exposure_beta.png"), p2, width=8, height=5, dpi=150)

#========================
# 11) EVENT STUDY PLOTS (exporters vs controls)
#========================
event_study <- function(dt_in, crash_var, grp_var, window=-2:2) {
  dtx <- copy(dt_in)
  setorder(dtx, country, date)
  dtx[, t := .I, by=country]
  ev <- dtx[get(crash_var)==1, .(country, grp=get(grp_var), t0=t)]
  if (nrow(ev)==0) return(NULL)
  es <- ev[, .(rel=window), by=.(country, grp, t0)]
  es[, t := t0 + rel]
  es <- merge(es, dtx[, .(country, t, r_embi)], by=c("country","t"), all.x=TRUE)
  es[, .(avg_r_embi = mean(r_embi, na.rm=TRUE)), by=.(grp, rel)][]
}

# Oil event study (oil exporters vs controls)
dt_oil <- panel[country %in% c(OIL_SET, CTRL_SET)]
dt_oil[, grp := ifelse(country %in% OIL_SET, "Oil exporters", "Controls")]
es_oil <- event_study(dt_oil, "oil_crash", "grp", EVENT_WINDOW)
if (!is.null(es_oil)) {
  p_oil <- ggplot(es_oil, aes(rel, avg_r_embi, color=grp)) +
    geom_line() + geom_point() + theme_minimal() +
    labs(title=paste0("Event study: EMBI returns around OIL crash weeks (q=",SHOCK_Q,")"),
         x="Weeks relative to crash (t=0)", y="Avg weekly EMBI return", color="")
  ggsave(file.path(OUT_DIR, "eventstudy_oil.png"), p_oil, width=9, height=5, dpi=150)
}

# Copper event study (copper exporters vs controls)
dt_cu <- panel[country %in% c(COPPER_SET, CTRL_SET)]
dt_cu[, grp := ifelse(country %in% COPPER_SET, "Copper exporters", "Controls")]
es_cu <- event_study(dt_cu, "copper_crash", "grp", EVENT_WINDOW)
if (!is.null(es_cu)) {
  p_cu <- ggplot(es_cu, aes(rel, avg_r_embi, color=grp)) +
    geom_line() + geom_point() + theme_minimal() +
    labs(title=paste0("Event study: EMBI returns around COPPER crash weeks (q=",SHOCK_Q,")"),
         x="Weeks relative to crash (t=0)", y="Avg weekly EMBI return", color="")
  ggsave(file.path(OUT_DIR, "eventstudy_copper.png"), p_cu, width=9, height=5, dpi=150)
}

#========================
# 12) POOLED INTERACTION MODELS (links exposure to sensitivity)
#========================
# NOTE: without date FE so r_oil/r_copper are identifiable (still clustered by date)
panel[, is_oil_exporter := as.integer(country %in% OIL_SET)]
panel[, is_cu_exporter  := as.integer(country %in% COPPER_SET)]

m_pool_oil <- feols(
  r_embi ~ r_oil*oil_exp + r_spx + r_vix + d_us10y + r_fx | country,
  data = panel[country %in% c(OIL_SET, CTRL_SET)],
  vcov = ~date
)

m_pool_cu <- feols(
  r_embi ~ r_copper*copper_exp + r_spx + r_vix + d_us10y + r_fx | country,
  data = panel[country %in% c(COPPER_SET, CTRL_SET)],
  vcov = ~date
)

# Crash-only slope pooled (state dependence)
panel[, oil_ret_crash_exp := r_oil*oil_exp*oil_crash]
panel[, oil_ret_norm_exp  := r_oil*oil_exp*(1-oil_crash)]
panel[, cu_ret_crash_exp  := r_copper*copper_exp*copper_crash]
panel[, cu_ret_norm_exp   := r_copper*copper_exp*(1-copper_crash)]

m_pool_oil_piece <- feols(
  r_embi ~ oil_ret_norm_exp + oil_ret_crash_exp + r_spx + r_vix + d_us10y + r_fx | country,
  data = panel[country %in% c(OIL_SET, CTRL_SET)],
  vcov = ~date
)

m_pool_cu_piece <- feols(
  r_embi ~ cu_ret_norm_exp + cu_ret_crash_exp + r_spx + r_vix + d_us10y + r_fx | country,
  data = panel[country %in% c(COPPER_SET, CTRL_SET)],
  vcov = ~date
)

sink(file.path(OUT_DIR, "pooled_models.txt"))
cat("=== POOLED MODELS (country FE, cluster by date) ===\n\n")
cat("\n--- Oil: r_embi ~ r_oil*oil_exp + controls ---\n")
print(summary(m_pool_oil))
cat("\n--- Copper: r_embi ~ r_copper*copper_exp + controls ---\n")
print(summary(m_pool_cu))

cat("\n\n=== POOLED PIECEWISE (state dependence) ===\n\n")
cat("\n--- Oil: norm vs crash exposure-weighted slopes ---\n")
print(summary(m_pool_oil_piece))
cat("\n--- Copper: norm vs crash exposure-weighted slopes ---\n")
print(summary(m_pool_cu_piece))
sink()

#========================
# 13) ROBUSTNESS: EXCLUDE 2025 (exposure imputed)
#========================
panel_pre2025 <- panel[year <= 2024]

m_pre_oil <- feols(
  r_embi ~ r_oil*oil_exp + r_spx + r_vix + d_us10y + r_fx | country,
  data = panel_pre2025[country %in% c(OIL_SET, CTRL_SET)],
  vcov = ~date
)

m_pre_cu <- feols(
  r_embi ~ r_copper*copper_exp + r_spx + r_vix + d_us10y + r_fx | country,
  data = panel_pre2025[country %in% c(COPPER_SET, CTRL_SET)],
  vcov = ~date
)

sink(file.path(OUT_DIR, "robustness_pre2025.txt"))
cat("=== ROBUSTNESS: pre-2025 only ===\n\n")
cat("\n--- Oil pooled interaction pre-2025 ---\n")
print(summary(m_pre_oil))
cat("\n--- Copper pooled interaction pre-2025 ---\n")
print(summary(m_pre_cu))
sink()

cat("\nDONE. Files in: ", normalizePath(OUT_DIR), "\n")
cat("Key outputs:\n",
    "- per_country_results.csv\n",
    "- tail_prob_results.csv\n",
    "- exposure_beta_oil.csv / exposure_beta_copper.csv\n",
    "- pooled_models.txt\n",
    "- per_country_summaries.txt\n",
    "- eventstudy_oil.png / eventstudy_copper.png\n",
    "- robustness_pre2025.txt\n")

