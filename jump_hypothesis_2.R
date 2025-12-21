############################################################
# claim2_jump_justification_pipeline.R
# Purpose: Claim 2 = Commodity crashes increase tail/jump sovereign moves
# DV: EMBI index returns (log returns)
# Models:
#   (A) Jump probability: feglm(jump_neg ~ crash*exposure + controls | country+date)
#   (B) Jump severity (jump weeks only): feols(jump_mag ~ crash*exposure + controls | country+date)
############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(fixest)
  library(ggplot2)
})

#========================
# 0) Load weekly Bloomberg data (final_dt) + exposures
#========================
source("data_cleaning.R")              # must create final_dt
dt <- copy(final_dt); setDT(dt)

# annual exposure 2013-2024
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
# 1) User inputs
#========================
OUT_DIR <- "output_claim2"
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR)

# Choose your sets (EDIT)
OIL_SET    <- c("Angola","Nigeria","Ecuador","Kazakhstan","Colombia")
COPPER_SET <- c("Chile","Peru","Zambia")
CTRL_SET   <- c("Brazil","Egypt","Morocco")

KEEP_COUNTRIES <- unique(c(OIL_SET, COPPER_SET, CTRL_SET))

# Quantiles
SHOCK_Q <- 0.10  # crash = bottom 10% commodity returns
JUMP_Q  <- 0.05  # jump_neg = bottom 5% EMBI returns (within each country)

EVENT_WINDOW <- -4:4

#========================
# 2) Clean junk columns that are literally "Requesting Data" etc
#========================
bad_col_pattern <- "Requesting Data|#N/A|N/A"
bad_cols <- names(dt)[grepl(bad_col_pattern, names(dt), ignore.case = TRUE)]
if (length(bad_cols) > 0) {
  message("Dropping junk columns: ", paste(bad_cols, collapse = ", "))
  dt[, (bad_cols) := NULL]
}

#========================
# 3) Map EMBI index levels (EDIT to match your columns)
#========================
embi_map <- data.table(
  country = c("Angola","Nigeria","Ecuador","Chile","Peru","Morocco","Egypt","Brazil","Colombia","Mexico","Zambia"),
  embi_col = c("JPGCAO_Index","JPGCNG_Index","JPGCEC_Index","JPGCCH_Index","JPGCPR_Index",
               "JPGCMOR_Index","JPGCEG_Index","JPGCBZ_Index","JPGCCO_Index","JPGCMX_Index","JPGCZM_Index")
)
embi_map <- embi_map[country %in% KEEP_COUNTRIES]

missing_embi_cols <- setdiff(embi_map$embi_col, names(dt))
if (length(missing_embi_cols) > 0) {
  stop("Missing EMBI index columns in final_dt:\n", paste(missing_embi_cols, collapse = ", "))
}

#========================
# 4) Global series + FX controls (optional)
#========================
global_cols <- c("date", "USGG10YR_Index", "SPX_Index", "VIXCLS", "CO1_Comdty", "LP1_Comdty")
miss_global <- setdiff(global_cols, names(dt))
if (length(miss_global) > 0) stop("Missing global columns:\n", paste(miss_global, collapse=", "))

USE_FX <- TRUE
dt[, ECU_FX_DUMMY := 1]  # Ecuador dollarized

fx_map <- data.table(
  country = c("Angola","Nigeria","Ecuador","Chile","Peru","Morocco","Egypt","Brazil","Colombia","Zambia","Kazakhstan"),
  fx_col  = c("USDAOA_Curncy","USDNGN_Curncy","ECU_FX_DUMMY","USDCLP_Curncy","USDPEN_Curncy",
              "USDMAD_Curncy","USDEGP_Curncy","USDBRL_Curncy","USDCOP_Curncy","USDZMW_Curncy","USDKZT_Curncy")
)
fx_map <- fx_map[country %in% KEEP_COUNTRIES]

if (USE_FX) {
  missing_fx_cols <- setdiff(fx_map$fx_col, names(dt))
  if (length(missing_fx_cols) > 0) stop("Missing FX columns:\n", paste(missing_fx_cols, collapse=", "))
}

#========================
# 5) Build long panel: EMBI index levels -> returns
#========================
dt[, date := as.Date(date)]
setorder(dt, date)

need_cols <- unique(c(global_cols, embi_map$embi_col, if (USE_FX) fx_map$fx_col))
panel_w <- dt[, ..need_cols]

embi_long <- melt(
  panel_w,
  id.vars = global_cols,
  measure.vars = embi_map$embi_col,
  variable.name = "embi_col",
  value.name = "embi_level"
)
embi_long <- merge(embi_long, embi_map, by = "embi_col", all.x = TRUE)

if (USE_FX) {
  fx_long <- melt(
    panel_w,
    id.vars = "date",
    measure.vars = fx_map$fx_col,
    variable.name = "fx_col",
    value.name = "fx_usdxxx"
  )
  fx_long <- merge(fx_long, fx_map, by = "fx_col", all.x = TRUE)
  
  panel <- merge(
    embi_long[, .(date, country, embi_level, USGG10YR_Index, SPX_Index, VIXCLS, CO1_Comdty, LP1_Comdty)],
    fx_long[, .(date, country, fx_usdxxx)],
    by = c("date","country"),
    all.x = TRUE
  )
} else {
  panel <- embi_long[, .(date, country, embi_level, USGG10YR_Index, SPX_Index, VIXCLS, CO1_Comdty, LP1_Comdty)]
  panel[, fx_usdxxx := NA_real_]
}
panel <- panel[country %in% KEEP_COUNTRIES]
setorder(panel, country, date)
panel[, embi_level := as.numeric(embi_level)]

# returns + controls
panel[, r_embi := log(embi_level) - log(shift(embi_level, 1)), by = country]
panel[, r_spx  := log(SPX_Index)   - log(shift(SPX_Index, 1)), by = country]
panel[, r_vix  := log(VIXCLS)      - log(shift(VIXCLS, 1)), by = country]
panel[, d_us10y := USGG10YR_Index - shift(USGG10YR_Index, 1), by = country]
panel[, r_oil   := log(CO1_Comdty) - log(shift(CO1_Comdty, 1)), by = country]
panel[, r_copper:= log(LP1_Comdty) - log(shift(LP1_Comdty, 1)), by = country]

if (USE_FX) {
  panel[, r_fx := log(fx_usdxxx) - log(shift(fx_usdxxx, 1)), by = country]
} else {
  panel[, r_fx := 0]
}

panel <- panel[!is.na(r_embi) & !is.na(r_oil) & !is.na(r_copper)]
panel[, year := as.integer(format(date, "%Y"))]

# merge exposures (2013-2024), carry 2024 into 2025 if needed
panel <- merge(panel, exp_use, by = c("country","year"), all.x = TRUE)
setorder(panel, country, date)
panel[, oil_exp := nafill(oil_exp, type="locf"), by=country]
panel[, oil_exp := nafill(oil_exp, type="nocb"), by=country]
panel[, copper_exp := nafill(copper_exp, type="locf"), by=country]
panel[, copper_exp := nafill(copper_exp, type="nocb"), by=country]

#========================
# 6) Define crash weeks + jump weeks
#========================
q_oil <- quantile(panel$r_oil, SHOCK_Q, na.rm=TRUE)
q_cu  <- quantile(panel$r_copper, SHOCK_Q, na.rm=TRUE)

panel[, oil_crash    := as.integer(r_oil    <= q_oil)]
panel[, copper_crash := as.integer(r_copper <= q_cu)]

# exposure-weighted crash variables (THIS is what survives date FE)
panel[, oil_crash_exp    := oil_exp    * oil_crash]
panel[, copper_crash_exp := copper_exp * copper_crash]

# "jump week" = very negative sovereign return, within-country
panel[, jump_neg := as.integer(r_embi <= quantile(r_embi, JUMP_Q, na.rm=TRUE)), by = country]

# severity measure on jump weeks (positive = worse)
panel[, jump_mag := -r_embi]
panel[, jump_mag := fifelse(jump_neg == 1, jump_mag, NA_real_)]

#========================
# 7) Core Claim 2 models (date FE design)
#========================
run_claim2_models <- function(dtx, crash_exp_var, label){
  # jump probability
  m_prob <- feglm(
    as.formula(paste0("jump_neg ~ ", crash_exp_var, " + r_fx | country + date")),
    data = dtx,
    family = "binomial",
    vcov = ~country
  )
  
  # jump severity (jump weeks only)
  d_jump <- dtx[jump_neg == 1]
  m_sev <- feols(
    as.formula(paste0("jump_mag ~ ", crash_exp_var, " + r_fx | country + date")),
    data = d_jump,
    vcov = ~country
  )
  
  # save
  sink(file.path(OUT_DIR, paste0("claim2_", label, "_models.txt")))
  cat("=== CLAIM 2 MODELS: ", label, " ===\n\n")
  cat("\n--- Jump probability (country+date FE) ---\n")
  print(summary(m_prob))
  cat("\n--- Jump severity (jump weeks only; country+date FE) ---\n")
  print(summary(m_sev))
  sink()
  
  list(prob = m_prob, sev = m_sev)
}

dt_oil <- panel[country %in% c(OIL_SET, CTRL_SET)]
dt_cu  <- panel[country %in% c(COPPER_SET, CTRL_SET)]

mods_oil <- run_claim2_models(dt_oil, "oil_crash_exp", "OIL_exporters_vs_controls")
mods_cu  <- run_claim2_models(dt_cu,  "copper_crash_exp", "COPPER_exporters_vs_controls")

#========================
# 8) Diagnostics: tail probability ratio by country
#========================
tail_diag <- function(dtx, crash_var){
  dtx[, .(
    n = .N,
    n_crash = sum(get(crash_var)==1, na.rm=TRUE),
    p_jump_crash = mean(jump_neg[get(crash_var)==1], na.rm=TRUE),
    p_jump_norm  = mean(jump_neg[get(crash_var)==0], na.rm=TRUE),
    ratio = mean(jump_neg[get(crash_var)==1], na.rm=TRUE) /
      mean(jump_neg[get(crash_var)==0], na.rm=TRUE)
  ), by = country][order(-ratio)]
}

diag_oil <- tail_diag(dt_oil, "oil_crash")
diag_cu  <- tail_diag(dt_cu,  "copper_crash")
fwrite(diag_oil, file.path(OUT_DIR, "tail_diag_oil.csv"))
fwrite(diag_cu,  file.path(OUT_DIR, "tail_diag_copper.csv"))

#========================
# 9) Event study around crash weeks (exporters vs controls)
#========================
event_study <- function(dtx, crash_var, group_var, window=-4:4){
  dd <- copy(dtx)
  setorder(dd, country, date)
  dd[, t := .I, by=country]
  
  ev <- dd[get(crash_var)==1, .(country, grp=get(group_var), t0=t)]
  if (nrow(ev)==0) stop("No events for ", crash_var)
  
  es <- ev[, .(rel=window), by=.(country, grp, t0)]
  es[, t := t0 + rel]
  es <- merge(es, dd[, .(country, t, r_embi)], by=c("country","t"), all.x=TRUE)
  es[, .(avg_r_embi = mean(r_embi, na.rm=TRUE)), by=.(grp, rel)][]
}

# group labels
dt_oil[, grp := ifelse(country %in% OIL_SET, "Oil exporters", "Controls")]
dt_cu[,  grp := ifelse(country %in% COPPER_SET, "Copper exporters", "Controls")]

es_oil <- event_study(dt_oil, "oil_crash", "grp", EVENT_WINDOW)
p_oil <- ggplot(es_oil, aes(rel, avg_r_embi, color=grp)) +
  geom_line() + geom_point() +
  theme_minimal() +
  labs(title=paste0("Event study: oil crash weeks (SHOCK_Q=",SHOCK_Q,")"),
       x="Weeks relative to crash (t=0)", y="Avg weekly EMBI return", color="")
ggsave(file.path(OUT_DIR, "eventstudy_oil.png"), p_oil, width=9, height=5, dpi=150)

es_cu <- event_study(dt_cu, "copper_crash", "grp", EVENT_WINDOW)
p_cu <- ggplot(es_cu, aes(rel, avg_r_embi, color=grp)) +
  geom_line() + geom_point() +
  theme_minimal() +
  labs(title=paste0("Event study: copper crash weeks (SHOCK_Q=",SHOCK_Q,")"),
       x="Weeks relative to crash (t=0)", y="Avg weekly EMBI return", color="")
ggsave(file.path(OUT_DIR, "eventstudy_copper.png"), p_cu, width=9, height=5, dpi=150)

#========================
# 10) Robustness: exclude 2025 (exposure ends in 2024)
#========================
panel_pre2025 <- panel[year <= 2024]
dt_oil_pre <- panel_pre2025[country %in% c(OIL_SET, CTRL_SET)]
dt_cu_pre  <- panel_pre2025[country %in% c(COPPER_SET, CTRL_SET)]

mods_oil_pre <- run_claim2_models(dt_oil_pre, "oil_crash_exp", "OIL_pre2025")
mods_cu_pre  <- run_claim2_models(dt_cu_pre,  "copper_crash_exp", "COPPER_pre2025")

#========================
# 11) Sensitivity loop over (SHOCK_Q, JUMP_Q)
#========================
run_grid <- function(SHOCK_GRID=c(0.05,0.10), JUMP_GRID=c(0.10,0.05,0.02)){
  out <- list()
  for (sq in SHOCK_GRID) for (jq in JUMP_GRID){
    tmp <- copy(panel)
    q_o <- quantile(tmp$r_oil, sq, na.rm=TRUE)
    q_c <- quantile(tmp$r_copper, sq, na.rm=TRUE)
    tmp[, oil_crash := as.integer(r_oil <= q_o)]
    tmp[, copper_crash := as.integer(r_copper <= q_c)]
    tmp[, oil_crash_exp := oil_exp * oil_crash]
    tmp[, copper_crash_exp := copper_exp * copper_crash]
    tmp[, jump_neg := as.integer(r_embi <= quantile(r_embi, jq, na.rm=TRUE)), by=country]
    tmp[, jump_mag := fifelse(jump_neg==1, -r_embi, NA_real_)]
    
    d_o <- tmp[country %in% c(OIL_SET, CTRL_SET)]
    d_c <- tmp[country %in% c(COPPER_SET, CTRL_SET)]
    
    m_o <- feglm(jump_neg ~ oil_crash_exp + r_fx | country + date,
                 data=d_o, family="binomial", vcov=~country)
    m_c <- feglm(jump_neg ~ copper_crash_exp + r_fx | country + date,
                 data=d_c, family="binomial", vcov=~country)
    
    out[[paste0("SH",sq,"_J",jq)]] <- data.table(
      shock_q = sq, jump_q = jq,
      oil_beta = coef(m_o)["oil_crash_exp"],
      oil_p    = summary(m_o)$coeftable["oil_crash_exp","Pr(>|z|)"],
      cu_beta  = coef(m_c)["copper_crash_exp"],
      cu_p     = summary(m_c)$coeftable["copper_crash_exp","Pr(>|z|)"]
    )
  }
  rbindlist(out)
}
grid_res <- run_grid()
fwrite(grid_res, file.path(OUT_DIR, "sensitivity_grid.csv"))

cat("\nDONE. Outputs in:", normalizePath(OUT_DIR), "\n")

