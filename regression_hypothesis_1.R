############################################################
# per_country_embi_simple.R
# One regression per country, compare betas across countries
############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(lmtest)
  library(sandwich)
  library(ggplot2)
})

#========================
# 0) Load data
#========================
source("data_cleaning.R")
dt <- copy(final_dt)
setDT(dt)


dt[, date := as.Date(date)]
setorder(dt, date)

#========================
# 1) Define EMBI index columns (LEVELS)
#========================
embi_map <- data.table(
  country = c("Angola","Colombia","Ecuador","Chile","Peru","Zambia","Morocco","Egypt","Brazil"),
  embi_col = c("JPGCAO_Index","JPGCCO_Index","JPGCEC_Index","JPGCCH_Index","JPGCPR_Index", "JPGCZM_Index",
               "JPGCMOR_Index","JPGCEG_Index","JPGCBZ_Index")
)

# FX columns (Ecuador dollarized)
dt[, ECU_FX_DUMMY := 1]
fx_map <- data.table(
  country = c("Angola","Colombia","Ecuador","Chile","Peru","Zambia", "Morocco","Egypt","Brazil"),
  fx_col  = c("USDAOA_Curncy","USDCOP_Curncy","ECU_FX_DUMMY","USDCLP_Curncy","USDPEN_Curncy","USDZMW_Curncy","USDMAD_Curncy","USDEGP_Curncy","USDBRL_Curncy")
)

# Choose which countries you want in this run
OIL_SET    <- c("Angola","Colombia","Ecuador")
COPPER_SET <- c("Chile","Peru","Zambia")
CTRL_SET   <- c("Brazil","Egypt","Morocco")
COUNTRIES  <- unique(c(OIL_SET, COPPER_SET, CTRL_SET))

embi_map <- embi_map[country %in% COUNTRIES]
fx_map   <- fx_map[country %in% COUNTRIES]

# Globals
global_cols <- c("date","USGG10YR_Index","SPX_Index","VIXCLS","CO1_Comdty","LP1_Comdty")

need <- unique(c(global_cols, embi_map$embi_col, fx_map$fx_col))
miss <- setdiff(need, names(dt))
if (length(miss) > 0) stop("Missing columns: ", paste(miss, collapse=", "))

#========================
# 2) Build long panel with returns
#========================
panel_w <- dt[, ..need]

embi_long <- melt(panel_w, id.vars = global_cols,
                  measure.vars = embi_map$embi_col,
                  variable.name="embi_col", value.name="embi_level")
embi_long <- merge(embi_long, embi_map, by="embi_col", all.x=TRUE)

fx_long <- melt(panel_w, id.vars="date",
                measure.vars = fx_map$fx_col,
                variable.name="fx_col", value.name="fx_usdxxx")
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

# EMBI weekly log return (DV)
panel[, r_embi := log(embi_level) - log(shift(embi_level, 1)), by=country]

# Commodity returns (same across countries but fine)
panel[, r_oil := log(CO1_Comdty) - log(shift(CO1_Comdty, 1)), by=country]
panel[, r_copper := log(LP1_Comdty) - log(shift(LP1_Comdty, 1)), by=country]

# Controls
panel[, r_spx := log(SPX_Index) - log(shift(SPX_Index, 1)), by=country]
panel[, r_vix := log(VIXCLS) - log(shift(VIXCLS, 1)), by=country]
panel[, d_us10y := USGG10YR_Index - shift(USGG10YR_Index, 1), by=country]
panel[, r_fx := log(fx_usdxxx) - log(shift(fx_usdxxx, 1)), by=country]

panel <- panel[!is.na(r_embi) & !is.na(r_oil) & !is.na(r_copper) &
                 !is.na(r_spx) & !is.na(r_vix) & !is.na(d_us10y)]

#========================
# 3) Define crash weeks (global thresholds)
#========================
SHOCK_Q <- 0.1
q_oil <- quantile(unique(panel[, .(date, r_oil)])$r_oil, SHOCK_Q, na.rm=TRUE)
q_cu  <- quantile(unique(panel[, .(date, r_copper)])$r_copper, SHOCK_Q, na.rm=TRUE)

panel[, oil_crash := as.integer(r_oil < q_oil)]
panel[, copper_crash := as.integer(r_copper < q_cu)]

#========================
# 4) Per-country regressions (Newey-West SEs)
#========================
nw_lag <- 4 # weekly data, 4 weeks is a common choice

run_one_country <- function(cty) {
  d <- panel[country == cty]
  
  # Choose the commodity driver depending on the country group
  if (cty %in% COPPER_SET) {
    x <- "r_copper"
    crash <- "copper_crash"
  } else if (cty %in% OIL_SET) {
    x <- "r_oil"
    crash <- "oil_crash"
  } else {
    # controls: run BOTH oil and copper separately later;
    # for now pick oil as default
    x <- "r_oil"
    crash <- "oil_crash"
  }
  
  # Baseline regression
  f1 <- as.formula(paste0("r_embi ~ ", x, " + r_spx + r_vix + d_us10y + r_fx"))
  m1 <- lm(f1, data=d)
  se1 <- NeweyWest(m1, lag=nw_lag, prewhite=FALSE)
  ct1 <- coeftest(m1, vcov.=se1)
  
  # Crash-slope regression: separate slope during crash weeks
  d[, x_norm := get(x) * (1 - get(crash))]
  d[, x_crash := get(x) * get(crash)]
  f2 <- as.formula("r_embi ~ x_norm + x_crash + r_spx + r_vix + d_us10y + r_fx")
  m2 <- lm(f2, data=d)
  se2 <- NeweyWest(m2, lag=nw_lag, prewhite=FALSE)
  ct2 <- coeftest(m2, vcov.=se2)
  
  # Extract betas
  out <- data.table(
    country = cty,
    driver = x,
    n = nrow(d),
    
    beta = ct1[x, 1],
    t = ct1[x, 3],
    p = ct1[x, 4],
    
    beta_norm = ct2["x_norm", 1],
    p_norm = ct2["x_norm", 4],
    beta_crash = ct2["x_crash", 1],
    p_crash = ct2["x_crash", 4]
  )
  
  out
}

res <- rbindlist(lapply(COUNTRIES, run_one_country))

# Label group
res[, group := fifelse(country %in% OIL_SET, "Oil exporter",
                       fifelse(country %in% COPPER_SET, "Copper exporter", "Control"))]

# Save + view
fwrite(res, file="output/per_country_results.csv")
print(res[order(group, p)])

#========================
# 5) Controls: run BOTH oil and copper (so controls are not “oil by default”)
#========================
run_control_both <- function(cty, x, crash) {
  d <- panel[country == cty]
  f1 <- as.formula(paste0("r_embi ~ ", x, " + r_spx + r_vix + d_us10y + r_fx"))
  m1 <- lm(f1, data=d)
  ct1 <- coeftest(m1, vcov.=NeweyWest(m1, lag=nw_lag, prewhite=FALSE))
  
  d[, x_norm := get(x) * (1 - get(crash))]
  d[, x_crash := get(x) * get(crash)]
  m2 <- lm(r_embi ~ x_norm + x_crash + r_spx + r_vix + d_us10y + r_fx, data=d)
  ct2 <- coeftest(m2, vcov.=NeweyWest(m2, lag=nw_lag, prewhite=FALSE))
  
  data.table(country=cty, driver=x,
             beta=ct1[x,1], p=ct1[x,4],
             beta_crash=ct2["x_crash",1], p_crash=ct2["x_crash",4])
}

ctrl_oil <- rbindlist(lapply(CTRL_SET, run_control_both, x="r_oil", crash="oil_crash"))
ctrl_cu  <- rbindlist(lapply(CTRL_SET, run_control_both, x="r_copper", crash="copper_crash"))

fwrite(ctrl_oil, "output/controls_oil_results.csv")
fwrite(ctrl_cu,  "output/controls_copper_results.csv")

#========================
# 6) Quick comparison plots (betas by group)
#========================
ggplot(res, aes(x=country, y=beta, color=group)) +
  geom_point(size=3) +
  coord_flip() +
  labs(title="Per-country commodity beta (baseline regression)",
       y="beta on commodity return", x="") +
  theme_minimal()

ggsave("output/per_country_beta_plot.png", width=8, height=5, dpi=150)

ggplot(res, aes(x=country, y=beta_crash, color=group)) +
  geom_point(size=3) +
  coord_flip() +
  labs(title="Per-country commodity beta DURING CRASH weeks",
       y="beta_crash", x="") +
  theme_minimal()

ggsave("output/per_country_beta_crash_plot.png", width=8, height=5, dpi=150)

cat("\nDONE. Outputs in output/:\n",
    "- per_country_results.csv\n",
    "- controls_oil_results.csv, controls_copper_results.csv\n",
    "- per_country_beta_plot.png\n",
    "- per_country_beta_crash_plot.png\n")

