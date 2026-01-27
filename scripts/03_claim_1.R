suppressPackageStartupMessages({
  library(data.table)
  library(fixest)
  library(lmtest)
  library(sandwich)
  library(ggplot2)
})

# -----------------------
# 0) USER INPUTS (EDIT ONLY THESE)
# -----------------------
COMMOD_COL   <- "BCOMEN_Index"   # e.g. "BCOMEN_Index", "BCOMIN_Index", "BCOMPR_Index", "BCOMAG_Index"
EXPOSURE_COL <- "share_energy"   # e.g. "share_energy", "share_indmetals", "share_precious", "share_agriculture"
NW_LAG <- 4
# HARD-CODE GROUPS (EDIT)
TEST_GROUP <- c("Qatar","Angola","Iraq","Saudi Arabia","Nigeria","Kazakhstan","Ecuador","Colombia","Bahrain","United Arab Emirates")
#TEST_GROUP <- c("Zambia","Chile","Kazakhstan","Bahrain", "Peru")
#TEST_GROUP <- c("United Arab Emirates","Peru","Senegal","Ivory Coast", "South Africa")

CONTROLS   <- c("Brazil","Mexico","Poland","Morocco","Turkey")

KEEP <- unique(c(TEST_GROUP, CONTROLS))

# -----------------------
# 1) Load + filter
# -----------------------
panel <- readRDS("data/processed/panel_weekly_levels.rds")
setDT(panel)

stopifnot(COMMOD_COL %in% names(panel))
stopifnot(EXPOSURE_COL %in% names(panel))

panel <- panel[country %in% KEEP]
panel[, grp := fifelse(country %in% TEST_GROUP, "Test",
                       fifelse(country %in% CONTROLS, "Control", NA_character_))]
panel <- panel[!is.na(grp)]
panel[, treat := as.integer(country %in% TEST_GROUP)]   # <-- THIS IS THE DUMMY
setorder(panel, country, date)

cat("Kept countries: ", length(unique(panel$country)), "\n",
    "Test group: ", length(unique(panel[treat==1, country])), "\n",
    "Controls: ", length(unique(panel[treat==0, country])), "\n", sep="")

# -----------------------
# 2) Transforms (weekly)
# -----------------------
panel[, d_log_embi := log(embi_level) - log(shift(embi_level, 1)), by = country]
panel[, d_log_fx   := log(fx_level)   - log(shift(fx_level, 1)),   by = country]

# Global series ONCE per date
g <- unique(panel[, .(
  date,
  commod_level = get(COMMOD_COL),
  spx_level    = SPX_Index,
  u10y_level   = USGG10YR_Index,
)])
setorder(g, date)
g[, d_log_commod := log(commod_level) - log(shift(commod_level, 1))]
g[, d_log_spx    := log(spx_level)    - log(shift(spx_level, 1))]
g[, d_log_vix   := log(spx_level)    - log(shift(spx_level, 1))]

g[, diff_us10y   := u10y_level        - shift(u10y_level, 1)]

panel <- merge(panel, g[, .(date, d_log_commod, d_log_spx, diff_us10y)], by="date", all.x=TRUE)

# stable exposure name (continuous)
panel[, exposure := get(EXPOSURE_COL)]

# sample
panel <- panel[
  is.finite(d_log_embi) &
    is.finite(d_log_commod) &
    is.finite(d_log_fx) &
    is.finite(d_log_spx) &
    is.finite(diff_us10y) &
    is.finite(exposure)
]

# -----------------------
# 3) Pooled regressions
# -----------------------

# A) baseline (no date FE) – continuous exposure version
m_A_cont <- feols(
  d_log_embi ~ d_log_commod + d_log_commod:exposure +
    d_log_fx + diff_us10y + d_log_spx | country,
  data = panel,
  vcov = ~country
)


# A2) preferred (WITH date FE) – DUMMY version (cleanest for your "test vs control" setup)
m_A_dummy <- feols(
  d_log_embi ~ d_log_commod + d_log_commod:treat +
    d_log_fx + diff_us10y + d_log_spx | country,
  data = panel,
  vcov = ~country
)


# B) preferred (WITH date FE) – continuous exposure version
m_B_cont <- feols(
  d_log_embi ~ d_log_commod:exposure + d_log_fx | country + date,
  data = panel,
  vcov = ~country
)

# B2) preferred (WITH date FE) – DUMMY version (cleanest for your "test vs control" setup)
m_B_dummy <- feols(
  d_log_embi ~ d_log_commod:treat + d_log_fx | country + date,
  data = panel,
  vcov = ~country
)

res <- etable(m_A_cont, m_A_dummy, m_B_cont, m_B_dummy, digits=4, fitstat=c("r2","n"))

# -----------------------
# 4) Per-country betas (NO interaction)
# -----------------------
run_one_country <- function(cty, lag_nw = NW_LAG) {
  dtc <- panel[country == cty]
  
  m <- lm(d_log_embi ~ d_log_commod + d_log_fx + diff_us10y + d_log_spx, data = dtc)
  V <- NeweyWest(m, lag = lag_nw, prewhite = FALSE)
  ct <- coeftest(m, vcov. = V)
  
  data.table(
    country = cty,
    grp = unique(dtc$grp),
    n = nrow(dtc),
    beta_commod = ct["d_log_commod", 1],
    se_commod   = ct["d_log_commod", 2],
    p_commod    = ct["d_log_commod", 4],
    exp_mean    = mean(dtc$exposure, na.rm = TRUE)
  )
}

res_country <- rbindlist(lapply(sort(unique(panel$country)), run_one_country), fill=TRUE)

ggplot(res_country, aes(x = exp_mean, y = beta_commod, shape = grp)) +
  geom_point() +
  geom_smooth(method="lm", se=TRUE) +
  labs(
    x = paste0("Average ", EXPOSURE_COL, " (in sample)"),
    y = paste0("Country beta on ", COMMOD_COL),
    title = "Per-country betas vs exposure (test vs controls)"
  ) +
  theme_minimal()

# -----------------------
# SETTINGS
# -----------------------
SHOCK_Q <- 0.05   # bottom 5% for "crash" variables
JUMP_Q  <- 0.05   # bottom 5% EMBI returns within each country
EXPOSURE_COL <- "share_energy"  # keep simple

# -----------------------
# ASSUMES you already have:
# panel with columns:
# country, date, d_log_embi, d_log_fx,
# d_log_oil (or d_log_commod), d_log_spx, d_vix
# and exposure column (e.g., share_crude_oil)
# -----------------------

setDT(panel)
setorder(panel, country, date)
panel[, exposure := get(EXPOSURE_COL)]

# -----------------------
# 1) Define "jump weeks" in sovereign credit (within-country tail)
# -----------------------
panel[, jump_neg := as.integer(d_log_embi <= quantile(d_log_embi, JUMP_Q, na.rm=TRUE)), by=country]

# -----------------------
# 2) Define shocks (placebos use SAME tail frequency idea)
# -----------------------

# Oil crash weeks: bottom 5% of oil returns (global)
q_commodity <- quantile(panel$d_log_commod, SHOCK_Q, na.rm=TRUE)
panel[, oil_crash := as.integer(d_log_commod <= q_commodity)]

# Placebo 1: SPX crash weeks (bottom 5%)
q_spx <- quantile(panel$d_log_spx, SHOCK_Q, na.rm=TRUE)
panel[, spx_crash := as.integer(d_log_spx <= q_spx)]

# Clean sample
dt <- panel[
  is.finite(jump_neg) &
    is.finite(oil_crash) &
    is.finite(spx_crash) &
    is.finite(d_log_fx) &
    is.finite(exposure)
]

# -----------------------
# 3) Regressions (date FE kills "everything goes down" common component)
# -----------------------

# Main Claim 2
m_oil <- feols(
  jump_neg ~ oil_crash:exposure + d_log_fx | country + date,
  data = dt,
  vcov = ~country
)

# Placebo shocks
m_spx <- feols(
  jump_neg ~ spx_crash:exposure + d_log_fx | country + date,
  data = dt,
  vcov = ~country
)


# Horse-race: does oil survive once you allow "risk-off" heterogeneity?
m_all <- feols(
  jump_neg ~ oil_crash:exposure + spx_crash:exposure + d_log_fx | country + date,
  data = dt,
  vcov = ~country
)

res_3 <- etable(m_oil, m_spx, m_all, digits=4, fitstat=c("n","r2"))
