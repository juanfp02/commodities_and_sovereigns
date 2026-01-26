suppressPackageStartupMessages({
  library(data.table)
  library(fixest)
})

############################################################
# 0) HELPERS
############################################################

# safe log return
logret <- function(x) c(NA_real_, diff(log(x)))

# safe diff
diff1 <- function(x) c(NA_real_, diff(x))

# Add log returns for given cols (global series: no by-group)
add_global_log_returns <- function(dt, level_cols, prefix = "r_") {
  for (cc in level_cols) {
    newc <- paste0(prefix, gsub("[^A-Za-z0-9]+", "_", cc))
    dt[, (newc) := logret(get(cc))]
  }
  dt
}

# Add log returns for given cols by country (panel series)
add_panel_log_returns <- function(dt, level_cols, by = "country", prefix = "r_") {
  setorder(dt, get(by), date)
  for (cc in level_cols) {
    newc <- paste0(prefix, gsub("[^A-Za-z0-9]+", "_", cc))
    dt[, (newc) := logret(get(cc)), by = ..by]
  }
  dt
}

# Add diffs for global rates (like US10Y)
add_global_diffs <- function(dt, level_cols, prefix = "d_") {
  for (cc in level_cols) {
    newc <- paste0(prefix, gsub("[^A-Za-z0-9]+", "_", cc))
    dt[, (newc) := diff1(get(cc))]
  }
  dt
}

############################################################
# 1) BUILD MASTER PANEL (LEVELS) — DO THIS ONCE
############################################################
# Inputs you already have in memory:
# EMBIs_dt, FXRates_dt, Global_Controls_dt, Commodities_dt
# Each must have: date + wide columns

build_master_panel <- function(EMBIs_dt, FXRates_dt, Global_Controls_dt, Commodities_dt,
                               embi_map, fx_map,
                               commodity_cols,
                               ensure_ecu_dummy = TRUE) {
  
  # --- enforce data.table + date
  for (nm in c("EMBIs_dt","FXRates_dt","Global_Controls_dt","Commodities_dt")) {
    DT <- get(nm)
    setDT(DT)
    if (!"date" %in% names(DT)) stop(nm, " must have a 'date' column.")
    DT[, date := as.Date(date)]
    setkey(DT, date)
    assign(nm, DT, inherits = TRUE)
  }
  
  # Ecuador USD dummy (optional)
  if (ensure_ecu_dummy && !("ECU_FX_DUMMY" %in% names(FXRates_dt))) {
    FXRates_dt[, ECU_FX_DUMMY := 1]
  }
  
  # --- check required columns exist
  stopifnot(all(embi_map$embi_col %in% names(EMBIs_dt)))
  stopifnot(all(fx_map$fx_col %in% names(FXRates_dt)))
  stopifnot(all(commodity_cols %in% names(Commodities_dt)))
  
  # --- EMBI long (levels)
  embi_long <- melt(
    EMBIs_dt[, c("date", embi_map$embi_col), with = FALSE],
    id.vars = "date",
    measure.vars = embi_map$embi_col,
    variable.name = "embi_col",
    value.name = "embi_level"
  )
  embi_long <- merge(embi_long, embi_map, by = "embi_col", all.x = TRUE)
  
  # --- FX long (levels)
  fx_long <- melt(
    FXRates_dt[, c("date", fx_map$fx_col), with = FALSE],
    id.vars = "date",
    measure.vars = fx_map$fx_col,
    variable.name = "fx_col",
    value.name = "fx_level"
  )
  fx_long <- merge(fx_long, fx_map, by = "fx_col", all.x = TRUE)
  
  # --- merge EMBI + FX
  panel <- merge(
    embi_long[, .(date, country, embi_level)],
    fx_long[, .(date, country, fx_level)],
    by = c("date","country"),
    all.x = TRUE
  )
  
  # --- merge global controls (wide by date)
  panel <- merge(panel, Global_Controls_dt, by = "date", all.x = TRUE)
  
  # --- merge commodities (wide by date, keep only needed)
  comm_keep <- c("date", commodity_cols)
  panel <- merge(panel, Commodities_dt[, ..comm_keep], by = "date", all.x = TRUE)
  
  setorder(panel, country, date)
  panel
}

############################################################
# 2) TRANSFORM PANEL (RETURNS) — DO THIS ONCE
############################################################
transform_master_panel <- function(panel,
                                   embi_level_col = "embi_level",
                                   fx_level_col   = "fx_level",
                                   # controls (levels)
                                   spx_col  = "SPX_Index",
                                   vix_col  = "VIXCLS",
                                   us10y_col = "USGG10YR_Index",
                                   # commodities already in panel; compute returns for these
                                   commodity_cols,
                                   ecu_country_name = "Ecuador") {
  
  setDT(panel)
  panel[, date := as.Date(date)]
  setorder(panel, country, date)
  
  # Panel returns
  panel[, r_embi := logret(get(embi_level_col)), by = country]
  panel[, r_fx   := logret(get(fx_level_col)),   by = country]
  
  # Ecuador FX return = 0 (if FX is dummy or USDized)
  panel[country == ecu_country_name, r_fx := 0]
  
  # Global controls: returns/diffs (same for all countries; compute once per date)
  # We'll compute at date-level and merge back (avoids slight differences)
  g <- unique(panel[, .(date,
                        spx = if (spx_col %in% names(panel)) get(spx_col) else NA_real_,
                        vix = if (vix_col %in% names(panel)) get(vix_col) else NA_real_,
                        us10y = if (us10y_col %in% names(panel)) get(us10y_col) else NA_real_)])
  setorder(g, date)
  
  if (spx_col %in% names(panel)) g[, r_spx := logret(spx)]
  if (vix_col %in% names(panel)) g[, r_vix := logret(vix)]
  if (us10y_col %in% names(panel)) g[, d_us10y := diff1(us10y)]
  
  panel <- merge(panel, g[, .(date, r_spx, r_vix, d_us10y)], by = "date", all.x = TRUE)
  
  # Commodity returns: create r_<colname>
  # e.g. CO1_Comdty -> r_CO1_Comdty (cleaned to r_CO1_Comdty)
  for (cc in commodity_cols) {
    newc <- paste0("r_", gsub("[^A-Za-z0-9]+", "_", cc))
    # compute once per date (same for all countries)
    tmp <- unique(panel[, .(date, x = get(cc))])
    setorder(tmp, date)
    tmp[, (newc) := logret(x)]
    panel <- merge(panel, tmp[, .(date, get(newc))], by = "date", all.x = TRUE)
  }
  
  # Drop first rows where returns NA
  panel <- panel[is.finite(r_embi)]
  
  panel
}

############################################################
# 3) REGRESSION RUNNER — VERY SIMPLE TO SWAP COUNTRIES/COMMODITIES
############################################################
run_claim1_reg <- function(panel,
                           countries,
                           x_return,          # e.g. "r_CO1_Comdty" cleaned name
                           add_controls = TRUE,
                           fe = c("country","date"),  # best default
                           cluster = c("country")) {
  
  dt <- panel[country %in% countries]
  
  # build RHS
  rhs <- c(x_return, "r_fx")
  if (add_controls) {
    rhs <- c(rhs, intersect(c("r_vix","d_us10y","r_spx"), names(dt)))
  }
  
  # formula
  fe_part <- paste(fe, collapse = " + ")
  fml <- as.formula(paste0("r_embi ~ ", paste(rhs, collapse = " + "), " | ", fe_part))
  
  # vcov
  vc <- if (length(cluster) == 1) {
    as.formula(paste0("~", cluster[1]))
  } else {
    as.formula(paste0("~", paste(cluster, collapse = " + ")))
  }
  
  feols(fml, data = dt, vcov = vc)
}

############################################################
# 4) EXAMPLE: OIL EXPORTERS (CLAIM 1) — EDIT ONLY THIS BLOCK
############################################################

# ---- 4A) Put your mappings here (edit column names as needed)

# EMBI columns -> country names
embi_map <- data.table(
  country = c("Angola","Nigeria","Ecuador","Kazakhstan"),
  embi_col = c("JPGCAO_Index","JPGCNG_Index","JPGCEC_Index","JPGCKZ_Index")
)

# FX columns -> country names
fx_map <- data.table(
  country = c("Angola","Nigeria","Ecuador","Kazakhstan"),
  fx_col  = c("USDAOA_Curncy","USDNGN_Curncy","ECU_FX_DUMMY","USDKZT_Curncy")
)

# Which commodities do you want available in the panel (you can add many here)
commodity_cols <- c("CO1_Comdty")  # oil, copper (example)

# ---- 4B) Build panel once
panel_levels <- build_master_panel(
  EMBIs_dt, FXRates_dt, Global_Controls_dt, Commodities_dt,
  embi_map = embi_map,
  fx_map   = fx_map,
  commodity_cols = commodity_cols
)

# ---- 4C) Transform once
panel <- transform_master_panel(
  panel_levels,
  commodity_cols = commodity_cols,
  ecu_country_name = "Ecuador"
)

# ---- 4D) Run oil regression (swap x_return later for copper etc.)
OIL_EXPORTERS <- c("Angola","Nigeria","Ecuador","Kazakhstan")
x_oil <- "r_CO1_Comdty"   # because CO1_Comdty -> r_CO1_Comdty after cleaning

m_oil <- run_claim1_reg(
  panel,
  countries = OIL_EXPORTERS,
  x_return  = x_oil,
  add_controls = TRUE,
  fe = c("country","date"),
  cluster = c("country")
)

print(summary(m_oil))

# Optional: easy loop over commodities later
# x_vars <- c("r_CO1_Comdty","r_LP1_Comdty")
# models <- lapply(x_vars, function(xx) run_claim1_reg(panel, OIL_EXPORTERS, xx))
# etable(models)