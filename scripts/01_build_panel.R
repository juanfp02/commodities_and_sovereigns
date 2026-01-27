############################################################
# 01_build_panel.R
# Build WEEKLY PANEL (LEVELS ONLY) from pre-cleaned RDS inputs
# - No transformations (no log-diffs / returns)
# - Long format: one row per (country, date)
# - Merges exposures (annual) and fills within country
############################################################

suppressPackageStartupMessages({
  library(data.table)
})

#========================
# 0) PATHS
#========================
RDS_DIR_MASTER <- "data/processed/masterfile_rds"
RDS_EXPOSURE   <- "data/processed/exposure/exposure_bcom_sectors_2000_2024.rds"

OUT_DIR <- "data/processed"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

OUT_RDS <- file.path(OUT_DIR, "panel_weekly_levels.rds")

#========================
# 1) LOAD INPUTS (RDS ONLY)
#========================

EMBIs_dt           <- readRDS(file.path(RDS_DIR_MASTER, "EMBIs.rds"))
FXRates_dt         <- readRDS(file.path(RDS_DIR_MASTER, "FXRates.rds"))
Global_Controls_dt <- readRDS(file.path(RDS_DIR_MASTER, "Global_Controls.rds"))
Commodities_dt     <- readRDS(file.path(RDS_DIR_MASTER, "Commodities.rds"))
country_map   <- readRDS(file.path(RDS_DIR_MASTER, "country_map.rds"))
exposure_clean <- readRDS(RDS_EXPOSURE)

setDT(EMBIs_dt); setDT(FXRates_dt); setDT(Global_Controls_dt); setDT(Commodities_dt)
setDT(country_map); setDT(exposure_clean)

#========================
# 2) USER SETTINGS (EDIT HERE)
#========================
COUNTRIES <- country_map$country
#COUNTRIES <- c("Angola","Iraq","Saudi Arabia","UAE","Qatar","Kazakhstan","Colombia","Ecuador", "Poland","Morocco","Egypt","Mexico","Brazil","Chile","Peru","Zambia","Senegal")

# Commodity price columns to include (levels)
commodity_cols <- c(
  "CL1_Comdty",      # Brent crude (example)
  "CO1_Comdty",
  "HG1_Comdty",      # copper (example)
  "BCOMEN_Index",
  "BCOMAG_Index",     # Agricultural prices
  "BCOMPR_Index",
  "BCOMIN_Index"
)

# Control columns to include (levels)
control_cols <- c(
  "SPX_Index",
  "VIX_Index",
  "USGG10YR_Index",
  "DXY_Curncy"      # if present in your Global_Controls sheet
)

# Exposure columns to merge from exposure_clean (annual)
# Use EXACT column names from exposure_clean_annual.rds
# (common: share_crude_oil, expgdp_crude_oil, share_copper_total, expgdp_copper_total)
exposure_cols <- c(
  "share_energy",
  "expgdp_energy",
  "share_industrial_metals",
  "expgdp_industrial_metals",
  "share_precious",
  "expgdp_precious",
  "share_agriculture",
  "expgdp_agriculture"
  
)
exposure_cols <- intersect(exposure_cols, names(exposure_clean))

#========================
# 3) BUILD MAPS (from country_map)
#========================
cm <- copy(country_map)[country %in% COUNTRIES]
stopifnot(nrow(cm) > 0)

embi_map <- cm[, .(country, iso3, embi_col)]
fx_map   <- cm[, .(country, iso3, fx_col)]

if (!"USD_DUMMY" %in% names(FXRates_dt)) {
  FXRates_dt[, USD_DUMMY := 1]
}
#========================
# 4) WIDE -> LONG (LEVELS)
#========================
# EMBI long
embi_long <- melt(
  EMBIs_dt[, c("date", embi_map$embi_col), with = FALSE],
  id.vars = "date",
  measure.vars = embi_map$embi_col,
  variable.name = "embi_col",
  value.name = "embi_level"
)
embi_long <- merge(embi_long, embi_map, by = "embi_col", all.x = TRUE)

# FX long
# --- make sure USD_DUMMY exists (always)
if (!"USD_DUMMY" %in% names(FXRates_dt)) FXRates_dt[, USD_DUMMY := 1]

# fx_long without merge() mapping (robust even when fx_col repeats)
fx_long <- rbindlist(lapply(seq_len(nrow(fx_map)), function(i){
  data.table(
    date     = FXRates_dt$date,
    country  = fx_map$country[i],
    iso3     = fx_map$iso3[i],
    fx_level = FXRates_dt[[ fx_map$fx_col[i] ]]
  )
}), use.names = TRUE)

# Merge EMBI + FX
panel <- merge(
  embi_long[, .(date, country, iso3, embi_level)],
  fx_long[,   .(date, country, iso3, fx_level)],
  by = c("date","country","iso3"),
  all.x = TRUE
)

# Merge global controls (levels)
panel <- merge(
  panel,
  Global_Controls_dt[, c("date", control_cols), with = FALSE],
  by = "date",
  all.x = TRUE
)

# Merge commodities (levels)
panel <- merge(
  panel,
  Commodities_dt[, c("date", commodity_cols), with = FALSE],
  by = "date",
  all.x = TRUE
)

#========================
# 5) YEAR + MERGE EXPOSURES (annual)
#========================
panel[, year := as.integer(format(date, "%Y"))]

# Keep only what we need from exposure_clean
exp_keep <- c("iso3","year", exposure_cols)
exp_keep <- intersect(exp_keep, names(exposure_clean))
exp_dt <- exposure_clean[, ..exp_keep]

setkey(exp_dt, iso3, year)
setkey(panel, iso3, year)

panel <- exp_dt[panel]  # left join: keeps all panel rows

# Fill exposure within country (covers missing years like 2025)
setorder(panel, iso3, date)
if (length(exposure_cols)) {
  for (cc in exposure_cols) {
    panel[, (cc) := nafill(get(cc), type = "locf"), by = iso3]
    panel[, (cc) := nafill(get(cc), type = "nocb"), by = iso3]
  }
}

#========================
# 6) FINAL CLEAN + SAVE
#========================
# Keep only selected countries (safety)
panel <- panel[country %in% COUNTRIES]

# Sort
setorder(panel, country, date)

saveRDS(panel, OUT_RDS)
