############################################################
# 01_ingest_clean.R
# FULL PIPELINE:
#  - Reads master Excel (XLSX) -> cleans -> saves sheet RDS
#  - Loads cleaned sheets into: EMBIs_dt, FXRates_dt, Commodities_dt, Global_Controls_dt
#  - Defines country_map once (33 EMBI countries)
#  - Builds exposure_dt + avg_exposure_dt from Comtrade + IMF GDP
#  - Saves merge-ready exposures:
#      data/processed/exposure/exposure_clean_annual.rds
#      data/processed/exposure/exposure_clean_avg.rds
############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(lubridate)
})

#========================
# 0) PATHS (EDIT ONLY THESE)
#========================
MASTERFILE <- "data/data_masterfile.xlsx"

RDS_DIR <- "data/processed/masterfile_rds"
dir.create(RDS_DIR, recursive = TRUE, showWarnings = FALSE)

GDP_FILE   <- "data/imf/weo_gdp_ngdpd_annual_2000_2024.csv"
COMTRADE_TOTAL_FILE <- "data/comtrade/combined_total_exports_2000_2024.csv"
COMTRADE_SEL_FILE   <- "data/comtrade/combined_selected_products_2000_2024.csv"

OUT_EXP_DIR <- "data/processed/exposure"
dir.create(OUT_EXP_DIR, recursive = TRUE, showWarnings = FALSE)

# Rebuild cleaned RDS even if they exist?
FORCE_REBUILD_RDS <- FALSE

# Only clean the sheets you actually use
SHEETS_NEEDED <- c("EMBIs", "FXRates", "Commodities", "Global_Controls")

#========================
# 1) HELPERS: DATE + NUMERIC (robust)
#========================
is_bad_string <- function(x) {
  x <- trimws(toupper(as.character(x)))
  x %in% c("", "NA", "N/A", "NULL", "#N/A", "#N/A N/A", "REQUESTING DATA...", "#NAME?", "#VALUE!")
}

# Robust numeric converter:
# - leaves numeric alone
# - handles "1,234.56" and "1.234,56"
# - strips spaces
numeric_conversion <- function(x) {
  if (is.numeric(x)) return(x)
  
  x0 <- as.character(x)
  x0 <- trimws(x0)
  x0[is_bad_string(x0)] <- NA_character_
  
  # if it already looks like plain number with dot decimal and optional commas thousands, try direct
  suppressWarnings({
    direct <- as.numeric(gsub(",", "", x0))
  })
  
  # fallback: European style (thousands "." and decimal ",")
  # Strategy:
  #  - if string contains both "." and "," and comma is likely decimal -> remove "." then replace "," with "."
  x1 <- x0
  has_dot <- grepl("\\.", x1)
  has_com <- grepl(",", x1)
  
  euro <- x1
  euro[has_dot & has_com] <- gsub("\\.", "", euro[has_dot & has_com])
  euro <- gsub(",", ".", euro)
  
  suppressWarnings({
    euro_num <- as.numeric(euro)
  })
  
  # choose better (fewer NAs)
  out <- direct
  out[is.na(out)] <- euro_num[is.na(out)]
  out
}

# Robust date parser:
# - works for Date / POSIX
# - works for Excel serial numbers
# - works for "dd.mm.yyyy", ISO, etc.
parse_any_date <- function(x) {
  if (inherits(x, "Date")) return(x)
  if (inherits(x, "POSIXct") || inherits(x, "POSIXt")) return(as.Date(x))
  
  # Excel serial date?
  if (is.numeric(x)) {
    # Excel origin (Windows): 1899-12-30
    return(as.Date(x, origin = "1899-12-30"))
  }
  
  xs <- trimws(as.character(x))
  xs[is_bad_string(xs)] <- NA_character_
  
  d <- suppressWarnings(as.Date(xs, format = "%d.%m.%Y"))
  if (all(is.na(d))) d <- suppressWarnings(as.Date(xs, format = "%Y-%m-%d"))
  if (all(is.na(d))) d <- suppressWarnings(as.Date(xs, format = "%d/%m/%Y"))
  if (all(is.na(d))) d <- suppressWarnings(as.Date(xs, format = "%m/%d/%Y"))
  
  if (all(is.na(d))) {
    d <- suppressWarnings(lubridate::as_date(parse_date_time(xs, orders = c("dmy", "ymd", "mdy"))))
  }
  d
}

clean_sheet_masterfile <- function(path, sheet) {
  dt <- as.data.table(readxl::read_excel(path, sheet = sheet, .name_repair = "unique"))
  
  # Standardize date column name
  if ("Dates" %in% names(dt)) setnames(dt, "Dates", "date")
  if ("observation_date" %in% names(dt)) setnames(dt, "observation_date", "date")
  if (!"date" %in% names(dt)) setnames(dt, names(dt)[1], "date")
  
  # Parse date robustly (this fixes the “0 rows” problem)
  dt[, date := parse_any_date(date)]
  dt <- dt[!is.na(date)]
  
  # Clean column names (non-date)
  value_cols <- setdiff(names(dt), "date")
  if (length(value_cols)) {
    clean_names <- gsub("\\s+", "_", value_cols)
    clean_names <- gsub("[^A-Za-z0-9_]", "", clean_names)
    setnames(dt, value_cols, clean_names)
    
    # Convert to numeric where possible (keep characters if you ever have non-numeric columns—rare here)
    value_cols <- setdiff(names(dt), "date")
    dt[, (value_cols) := lapply(.SD, numeric_conversion), .SDcols = value_cols]
  }
  
  setorder(dt, date)
  setkey(dt, date)
  dt
}

rds_path <- function(name) file.path(RDS_DIR, paste0(name, ".rds"))

#========================
# 2) CLEAN XLSX -> SAVE SHEET RDS (ONLY ONCE)
#========================
if (!file.exists(MASTERFILE)) stop("Masterfile not found: ", MASTERFILE)

need_rebuild <- FORCE_REBUILD_RDS || any(!file.exists(rds_path(SHEETS_NEEDED)))
if (need_rebuild) {
  message("[ingest] Cleaning masterfile sheets -> RDS: ", RDS_DIR)
  
  available <- excel_sheets(MASTERFILE)
  missing_sheets <- setdiff(SHEETS_NEEDED, available)
  if (length(missing_sheets)) stop("Missing sheets in masterfile: ", paste(missing_sheets, collapse = ", "))
  
  for (s in SHEETS_NEEDED) {
    dt_clean <- clean_sheet_masterfile(MASTERFILE, s)
    saveRDS(dt_clean, rds_path(s))
    message("  - saved ", s, " (rows=", nrow(dt_clean), ", cols=", ncol(dt_clean), ")")
  }
}

#========================
# 3) LOAD CLEANED SHEETS (RDS)
#========================
EMBIs_dt           <- readRDS(rds_path("EMBIs"))
FXRates_dt         <- readRDS(rds_path("FXRates"))
Commodities_dt     <- readRDS(rds_path("Commodities"))
Global_Controls_dt <- readRDS(rds_path("Global_Controls"))

setDT(EMBIs_dt); setDT(FXRates_dt); setDT(Commodities_dt); setDT(Global_Controls_dt)

stopifnot("date" %in% names(EMBIs_dt),
          "date" %in% names(FXRates_dt),
          "date" %in% names(Commodities_dt),
          "date" %in% names(Global_Controls_dt))

# One dummy for USDized economies
if (!"USD_DUMMY" %in% names(FXRates_dt)) FXRates_dt[, USD_DUMMY := 1]

#========================
# 4) COUNTRY MAP (33 EMBI COUNTRIES) — DEFINE ONCE
#========================
country_map <- data.table(
  embi_ticker = c(
    "JPGCAR","JPGCBZ","JPGCCH","JPGCCO","JPGCEC","JPGCSV","JPGCMX","JPGCPN","JPGCPR","JPGCUR","JPGCVN",
    "JPGCAO","JPGCEG","JPGCMOR","JPGCNG","JPGCSN","JPGCSAF","JPGCZM","JPGCKZ","JPGCPL","JPGCRO","JPGCTU",
    "JPGCHU","JPGCBH","JPGCIQ","JPGCQA","JPGCSA","JPGCAE","JPGCID","JPGCPK","JPGCPH","JPGCLK","JPGCVI"
  ),
  country = c(
    "Argentina","Brazil","Chile","Colombia","Ecuador","El Salvador","Mexico","Panama","Peru","Uruguay","Vietnam",
    "Angola","Egypt","Morocco","Nigeria","Senegal","South Africa","Zambia","Kazakhstan","Poland","Romania","Turkey",
    "Hungary","Bahrain","Iraq","Qatar","Saudi Arabia","United Arab Emirates","Indonesia","Pakistan","Philippines","Sri Lanka","Ivory Coast"
  ),
  iso3 = c(
    "ARG","BRA","CHL","COL","ECU","SLV","MEX","PAN","PER","URY","VNM",
    "AGO","EGY","MAR","NGA","SEN","ZAF","ZMB","KAZ","POL","ROU","TUR",
    "HUN","BHR","IRQ","QAT","SAU","ARE","IDN","PAK","PHL","LKA","CIV"
  )
)
country_map[, embi_col := paste0(embi_ticker, "_Index")]

ccy_by_iso3 <- c(
  ARG="ARS", BRA="BRL", CHL="CLP", COL="COP", ECU="USD", SLV="USD", MEX="MXN",
  PAN="USD", PER="PEN", URY="UYU", VNM="VND", AGO="AOA", EGY="EGP", MAR="MAD",
  NGA="NGN", SEN="XOF", ZAF="ZAR", ZMB="ZMW", KAZ="KZT", POL="PLN", ROU="EUR",
  TUR="TRY", HUN="HUF", BHR="BHD", IRQ="IQD", QAT="QAR", SAU="SAR", ARE="AED",
  IDN="IDR", PAK="PKR", PHL="PHP", LKA="LKR", CIV="XOF"
)

country_map[, ccy := ccy_by_iso3[iso3]]
country_map[, fx_col := ifelse(ccy == "USD", "USD_DUMMY", paste0("USD", ccy, "_Curncy"))]

#========================
# 5) EXPOSURES (Comtrade + GDP)
#========================
pick_value_col <- function(dt) {
  for (cc in c("primary_value","fobvalue","cifvalue","Trade.Value..US..","trade_value_usd")) {
    if (cc %in% names(dt)) return(cc)
  }
  stop("No recognized value column in Comtrade file.")
}

basket_sum <- function(x_codes, x_vals, prefixes) {
  if (length(x_codes) == 0) return(0)
  hit <- Reduce(`|`, lapply(prefixes, function(p) startsWith(x_codes, p)))
  sum(x_vals[hit], na.rm = TRUE)
}

# --- GDP
gdp <- fread(GDP_FILE)
setnames(gdp, tolower(names(gdp)))
stopifnot(all(c("country","year") %in% names(gdp)))
if (!"value" %in% names(gdp)) stop("GDP file must contain column 'value'.")

setnames(gdp, "value", "gdp_usd")
gdp[, iso3 := toupper(country)]
gdp[, year := as.integer(year)]
gdp[, gdp_usd := as.numeric(gdp_usd)]
gdp <- gdp[, .(iso3, year, gdp_usd)]

# --- Comtrade totals + selected
tot <- fread(COMTRADE_TOTAL_FILE)
sel <- fread(COMTRADE_SEL_FILE)

tot_val <- pick_value_col(tot)
sel_val <- pick_value_col(sel)

tot2 <- tot[, .(
  iso3 = toupper(reporter_iso),
  year = as.integer(ref_year),
  total_exports_usd = as.numeric(get(tot_val))
)]
sel2 <- sel[, .(
  iso3 = toupper(reporter_iso),
  year = as.integer(ref_year),
  hs_code = gsub("\\D", "", as.character(cmd_code)),
  exports_usd = as.numeric(get(sel_val))
)]

tot2 <- tot2[!is.na(iso3) & nchar(iso3) == 3 & !is.na(year)]
sel2 <- sel2[!is.na(iso3) & nchar(iso3) == 3 & !is.na(year) & !is.na(hs_code) & hs_code != ""]

# --- HS baskets (you can expand later)
HS_MAP <- list(
  crude_oil      = c("2709"),
  refined_oil    = c("2710"),
  nat_gas        = c("2711"),
  energy_total   = c("27"),
  
  copper_ore     = c("2603"),
  copper_total   = c("74"),
  aluminum       = c("76"),
  nickel         = c("75"),
  iron_ore       = c("2601"),
  iron_steel     = c("72","73"),
  
  gold           = c("7108"),
  precious_total = c("71"),
  
  wheat          = c("1001"),
  corn           = c("1005"),
  soybeans       = c("1201"),
  cereals_total  = c("10"),
  oilseeds_total = c("12"),
  
  live_cattle     = c("0102"),
  beef_total      = c("0201","0202"),
  livestock_total = c("01","02")
)

base <- unique(tot2[, .(iso3, year, total_exports_usd)])
setorder(base, iso3, year)

setkey(sel2, iso3, year)

exposure_dt <- base[, {
  dty <- sel2[.BY, nomatch = 0L]
  codes <- dty$hs_code
  vals  <- dty$exports_usd
  
  row <- vector("list", length(HS_MAP))
  names(row) <- paste0("exports_", names(HS_MAP), "_usd")
  
  for (nm in names(HS_MAP)) {
    row[[paste0("exports_", nm, "_usd")]] <- basket_sum(codes, vals, HS_MAP[[nm]])
  }
  as.data.table(row)
}, by = .(iso3, year, total_exports_usd)]

# Merge GDP
exposure_dt <- merge(exposure_dt, gdp, by = c("iso3","year"), all.x = TRUE)

# Shares + exports/GDP
export_cols <- grep("^exports_.*_usd$", names(exposure_dt), value = TRUE)
for (ec in export_cols) {
  short <- sub("^exports_", "", ec)
  short <- sub("_usd$", "", short)
  
  exposure_dt[, (paste0("share_", short)) :=
                fifelse(total_exports_usd > 0, get(ec) / total_exports_usd, NA_real_)]
  
  exposure_dt[, (paste0("expgdp_", short)) :=
                fifelse(gdp_usd > 0, get(ec) / gdp_usd, NA_real_)]
}

# Attach full name from country_map (iso3 stays key)
exposure_dt <- merge(exposure_dt,
                     unique(country_map[, .(iso3, country)]),
                     by = "iso3", all.x = TRUE)

setcolorder(exposure_dt, c("country","iso3","year", setdiff(names(exposure_dt), c("country","iso3","year"))))
setorder(exposure_dt, iso3, year)

# Average exposure table
share_cols  <- grep("^share_",  names(exposure_dt), value = TRUE)
expgdp_cols <- grep("^expgdp_", names(exposure_dt), value = TRUE)

avg_exposure_dt <- exposure_dt[, c(
  .(years_obs = .N),
  lapply(.SD, function(x) mean(x, na.rm = TRUE))
), by = iso3, .SDcols = c(share_cols, expgdp_cols)]

avg_exposure_dt <- merge(avg_exposure_dt,
                         unique(country_map[, .(iso3, country)]),
                         by = "iso3", all.x = TRUE)

setcolorder(avg_exposure_dt, c("country","iso3","years_obs",
                               setdiff(names(avg_exposure_dt), c("country","iso3","years_obs"))))
setorder(avg_exposure_dt, iso3)

#========================
# 6) SAVE MERGE-READY EXPOSURE TABLES (RDS + CSV)
#========================
# Keep only columns you actually need for regressions (small & stable)
keep_exp_cols <- c(
  "share_crude_oil", "expgdp_crude_oil",
  "share_copper_total", "expgdp_copper_total"
)
keep_exp_cols <- intersect(keep_exp_cols, names(exposure_dt))

exposure_clean <- exposure_dt[, c("iso3", "country", "year", keep_exp_cols), with = FALSE]
setkey(exposure_clean, iso3, year)

avg_keep_cols <- c("years_obs", grep("^share_|^expgdp_", names(avg_exposure_dt), value = TRUE))
avg_exposure_clean <- avg_exposure_dt[, c("iso3", "country", avg_keep_cols), with = FALSE]
setkey(avg_exposure_clean, iso3)

saveRDS(exposure_clean, file.path(OUT_EXP_DIR, "exposure_clean_annual.rds"))
saveRDS(avg_exposure_clean, file.path(OUT_EXP_DIR, "exposure_clean_avg.rds"))

# optional human-readable
fwrite(exposure_clean, file.path(OUT_EXP_DIR, "exposure_clean_annual.csv"))
fwrite(avg_exposure_clean, file.path(OUT_EXP_DIR, "exposure_clean_avg.csv"))

#========================
# 7) CLEANUP BIG RAW OBJECTS (optional but good)
#========================
rm(tot, sel, tot2, sel2, gdp, base)
invisible(gc())

#========================
# 8) FINAL MESSAGE
#========================
message("Loaded: EMBIs_dt, FXRates_dt, Global_Controls_dt, Commodities_dt, country_map")
message("Built: exposure_dt, avg_exposure_dt")
message("Saved: ", file.path(OUT_EXP_DIR, "exposure_clean_annual.rds"))
message("Saved: ", file.path(OUT_EXP_DIR, "exposure_clean_avg.rds"))
message("Exposure years: ", min(exposure_dt$year, na.rm=TRUE), " - ", max(exposure_dt$year, na.rm=TRUE))
message("Country map rows: ", nrow(country_map))