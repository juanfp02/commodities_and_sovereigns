############################################################
# load_masterfile_excel.R
# Reads an Excel workbook and creates one data.table per sheet
############################################################

library(data.table)
library(lubridate)
library(readxl)

# ---------- helpers (re-use your numeric conversion) ----------
numeric_conversion <- function(x) {
  if (is.numeric(x)) return(x)
  
  x <- as.character(x)
  x <- trimws(x)
  
  # remove thousands sep "." and switch decimal "," -> "."
  x <- gsub("\\.", "", x)
  x <- gsub(",", ".", x)
  
  # treat empty / Bloomberg NA-ish values
  x[x %in% c("", "NA", "#N/A", "#N/A N/A", "N/A")] <- NA
  
  suppressWarnings(as.numeric(x))
}

parse_any_date <- function(x) {
  if (inherits(x, "Date")) return(x)
  
  x <- as.character(x)
  x <- trimws(x)
  
  # Try common formats (Bloomberg EU + ISO)
  d <- suppressWarnings(as.Date(x, format = "%d.%m.%Y"))
  if (all(is.na(d))) d <- suppressWarnings(as.Date(x, format = "%Y-%m-%d"))
  if (all(is.na(d))) d <- suppressWarnings(as.Date(x, format = "%d/%m/%Y"))
  if (all(is.na(d))) d <- suppressWarnings(as.Date(x, format = "%m/%d/%Y"))
  
  # fallback: lubridate guess
  if (all(is.na(d))) {
    d <- suppressWarnings(lubridate::as_date(parse_date_time(x, orders = c("dmy", "ymd", "mdy"))))
  }
  d
}

clean_sheet_dt <- function(path, sheet) {
  dt <- as.data.table(read_excel(path, sheet = sheet, .name_repair = "unique"))
  
  # standardize date column name
  if ("Dates" %in% names(dt)) setnames(dt, "Dates", "date")
  if ("observation_date" %in% names(dt)) setnames(dt, "observation_date", "date")
  
  # if no "date", try first column
  if (!"date" %in% names(dt) && ncol(dt) >= 1) {
    setnames(dt, names(dt)[1], "date")
  }
  
  # parse date
  dt[, date := parse_any_date(date)]
  dt <- dt[!is.na(date)]
  setkey(dt, date)
  
  # clean column names (except date)
  value_cols <- setdiff(names(dt), "date")
  clean_names <- gsub("\\s+", "_", value_cols)
  clean_names <- gsub("[^A-Za-z0-9_]", "", clean_names)
  setnames(dt, value_cols, clean_names)
  
  # numeric conversion on non-date columns
  value_cols <- setdiff(names(dt), "date")
  dt[, (value_cols) := lapply(.SD, numeric_conversion), .SDcols = value_cols]
  
  dt
}

# ---------- main ----------
MASTERFILE <- "data/data_masterfile.xlsx"   # change if you want local path

sheets <- excel_sheets(MASTERFILE)

# list of data.tables, named by sheet
dt_list <- setNames(lapply(sheets, function(s) clean_sheet_dt(MASTERFILE, s)),
                    make.names(sheets))

# also create individual objects in your environment: <sheet>_dt
for (nm in names(dt_list)) {
  assign(paste0(nm, "_dt"), dt_list[[nm]], envir = .GlobalEnv)
}

# Optional: save each cleaned sheet as .rds files
# dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
# for (nm in names(dt_list)) {
#   saveRDS(dt_list[[nm]], file = file.path("data/processed", paste0(nm, ".rds")))
# }

message("Loaded sheets: ", paste(names(dt_list), collapse = ", "))



#------ exposure
#############################################################
suppressPackageStartupMessages({
  library(data.table)
})

#========================
# PATHS
#========================
GDP_FILE   <- "data/imf/weo_gdp_ngdpd_annual_2000_2024.csv"
COMTRADE_TOTAL_FILE <- "data/comtrade/combined_total_exports_2000_2024.csv"
COMTRADE_SEL_FILE   <- "data/comtrade/combined_selected_products_2000_2024.csv"

OUT_DIR <- "data/exposure"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

#========================
# 1) LOAD GDP (GDP already in USD!)
#========================
gdp <- fread(GDP_FILE)
setnames(gdp, tolower(names(gdp)))

stopifnot(all(c("country","year") %in% names(gdp)))
val_col <- intersect(names(gdp), c("value","gdp_raw","gdp","obs_value","obsvalue"))[1]
if (is.na(val_col)) stop("GDP file missing GDP numeric column.")

setnames(gdp, val_col, "gdp_usd")
gdp[, country := toupper(country)]
gdp[, year := as.integer(year)]
gdp[, gdp_usd := as.numeric(gdp_usd)]
gdp <- gdp[, .(country, year, gdp_usd)]

#========================
# 2) LOAD COMTRADE TOTAL + SELECTED (your column names)
#========================
tot <- fread(COMTRADE_TOTAL_FILE)
sel <- fread(COMTRADE_SEL_FILE)

tot2 <- tot[, .(
  country = toupper(reporter_iso),
  year    = as.integer(ref_year),
  total_exports_usd = as.numeric(primary_value)
)]

sel2 <- sel[, .(
  country = toupper(reporter_iso),
  year    = as.integer(ref_year),
  hs_code = gsub("\\D", "", as.character(cmd_code)),  # digits only
  exports_usd = as.numeric(primary_value)
)]

#========================
# 3) DEFINE HS BASKETS (match your Bloomberg commodities)
#========================
# IMPORTANT: Comtrade product codes are HS. For futures like HRC/FeederCattle,
# we map to reasonable HS categories that represent those commodities in trade.

HS_MAP <- list(
  
  # ---- Energy (CO1, CL1, BCOMEN)
  oil_crude      = c("2709"),
  oil_refined    = c("2710"),
  nat_gas_lng    = c("2711"),
  energy_total   = c("27"),      # broad energy chapter
  
  # ---- Industrial metals (HG1, LA1, LN1, HRC1, IOE1)
  copper_ore     = c("2603"),
  copper_metal   = c("74"),      # copper & articles
  aluminum       = c("76"),
  nickel         = c("75"),
  steel_primary  = c("72"),      # iron & steel
  steel_articles = c("73"),
  iron_ore       = c("2601"),
  
  # ---- Precious metals (GC1, BCOMPR)
  gold           = c("7108"),
  precious_total = c("71"),
  
  # ---- Agriculture (W1, S1, C1, BCOMAG)
  wheat          = c("1001"),
  corn           = c("1005"),
  cereals_total  = c("10"),
  soybeans       = c("1201"),
  oilseeds_total = c("12"),
  
  # ---- Livestock (FC1, LC1, CVT1)
  live_cattle    = c("0102"),
  meat_beef_fresh= c("0201"),
  meat_beef_frozen=c("0202"),
  livestock_total= c("01","02")
)

# Helper to compute exports for each basket by startsWith logic
basket_exports <- function(dt, code_prefixes) {
  # dt: sel2 for a given country-year
  # code_prefixes: vector of HS prefixes like "27", "2709"
  dt[Reduce(`|`, lapply(code_prefixes, function(p) startsWith(hs_code, p))),
     sum(exports_usd, na.rm = TRUE)]
}

#========================
# 4) BUILD COUNTRY-YEAR TABLE OF EXPORTS PER BASKET
#========================
# First create the list of all country-years that appear in totals
base <- unique(tot2[, .(country, year, total_exports_usd)])

# Compute basket exports by country-year
# (This is efficient enough for your country set + 2000-2024)
out_list <- base[, {
  dty <- sel2[country == .BY$country & year == .BY$year]
  row <- list()
  for (nm in names(HS_MAP)) {
    row[[paste0("exports_", nm, "_usd")]] <- basket_exports(dty, HS_MAP[[nm]])
  }
  as.data.table(row)
}, by = .(country, year)]

# Merge total exports back in
out <- merge(base, out_list, by = c("country","year"), all.x = TRUE)
for (cc in names(out)) if (is.numeric(out[[cc]])) out[is.na(get(cc)), (cc) := 0]

# Merge GDP
out <- merge(out, gdp, by = c("country","year"), all.x = TRUE)

#========================
# 5) EXPOSURES: share of exports + exports/GDP
#========================
export_cols <- grep("^exports_.*_usd$", names(out), value = TRUE)

# share of total exports
for (ec in export_cols) {
  nm <- sub("^exports_", "share_", ec)
  nm <- sub("_usd$", "", nm)
  out[, (nm) := fifelse(total_exports_usd > 0, get(ec) / total_exports_usd, NA_real_)]
}

# exports / GDP
for (ec in export_cols) {
  nm <- sub("^exports_", "expgdp_", ec)
  nm <- sub("_usd$", "", nm)
  out[, (nm) := fifelse(gdp_usd > 0, get(ec) / gdp_usd, NA_real_)]
}

setorder(out, country, year)

#========================
# 6) SAVE
#========================
fwrite(out, file.path(OUT_DIR, "exposure_country_year_2000_2024_full.csv"))
saveRDS(out, file.path(OUT_DIR, "exposure_country_year_2000_2024_full.rds"))

cat("DONE: exposures written to:\n",
    file.path(OUT_DIR, "exposure_country_year_2000_2024_full.csv"), "\n")


rm(c(avg_expgdp, avg_share, avg_share_w, base, embi_map, expo, fx_map, gdp, HS_MAP, out_list, panel_levels, sel,sel2, tot, tot2))
