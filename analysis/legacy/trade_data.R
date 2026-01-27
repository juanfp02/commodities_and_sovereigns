############################################################
# comtrade_pull_2000_2024.R
# Pull UN Comtrade exports for SELECTED products only (HS4)
#
# Output:
#  - data/comtrade/raw/selected_products_<year>.csv
#  - data/comtrade/combined_selected_products_2000_2024.csv
############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(comtradr)
})

#========================
# 0) SETUP
#========================

COMTRADE_API_KEY <- Sys.getenv("COMTRADE_PRIMARY")
if (COMTRADE_API_KEY != "") {
  comtradr::set_primary_comtrade_key(COMTRADE_API_KEY)
}

YEARS <- 2016:2024

REPORTER <- c(
  "ARG","BRA","CHL","COL","ECU","SLV","MEX","PAN","PER","URY","VNM",
  "AGO","EGY","MAR","NGA","SEN","ZAF","ZMB",
  "KAZ","POL","ROU","TUR","HUN",
  "BHR","IRQ","QAT","SAU","ARE",
  "IDN","PAK","PHL","LKA",
  "CIV"
)

PARTNER <- "World"
FLOW <- "Export"
FREQ <- "A"
CLASSIFICATION <- "HS"

#========================
# 1) HS4 CODE LISTS (NO DOUBLE COUNTING)
#========================

# Energy (BCOM Energy)
BCOM_ENERGY_CODES_HS4 <- c("2709","2710","2711")

# Agriculture (BCOM Agriculture)
BCOM_AGRICULTURE_CODES_HS4 <- c("0901","1001","1005","1201","1507","1701","2304","5201")

# Industrial Metals (BCOM Industrial Metals) — ore + refined/unwrought
BCOM_INDUSTRIAL_METALS_CODES_HS4 <- c(
  # Copper
  "2603","7402","7403",
  # Aluminum
  "2606","7601",
  # Nickel
  "2604","7502",
  # Zinc
  "2608","7901",
  # Lead
  "2607","7801"
)

# Precious Metals (BCOM Precious Metals)
BCOM_PRECIOUS_METALS_CODES_HS4 <- c("7106","7108")

# Combine + deduplicate
PRODUCT_CODES <- unique(c(
  BCOM_ENERGY_CODES_HS4,
  BCOM_AGRICULTURE_CODES_HS4,
  BCOM_INDUSTRIAL_METALS_CODES_HS4,
  BCOM_PRECIOUS_METALS_CODES_HS4
))

# Output folders
RAW_DIR <- "data/comtrade/raw"
OUT_DIR <- "data/comtrade"
dir.create(RAW_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

#========================
# 2) HELPER: Pull selected product codes for a year
#========================
pull_selected_for_year <- function(year, product_codes) {
  
  message("Pulling selected products for year ", year, " ...")
  
  pieces <- lapply(product_codes, function(cc) {
    message("  - HS ", cc)
    out <- ct_get_data(
      reporter = REPORTER,
      partner  = PARTNER,
      flow_direction = FLOW,
      start_date = year,
      end_date   = year,
      frequency  = FREQ,
      commodity_classification = CLASSIFICATION,
      commodity_code = cc
    )
    d <- as.data.table(out)
    d[, query_type := paste0("HS_", cc)]
    d
  })
  
  rbindlist(pieces, fill = TRUE)
}

#========================
# 3) PULL + SAVE PER YEAR
#========================
all_sel <- list()

for (yr in YEARS) {
  
  dt_sel <- tryCatch(
    pull_selected_for_year(year = yr, product_codes = PRODUCT_CODES),
    error = function(e) {
      message("ERROR selected year ", yr, ": ", e$message)
      NULL
    }
  )
  
  if (!is.null(dt_sel)) {
    fwrite(dt_sel, file = sprintf("%s/selected_products_%d.csv", RAW_DIR, yr))
    all_sel[[as.character(yr)]] <- dt_sel
  }
  
  # gentle rate-limit protection
  Sys.sleep(0.5)
}

#========================
# 4) COMBINE + SAVE
#========================
if (length(all_sel) > 0) {
  dt_sel_all <- rbindlist(all_sel, fill = TRUE)
  fwrite(dt_sel_all, file.path(OUT_DIR, "combined_selected_products_2000_2024.csv"))
  message("DONE. Wrote: ", file.path(OUT_DIR, "combined_selected_products_2000_2024.csv"))
} else {
  message("No data pulled. Check API key / limits / parameters.")
}

#========================
# 5) OPTIONAL: merge raw yearly files (if you ever re-run partially)
#========================
merge_pattern <- function(pattern, out_file) {
  files <- list.files(RAW_DIR, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) stop("No files found for pattern: ", pattern)
  
  get_year <- function(fp) as.integer(gsub("\\D", "", basename(fp)))
  files <- files[order(sapply(files, get_year))]
  
  message("Merging ", length(files), " files into: ", out_file)
  
  dt_all <- rbindlist(lapply(files, function(f) {
    dt <- fread(f)
    if (!("ref_year" %in% names(dt))) dt[, ref_year := get_year(f)]
    dt
  }), fill = TRUE)
  
  fwrite(dt_all, file.path(OUT_DIR, out_file))
  invisible(dt_all)
}

# Uncomment if you want to rebuild from raw folder later:
# dt_sel_all2 <- merge_pattern("^selected_products_\\d{4}\\.csv$", "combined_selected_products_2000_2024.csv")