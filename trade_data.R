############################################################
# comtrade_pull_2000_2024.R
# Pull UN Comtrade exports for:
#  - TOTAL exports (all products)
#  - Selected products (HS codes placeholder)
#
# Output:
#  - data/comtrade/raw/total_exports_<year>.csv
#  - data/comtrade/raw/selected_products_<year>.csv
#  - data/comtrade/combined_total_exports_2000_2024.csv
#  - data/comtrade/combined_selected_products_2000_2024.csv
############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(comtradr)
  library(purrr)
})

#========================
# 0) SETUP
#========================

# Put your Comtrade API key here (recommended if you have one)
# If you don't have a key, leave it NULL (you will hit stricter limits)
COMTRADE_API_KEY <- Sys.getenv("COMTRADE_PRIMARY")
if (COMTRADE_API_KEY != "") {
  comtradr::set_primary_comtrade_key(COMTRADE_API_KEY)
}

YEARS <- 2017:2024

# Reporter = "all" means all reporters; can be very large.
# Consider using a reporter list if you hit limits.
REPORTER <- c(
  "ARG", # JPGCAR  Argentina
  "BRA", # JPGCBZ  Brazil
  "CHL", # JPGCCH  Chile
  "COL", # JPGCCO  Colombia
  "ECU", # JPGCEC  Ecuador
  "SLV", # JPGCSV  El Salvador
  "MEX", # JPGCMX  Mexico
  "PAN", # JPGCPN  Panama
  "PER", # JPGCPR  Peru
  "URY", # JPGCUR  Uruguay
  "VNM", # JPGCVN  Vietnam
  
  "AGO", # JPGCAO  Angola
  "EGY", # JPGCEG  Egypt
  "MAR", # JPGCMOR Morocco
  "NGA", # JPGCNG  Nigeria
  "SEN", # JPGCSN  Senegal
  "ZAF", # JPGCSAF South Africa
  "ZMB", # JPGCZM  Zambia
  
  "KAZ", # JPGCKZ  Kazakhstan
  "POL", # JPGCPL  Poland
  "ROU", # JPGCRO  Romania
  "TUR", # JPGCTU  Turkey
  "HUN", # JPGCHU  Hungary
  
  "BHR", # JPGCBH  Bahrain
  "IRQ", # JPGCIQ  Iraq
  "QAT", # JPGCQA  Qatar
  "SAU", # JPGCSA  Saudi Arabia
  "ARE", # JPGCAE  United Arab Emirates
  
  "IDN", # JPGCID  Indonesia
  "PAK", # JPGCPK  Pakistan
  "PHL", # JPGCPH  Philippines
  "LKA", # JPGCLK  Sri Lanka
  
  "CIV"  # JPGCVI  Côte d’Ivoire (Ivory Coast)
)

# Partner: "World" (WLD) is typical when you want total exports.
PARTNER <- "World"

# Trade flow: "Export"
FLOW <- "Export"

# Frequency: Annual
FREQ <- "A"

# Classification: HS (choose version)
# If you want HS as reported each year, keep it simple: "HS"
CLASSIFICATION <- "HS"

# Product codes placeholder: put HS codes you care about
# Example for oil crude (2709) and refined (2710), copper ore (2603), copper cathodes (7403) etc.
PRODUCT_CODES <- c(
  
  # ========================
  # Energy (CO1, CL1, BCOMEN)
  # ========================
  "27",    # Mineral fuels, oils, waxes (broad energy chapter)
  "2709",  # Crude petroleum oils
  "2710",  # Petroleum oils (refined) + preparations
  
  # (Optional if you want gas explicitly; still usually within 27)
  "2711",  # Petroleum gases and other gaseous hydrocarbons (LNG/LPG etc.)
  
  # ========================
  # Industrial metals (HG1, LA1, LN1, HRC1, IOE1)
  # ========================
  "74",    # Copper and articles thereof (HG1 proxy basket)
  "2603",  # Copper ores and concentrates (important for exporter exposure)
  
  "76",    # Aluminum and articles thereof (LA1)
  "2606",
  
  "75",    # Nickel and articles thereof (LN1)
  "2603",
  
  "72",    # Iron and steel (HRC1 steel proxy)
  "73",    # Articles of iron or steel (optional but often relevant to steel exports)
  "2601",  # Iron ores and concentrates (IOE1)
  
  # ========================
  # Precious metals (GC1, BCOMPR)
  # ========================
  "71",    # Natural/cultured pearls, precious stones, precious metals
  "7108",  # Gold (incl. plated), unwrought or semi-manufactured
  
  # ========================
  # Agriculture (W1, S1, C1, BCOMAG)
  # ========================
  "10",    # Cereals (broad: wheat, corn, etc.)
  "1001",  # Wheat and meslin
  "1005",  # Maize (corn)
  
  "12",    # Oil seeds and oleaginous fruits (broad)
  "1201",  # Soybeans
  
  # ========================
  # Livestock / meat (FC1, LC1, CVT1)
  # ========================
  "01",    # Live animals
  "0102",  # Live bovine animals
  "02",    # Meat and edible meat offal
  "0201",  # Beef, fresh or chilled
  "0202"   # Beef, frozen
)

# Optionally: pull ALL products totals (TOTAL)
PULL_TOTAL_EXPORTS <- TRUE

# Output folders
dir.create("data/comtrade/raw", recursive = TRUE, showWarnings = FALSE)
dir.create("data/comtrade", recursive = TRUE, showWarnings = FALSE)

#========================
# 1) HELPER: safe query + pagination
#========================

pull_comtrade <- function(year,
                          product_codes = NULL,
                          pull_total = FALSE) {
  
  if (pull_total) {
    # TOTAL exports (all products)
    message("Pulling TOTAL exports for year ", year, " ...")
    res <- ct_get_data(
      reporter = REPORTER,
      partner  = PARTNER,
      flow_direction = FLOW,
      start_date = year,
      end_date   = year,
      frequency  = FREQ,
      commodity_classification = CLASSIFICATION,
      commodity_code = "TOTAL"
    )
    dt <- as.data.table(res)
    dt[, query_type := "TOTAL"]
    return(dt)
  }
  
  # Selected product codes
  if (is.null(product_codes) || length(product_codes) == 0) {
    stop("product_codes is empty. Provide HS codes or use pull_total=TRUE.")
  }
  
  # comtradr supports multiple commodity_code values, but large lists can fail.
  # We'll loop codes and bind (more robust).
  message("Pulling selected products for year ", year, " ...")
  
  pieces <- lapply(product_codes, function(cc) {
    message("  - code ", cc)
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
# 2) PULL + SAVE PER YEAR
#========================

all_total <- list()
all_sel   <- list()

for (yr in YEARS) {
  
  # --- TOTAL exports ---
  if (PULL_TOTAL_EXPORTS) {
    dt_total <- tryCatch(
      pull_comtrade(year = yr, pull_total = TRUE),
      error = function(e) {
        message("ERROR total year ", yr, ": ", e$message)
        NULL
      }
    )
    
    if (!is.null(dt_total)) {
      fwrite(dt_total, file = sprintf("data/comtrade/raw/total_exports_%d.csv", yr))
      all_total[[as.character(yr)]] <- dt_total
    }
  }
  
  # --- Selected products ---
  dt_sel <- tryCatch(
    pull_comtrade(year = yr, product_codes = PRODUCT_CODES, pull_total = FALSE),
    error = function(e) {
      message("ERROR selected year ", yr, ": ", e$message)
      NULL
    }
  )
  
  if (!is.null(dt_sel)) {
    fwrite(dt_sel, file = sprintf("data/comtrade/raw/selected_products_%d.csv", yr))
    all_sel[[as.character(yr)]] <- dt_sel
  }
  
  # Small sleep to reduce rate-limit risk
  Sys.sleep(0.5)
}

#========================
# 3) COMBINE + SAVE
#========================
if (length(all_total) > 0) {
  dt_total_all <- rbindlist(all_total, fill = TRUE)
  fwrite(dt_total_all, "data/comtrade/combined_total_exports_2000_2024.csv")
}

if (length(all_sel) > 0) {
  dt_sel_all <- rbindlist(all_sel, fill = TRUE)
  fwrite(dt_sel_all, "data/comtrade/combined_selected_products_2000_2024.csv")
}

message("DONE.")
RAW_DIR <- "data/comtrade/raw"     # where the yearly files are
OUT_DIR <- "data/comtrade"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Patterns for your yearly files
PAT_TOTAL <- "^total_exports_\\d{4}\\.csv$"
PAT_SEL   <- "^selected_products_\\d{4}\\.csv$"

#========================
# 1) Helper: merge all files matching a pattern
#========================
merge_pattern <- function(pattern, out_file) {
  files <- list.files(RAW_DIR, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) stop("No files found for pattern: ", pattern)
  
  # Sort by year extracted from filename
  get_year <- function(fp) as.integer(gsub("\\D", "", basename(fp)))
  files <- files[order(sapply(files, get_year))]
  
  message("Merging ", length(files), " files into: ", out_file)
  
  dt_all <- rbindlist(lapply(files, function(f) {
    dt <- fread(f)
    
    # If ref_year missing (shouldn't be), create from filename
    if (!("ref_year" %in% names(dt))) {
      dt[, ref_year := get_year(f)]
    }
    dt
  }), fill = TRUE)
  
  fwrite(dt_all, file.path(OUT_DIR, out_file))
  invisible(dt_all)
}

#========================
# 2) Run merges
#========================
dt_total <- merge_pattern(PAT_TOTAL, "combined_total_exports_2000_2024.csv")
dt_sel   <- merge_pattern(PAT_SEL,   "combined_selected_products_2000_2024.csv")

message("DONE.")
message("Rows TOTAL: ", nrow(dt_total))
message("Rows SELECTED: ", nrow(dt_sel))
