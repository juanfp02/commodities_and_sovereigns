############################################################
# imf_sdmx_weo_gdp_annual.R
# Pull IMF WEO GDP current prices (NGDPD) annual, 2000–2024
############################################################

suppressPackageStartupMessages({
  library(rsdmx)
  library(data.table)
})

# --- Your EMBI country set (ISO3)
REPORTERS_ISO3 <- c(
  "ARG","BRA","CHL","COL","ECU","SLV","MEX","PAN","PER","URY","VNM",
  "AGO","EGY","MAR","NGA","SEN","ZAF","ZMB",
  "KAZ","POL","ROU","TUR","HUN",
  "BHR","IRQ","QAT","SAU","ARE",
  "IDN","PAK","PHL","LKA","CIV"
)

# --- WEO SDMX settings
flowref <- "IMF.RES,WEO"    # WEO dataset in IMF Data portal
indicator <- "NGDPD"        # GDP, current prices (WEO)
freq <- "A"
start_year <- 2000
end_year   <- 2024

# SDMX key format for WEO is COUNTRY.INDICATOR.FREQUENCY
# Multiple countries in one dimension are separated by '+'
key <- paste0(paste(REPORTERS_ISO3, collapse = "+"), ".", indicator, ".", freq)

# --- Pull data (public)
sx <- readSDMX(
  providerId = "IMF_DATA",
  resource = "data",
  flowRef = flowref,
  key = key,
  start = start_year,
  end   = end_year
)

df <- as.data.frame(sx)
dt <- as.data.table(df)

# --- Standardize column names (rsdmx names can vary slightly)
# Typical columns include: COUNTRY, INDICATOR, FREQUENCY, TIME_PERIOD, obsValue
# We'll defensively map what exists.
cn <- names(dt)

# Find time and value columns
time_col <- cn[cn %in% c("TIME_PERIOD", "time", "Time", "obsTime")][1]
val_col  <- cn[cn %in% c("obsValue", "OBS_VALUE", "value", "Value")][1]

# Find country/indicator columns
cty_col  <- cn[cn %in% c("COUNTRY", "REF_AREA", "ref_area")][1]
ind_col  <- cn[cn %in% c("INDICATOR", "INDICATOR_CODE", "indicator")][1]

# Minimal cleanup
setnames(dt, c(cty_col, ind_col, time_col, val_col), c("country", "indicator", "year", "value"), skip_absent = TRUE)
dt[, year := as.integer(year)]
dt[, value := as.numeric(value)]

# Keep only what you need
gdp_dt <- dt[, .(country, year, indicator, value)][order(country, year)]

# Save
dir.create("data/imf", recursive = TRUE, showWarnings = FALSE)
fwrite(gdp_dt, "data/imf/weo_gdp_ngdpd_annual_2000_2024.csv")

