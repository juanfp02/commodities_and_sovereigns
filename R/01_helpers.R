library(data.table)
library(lubridate)


  FX        <- "data/FX_Rates.csv"
  EMBI      <- "data/EMBI.csv"
  DGS10     <- "data/DGS10.csv"
  DGS2      <- "data/DGS2.csv"
  VIX       <- "data/VIXCLS.csv"
  Indices   <- "data/Indices.csv"
  Commod   <- "data/Commodity_prices.csv"
  Brent    <- "data/DCOILBRENTEU.csv"
  CDS      <- "data/CDS.csv"
  USDIndex <- "data/DTWEXBGS.csv"
  
numeric_conversion <- function(x) {
    
    # If already numeric, return as is
    if (is.numeric(x)) return(x)
    
    # Work on character
    x <- as.character(x)
    
    # Remove non-numeric junk
    x <- trimws(x)
    
    # Convert European format:
    # 1. Remove thousands separator "."
    # 2. Replace decimal comma "," with "."
    x <- gsub("\\.", "", x)
    x <- gsub(",", ".", x)
    
    as.numeric(x)
}


clean_bloomberg_dt <- function(file) {
  
  dt <- fread(
    file,
    sep = ";",
    na.strings = c("", "NA", "#N/A", "#N/A N/A")
  )
  
  # Rename date column
  setnames(dt, "Dates", "date")
  
  # Parse European dates
  dt[, date := as.Date(date, format = "%d.%m.%Y")]
  dt <- dt[!is.na(date)]
  
  # Identify value columns
  value_cols <- setdiff(names(dt), "date")
  
  # ---- CORRECT NUMERIC HANDLING ----
  dt[, (value_cols) := lapply(.SD, numeric_conversion), .SDcols = value_cols]
  
  # Clean column names
  clean_names <- gsub(" ", "_", value_cols)
  setnames(dt, value_cols, clean_names)
  
  setkey(dt, date)
  return(dt)
}

clean_FRED_dt <- function(file) {
  
  dt <- fread(file, na.strings = c("", "NA", "."))
  
  # Rename date column
  if ("observation_date" %in% names(dt)) {
    setnames(dt, "observation_date", "date")
  } else {
    stop("No 'observation_date' column found")
  }
  
  # Ensure Date class (FRED uses ISO format)
  dt[, date := as.Date(date)]
  
  # Drop invalid dates (safety)
  dt <- dt[!is.na(date)]
  
  # Ensure numeric columns
  value_cols <- setdiff(names(dt), "date")
  dt[, (value_cols) := lapply(.SD, as.numeric), .SDcols = value_cols]
  
  setkey(dt, date)
  return(dt)
}

FX_dt <- clean_bloomberg_dt(FX)
EMBI_dt <- fread(EMBI)
CDS_dt <- clean_bloomberg_dt(CDS)
Commod_dt <- clean_bloomberg_dt(Commod)
Indices_dt <- clean_bloomberg_dt(Indices)

VIX_dt   <- clean_FRED_dt(VIX)
DGS10_dt <- clean_FRED_dt(DGS10)
DGS2_dt <- clean_FRED_dt(DGS2)
Brent_dt <- clean_FRED_dt(Brent)
USDIndex_dt <- clean_FRED_dt(USDIndex)

FX_dates <- FX_dt[, .(date)]


EMBI_dt        <- EMBI_dt[FX_dates, on = "date"]
DGS10_dt       <- DGS10_dt[FX_dates, on = "date"]
DGS2_dt        <- DGS2_dt[FX_dates, on = "date"]
CDS_dt        <- CDS_dt[FX_dates, on = "date"]
VIX_dt         <- VIX_dt[FX_dates, on = "date"]
Indices_dt     <- Indices_dt[FX_dates, on = "date"]
Commod_dt <- Commod_dt[FX_dates, on = "date"]
Brent_dt       <- Brent_dt[FX_dates, on = "date"]
USDIndex_dt   <- USDIndex_dt[FX_dates, on = "date"]

final_dt <- FX_dt
final_dt <- merge(final_dt, EMBI_dt, all = TRUE)
final_dt <- merge(final_dt, DGS10_dt, all = TRUE)
final_dt <- merge(final_dt, DGS2_dt, all = TRUE)
final_dt <- merge(final_dt, VIX_dt, all = TRUE)
final_dt <- merge(final_dt, Indices_dt, all = TRUE)
final_dt <- merge(final_dt, Commod_dt, all = TRUE)
final_dt <- merge(final_dt, Brent_dt, all = TRUE)
final_dt <- merge(final_dt, CDS_dt, all = TRUE)
final_dt <- merge(final_dt, USDIndex_dt, all = TRUE)

rm(list=c(EMBI_dt,DGS10_dt, DGS2_dt, CDS_dt, VIX_dt, Indices_dt, Commod_dt, Brent_dt, USDIndex_dt))
