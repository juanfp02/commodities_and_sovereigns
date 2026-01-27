############################################################
# 00_extract_masterfile_to_rds.R
# Read data/data_masterfile.xlsx and save one .rds per sheet
############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(lubridate)
})

# -------------------------
# USER PATHS
# -------------------------
MASTERFILE <- "data/data_masterfile.xlsx"
OUT_DIR    <- "data/processed/masterfile_rds"

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# -------------------------
# Helpers
# -------------------------
clean_colnames <- function(nms) {
  nms <- gsub("\\s+", "_", nms)
  nms <- gsub("[^A-Za-z0-9_]", "", nms)
  nms
}

parse_excel_date <- function(x) {
  # Handles: Date, POSIXct, Excel serial numeric, and character formats.
  if (inherits(x, "Date")) return(x)
  if (inherits(x, "POSIXct") || inherits(x, "POSIXt")) return(as.Date(x))
  
  if (is.numeric(x)) {
    # Excel serial dates (Windows): origin 1899-12-30 (R convention)
    return(as.Date(x, origin = "1899-12-30"))
  }
  
  x <- trimws(as.character(x))
  x[x == ""] <- NA_character_
  
  d <- suppressWarnings(as.Date(x, format = "%d.%m.%Y"))
  if (all(is.na(d))) d <- suppressWarnings(as.Date(x, format = "%Y-%m-%d"))
  if (all(is.na(d))) d <- suppressWarnings(as.Date(x, format = "%d/%m/%Y"))
  if (all(is.na(d))) d <- suppressWarnings(as.Date(x, format = "%m/%d/%Y"))
  
  if (all(is.na(d))) {
    d <- suppressWarnings(as_date(parse_date_time(x, orders = c("dmy", "ymd", "mdy"))))
  }
  d
}

to_numeric_safe <- function(x) {
  # Keep numeric as-is; for character, handle US/EU + scientific notation safely.
  if (is.numeric(x)) return(x)
  
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "N/A", "NULL", "#N/A", "#N/A N/A", "Requesting Data...", "#NAME?")] <- NA_character_
  if (all(is.na(x))) return(as.numeric(x))
  
  out <- rep(NA_real_, length(x))
  
  # scientific notation -> don't delete dots
  is_sci <- grepl("[eE][+-]?[0-9]+", x)
  xs <- x[is_sci]
  xs <- gsub(" ", "", xs)
  xs <- gsub(",", "", xs)
  out[is_sci] <- suppressWarnings(as.numeric(xs))
  
  xn <- x[!is_sci]
  if (length(xn)) {
    has_comma <- grepl(",", xn)
    has_dot   <- grepl("\\.", xn)
    
    # EU heuristic: comma exists and is likely decimal separator
    is_eu <- has_comma & ( (!has_dot) | (regexpr(",", xn) > regexpr("\\.", xn)) )
    
    # EU: remove thousands dots; comma -> dot
    xe <- xn[is_eu]
    xe <- gsub("\\.", "", xe)
    xe <- gsub(",", ".", xe)
    out[!is_sci][is_eu] <- suppressWarnings(as.numeric(xe))
    
    # US/plain: remove thousands commas only
    xu <- xn[!is_eu]
    xu <- gsub(",", "", xu)
    out[!is_sci][!is_eu] <- suppressWarnings(as.numeric(xu))
  }
  
  out
}

read_one_sheet <- function(path, sheet) {
  raw <- read_excel(path, sheet = sheet, .name_repair = "unique")
  dt <- as.data.table(raw)
  
  # Standardize date col name
  if ("Dates" %in% names(dt)) setnames(dt, "Dates", "date")
  if ("observation_date" %in% names(dt)) setnames(dt, "observation_date", "date")
  if (!"date" %in% names(dt)) setnames(dt, names(dt)[1], "date")
  
  # Parse date robustly
  dt[, date := parse_excel_date(date)]
  dt <- dt[!is.na(date)]
  setorder(dt, date)
  setkey(dt, date)
  
  # Clean column names (except date)
  value_cols <- setdiff(names(dt), "date")
  if (length(value_cols)) {
    new_names <- clean_colnames(value_cols)
    setnames(dt, value_cols, new_names)
  }
  
  # Convert non-date columns to numeric safely
  value_cols <- setdiff(names(dt), "date")
  if (length(value_cols)) {
    dt[, (value_cols) := lapply(.SD, to_numeric_safe), .SDcols = value_cols]
  }
  
  dt
}

# -------------------------
# Main
# -------------------------
sheets <- excel_sheets(MASTERFILE)
message("Found sheets: ", paste(sheets, collapse = ", "))

dt_list <- list()
manifest <- data.table(sheet = character(), file = character(), n = integer(),
                       date_min = as.Date(character()), date_max = as.Date(character()))

for (s in sheets) {
  message("Reading sheet: ", s)
  dt <- read_one_sheet(MASTERFILE, s)
  
  # Save
  s_clean <- make.names(s)
  out_file <- file.path(OUT_DIR, paste0(s_clean, ".rds"))
  saveRDS(dt, out_file)
  
  dt_list[[s_clean]] <- dt
  manifest <- rbind(
    manifest,
    data.table(
      sheet = s,
      file = out_file,
      n = nrow(dt),
      date_min = if (nrow(dt)) min(dt$date) else as.Date(NA),
      date_max = if (nrow(dt)) max(dt$date) else as.Date(NA)
    )
  )
}


print(manifest)
message("DONE. RDS saved to: ", normalizePath(OUT_DIR))

