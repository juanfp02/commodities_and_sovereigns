############################################################
# 02_build_exposure_bcom_sectors.R
# Build annual "BCOM sector" exposure from UN Comtrade pulls
# - Uses combined_selected_products_2000_2024.csv (HS2 codes)
# - Uses combined_total_exports_2000_2024.csv (TOTAL)
# - Uses IMF GDP (annual, USD)
#
# Output:
# - data/processed/exposure/exposure_bcom_sectors_2000_2024.rds
# - data/processed/exposure/exposure_bcom_sectors_2000_2024.csv
# - data/processed/exposure/avg_exposure_bcom_sectors.csv
############################################################

suppressPackageStartupMessages({
  library(data.table)
})

#========================
# 0) PATHS (EDIT IF NEEDED)
#========================
SEL_FILE <- "data/comtrade/combined_selected_products_2000_2024.csv"
TOT_FILE <- "data/comtrade/combined_total_exports_2000_2024.csv"
GDP_FILE <- "data/imf/weo_gdp_ngdpd_annual_2000_2024.csv"

OUT_DIR <- "data/processed/exposure"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

OUT_RDS <- file.path(OUT_DIR, "exposure_bcom_sectors_2000_2024.rds")
OUT_CSV <- file.path(OUT_DIR, "exposure_bcom_sectors_2000_2024.csv")
OUT_AVG <- file.path(OUT_DIR, "avg_exposure_bcom_sectors.csv")

#========================
# 1) BCOM HS2 GROUPS (NON-OVERLAPPING)
#========================
BCOM_HS2 <- list(
  # Energy
  energy = c("27"),
  
  # Industrial metals (base metals + iron/steel)
  industrial_metals = c("72","73","74","75","76","78","79","80"),
  
  # Precious metals / gems
  precious = c("71"),
  
  # Agriculture (broad)
  agriculture = c("07","08","09","10","11","12","15","17","18","52"),
  
  # Livestock
  livestock = c("01","02")
)

#========================
# 2) LOAD COMTRADE "SELECTED PRODUCTS"
#========================
sel <- fread(SEL_FILE)

need_sel <- c("ref_year", "reporter_iso", "cmd_code", "primary_value")
miss_sel <- setdiff(need_sel, names(sel))
if (length(miss_sel)) stop("SEL_FILE missing columns: ", paste(miss_sel, collapse = ", "))

sel <- sel[
  !is.na(reporter_iso) & nchar(reporter_iso) == 3 &
    !is.na(ref_year) &
    !is.na(cmd_code)
]

# Keep clean keys and values
sel2 <- sel[, .(
  iso3 = toupper(reporter_iso),
  year = as.integer(ref_year),
  hs_code = gsub("\\D", "", as.character(cmd_code)),
  exports_usd = as.numeric(primary_value)
)]

sel2 <- sel2[!is.na(iso3) & !is.na(year) & !is.na(hs_code) & hs_code != ""]
sel2[is.na(exports_usd), exports_usd := 0]

# Convert any HS4/HS6 to HS2 (first two digits)
sel2[, hs2 := substr(hs_code, 1, 2)]
sel2 <- sel2[nchar(hs2) == 2]

# Collapse to iso3-year-hs2 (prevents any accidental duplicates)
sel_hs2 <- sel2[, .(exports_usd = sum(exports_usd, na.rm = TRUE)), by = .(iso3, year, hs2)]

#========================
# 3) BUILD SECTOR EXPORTS BY iso3-year
#========================
# Base table with all country-years appearing in selected data
base <- unique(sel_hs2[, .(iso3, year)])
setorder(base, iso3, year)

# Add one column per sector: exports_<sector>_usd
for (nm in names(BCOM_HS2)) {
  hs_set <- BCOM_HS2[[nm]]
  tmp <- sel_hs2[hs2 %in% hs_set,
                 .(x = sum(exports_usd, na.rm = TRUE)),
                 by = .(iso3, year)]
  setnames(tmp, "x", paste0("exports_", nm, "_usd"))
  base <- merge(base, tmp, by = c("iso3","year"), all.x = TRUE)
  base[is.na(get(paste0("exports_", nm, "_usd"))),
       (paste0("exports_", nm, "_usd")) := 0]
}

#========================
# 4) LOAD TOTAL EXPORTS (TOTAL)
#========================
tot <- fread(TOT_FILE)

need_tot <- c("ref_year", "reporter_iso", "primary_value")
miss_tot <- setdiff(need_tot, names(tot))
if (length(miss_tot)) stop("TOT_FILE missing columns: ", paste(miss_tot, collapse = ", "))

tot2 <- tot[, .(
  iso3 = toupper(reporter_iso),
  year = as.integer(ref_year),
  total_exports_usd = as.numeric(primary_value)
)]
tot2 <- tot2[!is.na(iso3) & nchar(iso3) == 3 & !is.na(year)]
tot2 <- tot2[, .(total_exports_usd = sum(total_exports_usd, na.rm = TRUE)), by = .(iso3, year)]

# Merge totals onto base
exp_dt <- merge(base, tot2, by = c("iso3","year"), all.x = TRUE)

#========================
# 5) LOAD GDP (IMF) AND MERGE
#========================
gdp <- fread(GDP_FILE)
setnames(gdp, tolower(names(gdp)))

need_gdp <- c("country","year","value")
miss_gdp <- setdiff(need_gdp, names(gdp))
if (length(miss_gdp)) stop("GDP_FILE missing columns: ", paste(miss_gdp, collapse = ", "))

gdp2 <- gdp[, .(
  iso3 = toupper(as.character(country)),
  year = as.integer(year),
  gdp_usd = as.numeric(value)
)]
gdp2 <- gdp2[!is.na(iso3) & nchar(iso3) == 3 & !is.na(year)]

exp_dt <- merge(exp_dt, gdp2, by = c("iso3","year"), all.x = TRUE)

#========================
# 6) COMPUTE EXPOSURES
#========================
# share_sector = exports_sector / total_exports
# expgdp_sector = exports_sector / gdp
sector_export_cols <- grep("^exports_.*_usd$", names(exp_dt), value = TRUE)

for (ec in sector_export_cols) {
  nm <- sub("^exports_", "", ec)
  nm <- sub("_usd$", "", nm)
  
  exp_dt[, paste0("share_", nm) :=
           fifelse(total_exports_usd > 0, get(ec) / total_exports_usd, NA_real_)]
  
  exp_dt[, paste0("expgdp_", nm) :=
           fifelse(gdp_usd > 0, get(ec) / gdp_usd, NA_real_)]
}

setorder(exp_dt, iso3, year)

#========================
# 7) AVERAGE EXPOSURE TABLE (2000-2024 average)
#========================
avg_cols <- grep("^(share_|expgdp_)", names(exp_dt), value = TRUE)

avg_dt <- exp_dt[, c(
  .(years_obs = sum(!is.na(total_exports_usd))),
  lapply(.SD, function(x) mean(x, na.rm = TRUE))
), by = iso3, .SDcols = avg_cols]

setorder(avg_dt, iso3)

#========================
# 8) SAVE
#========================
saveRDS(exp_dt, OUT_RDS)
fwrite(exp_dt, OUT_CSV)
fwrite(avg_dt, OUT_AVG)

cat("DONE\n")
cat("Saved annual sector exposure: ", OUT_RDS, "\n", sep = "")
cat("Saved annual sector exposure CSV: ", OUT_CSV, "\n", sep = "")
cat("Saved average exposure table: ", OUT_AVG, "\n", sep = "")
