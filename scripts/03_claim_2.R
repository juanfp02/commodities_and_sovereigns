# ====================================
# SIMPLE GLOBAL STATISTICS METHOD
# ====================================

suppressPackageStartupMessages(library(data.table))

CDS_RDS   <- "data/processed/masterfile_rds/CDS_daily.rds"
BRENT_CSV <- "data/DCOILBRENTEU.CSV"

# Define groups
CONTROLS <- c("Morocco_USD","Egypt_USD","Poland_USD","Brazil_USD")
OIL_EXPORTERS <- c("SaudiArabia_USD","Kazakhstan_USD","Iraq_USD","Colombia_USD","Qatar_USD")

# Threshold
THRESHOLD <- 2  # 2.5 sigma

# Load CDS
cds_wide <- readRDS(CDS_RDS)
setDT(cds_wide)
if ("Dates" %in% names(cds_wide)) setnames(cds_wide, "Dates", "date")
cds_wide[, date := as.Date(date)]

cds_cols <- setdiff(names(cds_wide), "date")
cds_long <- melt(cds_wide, id.vars = "date", measure.vars = cds_cols,
                 variable.name = "cds_col", value.name = "cds_spread")
cds_long[, cds_spread := as.numeric(cds_spread)]

# Load Brent
brent <- fread(BRENT_CSV)
setnames(brent, c("observation_date","DCOILBRENTEU"), c("date","brent_level"))
brent[, date := as.Date(date)]
brent[, brent_level := as.numeric(brent_level)]
brent <- brent[is.finite(brent_level)]
setorder(brent, date)
brent[, r_brent := log(brent_level) - log(shift(brent_level))]

# Merge
dt <- merge(cds_long[, .(date, cds_col, cds_spread)],
            brent[, .(date, r_brent)],
            by = "date", all.x = TRUE)
setorder(dt, cds_col, date)
dt[, r_cds := log(cds_spread) - log(shift(cds_spread)), by = cds_col]

# Clean
dt <- dt[is.finite(r_brent) & is.finite(r_cds)]

# ====================================
# GLOBAL STATISTICS JUMP DETECTION
# ====================================

# Brent: calculate GLOBAL mean and SD
brent_mean <- mean(dt$r_brent, na.rm = TRUE)
brent_sd <- sd(dt$r_brent, na.rm = TRUE)

cat("Brent global mean:", brent_mean, "\n")
cat("Brent global SD:", brent_sd, "\n")

# Define jumps: negative moves beyond threshold
dt[, z_brent := (r_brent - brent_mean) / brent_sd]
dt[, brent_jump := as.integer(z_brent < -THRESHOLD)]

# CDS: calculate per-country mean and SD
dt[, `:=`(
  cds_mean = mean(r_cds, na.rm = TRUE),
  cds_sd = sd(r_cds, na.rm = TRUE)
), by = cds_col]

# Define jumps: positive moves beyond threshold
dt[, z_cds := (r_cds - cds_mean) / cds_sd, by = cds_col]
dt[, cds_jump := as.integer(z_cds > THRESHOLD), by = cds_col]

# Co-jump
dt[, cojump := brent_jump & cds_jump]

# Count total jumps
cat("\nTotal Brent jump observations:", sum(dt$brent_jump, na.rm=T))
cat("\nUnique Brent jump days:", length(unique(dt[brent_jump==1]$date)), "\n")

# Check March 2020
cat("\n=== MARCH 2020 VERIFICATION ===\n")
march2020 <- dt[date >= "2020-03-01" & date <= "2020-03-31"]
march_summary <- unique(march2020[, .(date, r_brent, z_brent, brent_jump)])
print(march_summary[order(r_brent)][1:10])

# Check famous dates
cat("\n=== FAMOUS CRASH DATES ===\n")
famous_dates <- as.Date(c("2020-03-09", "2020-03-18", "2016-01-20", "2015-01-05"))
for (d in famous_dates) {
  row <- unique(dt[date == d, .(date, r_brent, z_brent, brent_jump)])
  if (nrow(row) > 0) {
    cat(sprintf("%s: r=%.4f, z=%.2f, jump=%d\n", 
                d, row$r_brent, row$z_brent, row$brent_jump))
  }
}

# ====================================
# GROUP ANALYSIS
# ====================================

dt[, grp := fcase(
  cds_col %in% OIL_EXPORTERS, "OilExporter",
  cds_col %in% CONTROLS, "Control",
  default = NA_character_
)]
dt <- dt[!is.na(grp)]

# Country-level stats
co_by_country <- dt[, {
  n_brent <- sum(brent_jump)
  n_cds <- sum(cds_jump)
  n_both <- sum(cojump)
  
  list(
    days = .N,
    n_brent_jumps = n_brent,
    n_cds_jumps = n_cds,
    n_cojumps = n_both,
    p_cds_given_brent = if(n_brent > 0) n_both / n_brent else NA_real_
  )
}, by = .(grp, cds_col)][order(grp, -p_cds_given_brent)]

print(co_by_country)

# Group-level stats
grp_stats <- dt[, {
  A <- sum(cds_jump==1 & brent_jump==1)
  B <- sum(cds_jump==1 & brent_jump==0)
  C <- sum(cds_jump==0 & brent_jump==1)
  D <- sum(cds_jump==0 & brent_jump==0)
  
  list(
    A=A, B=B, C=C, D=D,
    p_cds_given_brent = A/(A+C),
    p_cds_baseline = (A+B)/(A+B+C+D),
    odds_ratio = (A*D)/(B*C)
  )
}, by = grp]

print(grp_stats)

# Fisher tests
tab_control <- with(dt[grp=="Control"], table(cds_jump, brent_jump))
tab_export <- with(dt[grp=="OilExporter"], table(cds_jump, brent_jump))

cat("\n=== FISHER TESTS ===\n")
cat("\nControl:\n")
print(fisher.test(tab_control))
cat("\nOil Exporters:\n")
print(fisher.test(tab_export))

