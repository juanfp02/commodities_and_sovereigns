suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(lmtest)
  library(sandwich)
  library(fixest)
  library(car)
})

source("data_ingestion.R")


expo <- out
expo[, country := toupper(country)]
expo[, year := as.integer(year)]

# Choose which exposure metric to average:
# - shares: share_<commodity>
# - exports/GDP: expgdp_<commodity>
share_cols  <- grep("^share_",  names(expo), value = TRUE)
expgdp_cols <- grep("^expgdp_", names(expo), value = TRUE)

# 1) Simple unweighted average over years (equal weight per year)
avg_share <- expo[, lapply(.SD, mean, na.rm = TRUE),
                  by = country, .SDcols = share_cols]

avg_expgdp <- expo[, lapply(.SD, mean, na.rm = TRUE),
                   by = country, .SDcols = expgdp_cols]

avg_all <- merge(avg_share, avg_expgdp, by = "country", all = TRUE)

# 2) Optional: export-weighted averages (weights = total_exports_usd)
# This is often better if export series are noisy early years.
wmean <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  sum(x[ok] * w[ok]) / sum(w[ok])
}

avg_share_w <- expo[, c(list(country = country[1]),
                        lapply(.SD, wmean, w = total_exports_usd)),
                    by = country, .SDcols = share_cols]

setnames(avg_share_w, share_cols, paste0(share_cols, "_wavg"))

# Merge weighted + unweighted
avg_all <- merge(avg_all, avg_share_w, by = "country", all = TRUE)

