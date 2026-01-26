source("data_cleaning.R")

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

final_dt

#################


cols_keep <- c(
  "date",
  "ANGOL_CDS_USD_SR_5Y_D14_Corp", #CDS
  "USDAOA_Curncy", #Fx Rate
  "USGG10YR_Index", #UST10y
  "SPX_Index",
  "CO1_Comdty",
  "VIXCLS"
)


dt_sub <- final_dt[, ..cols_keep]
dt_sub <- na.omit(dt_sub)

ret_cols <- c("ANGOL_CDS_USD_SR_5Y_D14_Corp",
              "USDAOA_Curncy",
              "USGG10YR_Index",
              "SPX_Index",
              "CO1_Comdty",
              "VIXCLS")

for (col in ret_cols) {
  dt_sub[, paste0("r_", col) := c(NA, diff(log(get(col))))]
}

dt_sub <- dt_sub[!is.na(r_ANGOL_CDS_USD_SR_5Y_D14_Corp)]

m1 <- lm(
  r_ANGOL_CDS_USD_SR_5Y_D14_Corp ~ 
    r_CO1_Comdty + r_VIXCLS + r_USGG10YR_Index + r_SPX_Index + r_USDAOA_Curncy,
  data = dt_sub
)

summary(m1)
