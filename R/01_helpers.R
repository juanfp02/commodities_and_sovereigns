suppressPackageStartupMessages({
  library(data.table)
})

logret <- function(x) {
  x <- as.numeric(x)
  out <- rep(NA_real_, length(x))
  ok <- is.finite(x) & x > 0
  out[2:length(x)] <- ifelse(ok[2:length(x)] & ok[1:(length(x)-1)], diff(log(x)), NA_real_)
  out
}

diff1 <- function(x) {
  x <- as.numeric(x)
  out <- rep(NA_real_, length(x))
  out[2:length(x)] <- diff(x)
  out
}

clean_cols <- function(dt) {
  # Bloomberg sometimes includes spaces; make them consistent
  setnames(dt, names(dt), gsub("\\s+", "_", names(dt)))
  dt
}