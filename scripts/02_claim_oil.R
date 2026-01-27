suppressPackageStartupMessages({
  library(data.table)
  library(fixest)
})

panel <- readRDS("data/processed/panel_weekly.rds")

# Pick oil exporters (edit freely)
OIL_EXPORTERS <- c("Angola","Nigeria","Ecuador","Kazakhstan")

dt <- panel[country %in% OIL_EXPORTERS]

# Oil return column generated from CO1_Comdty
x_oil <- "r_CO1_Comdty"

# Model 1: simple pooled within-country with date FE (cleanest)
m1 <- feols(
  r_embi ~ get(x_oil) + r_fx | country + date,
  data = dt,
  vcov = ~country
)

# Model 2: add controls (will be collinear if date FE, but keep as alternate spec w/o date FE)
m2 <- feols(
  r_embi ~ get(x_oil) + r_fx + r_vix + d_us10y + r_spx | country,
  data = dt,
  vcov = ~country + date
)

# Country-by-country (super simple, what you asked for earlier)
run_one <- function(cty) {
  dti <- dt[country == cty]
  m <- lm(r_embi ~ get(x_oil) + r_fx + r_vix + d_us10y + r_spx, data = dti)
  co <- summary(m)$coefficients
  data.table(
    country = cty,
    beta_oil = co["get(x_oil)", "Estimate"],
    p_oil    = co["get(x_oil)", "Pr(>|t|)"],
    n = nrow(dti)
  )
}

res <- rbindlist(lapply(OIL_EXPORTERS, run_one))
res <- res[order(p_oil)]

dir.create("output/claim1", showWarnings = FALSE, recursive = TRUE)

sink("output/claim1/claim1_oil_models.txt")
cat("=== Claim 1: Oil exporters ===\n\n")
cat("\n--- m1: country + date FE ---\n"); print(summary(m1))
cat("\n--- m2: country FE (controls included) ---\n"); print(summary(m2))
cat("\n--- Per-country OLS ---\n"); print(res)
sink()

fwrite(res, "output/claim1/claim1_oil_per_country.csv")

cat("Wrote:\n")
cat(" - output/claim1/claim1_oil_models.txt\n")
cat(" - output/claim1/claim1_oil_per_country.csv\n")