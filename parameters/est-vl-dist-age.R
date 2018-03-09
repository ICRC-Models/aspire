## Analyze data from Linkages study to estimate the distribution of viral load among males by age. Data were shared in email dated 2017-12-12 from Torin Schaafsma, and are stored on ICRC home drive.

setwd("C:/Users/kpeebles/Dropbox/UW Epi Program/RB/sources-eff-dilution-ring/ring-eff")

load("H:/linkages_vl.RDATA")

library(data.table)

## Subset data from males
m_vl <- as.data.table(linkages_viral_loads[!linkages_viral_loads$female, ])

## Groups ages in same age groups as are used in aspire model
age_mix_mat <- read_excel(paste0(getwd(), "/data/age-mixing-matrix.xlsx"),
                          range = "B1:I9",
                          col_names = T)

age_cats <- names(age_mix_mat)

m_vl[, age_cat := age_cats[findInterval(x = m_vl$age, vec = c(15, 20, 25, 30, 35, 40, 45, 50))]] # age categories (15-19, 20-24, 25-29, 30-34, 35-39, 40-44, 45-49, 50+)

## Try larger age groups
age_cats_wide <- c("15-24", "25-34", "35-44", "45-59")
m_vl[, age_cat_wide := age_cats_wide[findInterval(x = m_vl$age, vec = c(15, 25, 35, 45))]]

## Create viral load categories
vl_cats <- c("1-49", "50-1500", "1501-10000", "10001-50000", ">50000")

## Baseline viral load categories
m_vl[, vl_cat_bl := vl_cats[findInterval(x = m_vl$baseline_viral_load, vec = c(0, 50, 1501, 10001, 50001))]]

## Crosstab of age category and viral load category
with(m_vl, table(age_cat, vl_cat_bl))
with(m_vl, table(age_cat_wide, vl_cat_bl))
# Crosstab shows small numbers. May be better to estimate a linear relationship between age and viral load, and modify model function to use the linear function.

## Scatterplot of age and log10 viral load
m_vl[, log_vl := log10(baseline_viral_load)]
m_vl[, plot(x = age, y = log_vl, type = "p", pch = 19)]

## For regression estimate, reassign 0 values to 10
m_vl[, log_vl := ifelse(log_vl == -Inf, 1, log_vl)]
m_vl[, lm(log_vl ~ age)]

m_vl[age < 60, lm(log_vl ~ age)]

m_vl[, lines(x = age, y = 4.72286 - 0.02239 * age, col = "red", lwd = 2)]
m_vl[age < 60, lines(x = age, y = 4.80883 - 0.02491 * age, col = "blue", lwd = 2)]

abline(h = c(log10(50), log10(1501), log10(10001), log10(50001)), col = "grey")

## Set maximum viral load value
## Look at distribution of viral loads above 50,000 to determine reasonable maximum value/distribution from which to draw vl's above 50,000
m_vl[, max(c(baseline_viral_load, exit_viral_load), na.rm = T)]
