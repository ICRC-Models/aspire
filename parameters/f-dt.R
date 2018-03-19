## Purpose: Create data table of study participants with fixed characteristics.
## Date:    9 March 2018
## Author:  Kathryn Peebles

## Data from Elizabeth Brown, email dated 2018-02-08. Note that there are 2,614 observations here, while Baeten, 2016 reports 2,629 women in the trial. 15 women did not report a partner in the previous 3 months, and so are excluded from model analyses.

library(data.table)

setwd("~/Documents/code/aspire/data-private")

load(paste0(getwd(), "/hiv.RData"))
dt <- as.data.table(hiv)

f_dt <- dt[, .(site, country, arm, fu_days, nVI7base, age, basenumsp, baseSTD)]

f_dt[country == "Malawi",       country := "mal"]
f_dt[country == "South Africa", country := "sa"]
f_dt[country == "Uganda",       country := "uga"]
f_dt[country == "Zimbabwe",     country := "zim"]

f_dt[, site := tstrsplit(x = site, split = ": ", keep = 2)]

f_dt[, arm := ifelse(test = arm == "Placebo",
                     yes  = 0,
                     no   = 1)]

f_dt[, f_age_cat := c("18-21", "22-26", "27-45")[findInterval(x = f_dt$age, vec = c(18, 22, 27))]] # age categories (18-21, 22-26, 27-45)

f_dt[, sti := ifelse(test = baseSTD == T, yes = 1, no = 0)]
f_dt[, baseSTD := NULL]

setkey(x = f_dt, country, age)

f_dt[, id := 1:nrow(f_dt)]

setnames(x = f_dt, old = c("age", "nVI7base", "basenumsp"), new = c("f_age", "n_acts", "n_part"))

## Set max number of baseline partners to 6. 
f_dt[n_part > 6, n_part := 6]

## Per email from EB dated 16 Nov 2017, there were no women of 2614 with 0 partners.
f_dt[n_part == 0, n_part := 1]

## Set max number of concurrent partners in simulation to n_part + 1.
f_dt[, max_part := n_part + 1]

## Set max number of weekly acts to 35.
f_dt[n_acts > 35, n_acts := 35]

## Baeten, 2016 specifies that participants were sexually active. Therefore set minimum number of monthly acts to 1.
## TO DO: Ask EB how "sexually active" was defined for study eligibility.
f_dt[n_acts == 0, n_acts := 1]

## Model runs in 30-day increments. Scale weekly number of acts to monthly number of acts.
f_dt[, n_acts := n_acts * 30/7]

## Number of acts is total acts. Convert to per-partner number of monthly acts.
f_dt[, pp_acts := n_acts/n_part]
f_dt[, n_acts := NULL]

setcolorder(x = f_dt, neworder = c("id", "country", "site", "f_age", "f_age_cat", "n_part", "max_part", "pp_acts", "sti", "arm", "fu_days"))

save(f_dt, file = paste0(getwd(), "/f_dt.RData"))
