## Purpose: Create table with number of participants, number of infections, and woman-time by age group and arm
## Author:  Kathryn Peebles
## Date:    18 July 2018

library(data.table)
library(dplyr)

load("~/Documents/code/aspire/data-private/forKathrn.RData")
dt <- as.data.table(d1)

dt[, age_cat := c("18-21", "22-26", "27-45")[findInterval(x = dt[, Age], vec = c(18, 22, 27))]] # age categories (18-21, 22-26, 27-45)

# Create two tables, one that estimates incidence by age and arm, and another that includes placebo arm data only.
# Age and arm first
setnames(dt, old = "p-y", new = "py")

inf_obs_arm <- as.data.table(dt %>% group_by(arm, age_cat) %>% summarise(hiv_inf = sum(Infections), n = sum(N), py = sum(py)))

inf_obs_arm[, `:=`(inc_rate = hiv_inf/py * 100,
                   se       = sqrt(hiv_inf/py/py) * 100)]

save(inf_obs_arm, file = "~/Documents/code/aspire/data-public/inf_obs_arm.RData")

# Estimate age-specific incidence only in the placebo arm
dt <- dt[arm == "Placebo"]

inf_obs <- as.data.table(dt %>% group_by(age_cat) %>% summarise(hiv_inf = sum(Infections), n = sum(N), py = sum(py)))

inf_obs[, `:=`(inc_rate = hiv_inf/py * 100,
               se       = sqrt(hiv_inf/py/py) * 100)]

## TO DO: Include these estimates in supplement. Review with EB to make sure they're correct.

save(inf_obs, file = "~/Documents/code/aspire/data-public/inf_obs.RData")
