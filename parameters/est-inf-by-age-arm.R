## Purpose: Create table with number of participants, number of infections, and woman-time by age group and arm
## Author:  Kathryn Peebles
## Date:    21 March 2018

library(data.table)

load("~/Documents/code/aspire/data-private/forKathrn.RData")
dt <- as.data.table(d1)

dt[, age_cat := c("18-21", "22-26", "27-45")[findInterval(x = dt[, Age], vec = c(18, 22, 27))]] # age categories (18-21, 22-26, 27-45)

dt[, arm := ifelse(test = arm == "Placebo",
                   yes  = 0,
                   no   = 1)]

setnames(dt, old = "p-y", new = "py")

inf_obs <- as.data.table(dt %>% group_by(arm, age_cat) %>% summarise(hiv_inf = sum(Infections), n = sum(N), py = sum(py)))

inf_obs[, `:=`(inc_rate = hiv_inf/py * 100,
               se       = sqrt(hiv_inf/py/py) * 100)]

## TO DO: Ask EB about person-time in hiv data set. Incidence rates estimated here do not match those provided in Baeten, 2016.

save(inf_obs, file = "~/Documents/code/aspire/data-public/inf_obs.RData")
