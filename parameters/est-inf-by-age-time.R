## Purpose: Create table with number of participants, number of infections, and personn-time by age group and arm. Also estimate incidence in placebo arm from month 15 of the trial onward.
## Author:  Kathryn Peebles
## Date:    18 July 2018

library(data.table)
library(dplyr)
library(epiR)

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

# Add incidence from month 15 onward
load("~/Documents/code/aspire/data-private/hiv.RData")
hiv <- as.data.table(hiv)

hiv[, censordt := as.Date(censordt, format = "%d%b%Y")]

hiv[, month_15_date := enrolldt + 30 * 15]

# Limit dataset to women with at least 15 months of HIV-free follow-up and in placebo arm
dt_15 <- hiv[censordt > month_15_date & arm == "Placebo"]

# Calculate incidence
inf_15 <- as.data.table(dt_15 %>% summarise(hiv_inf = sum(censor), n = sum(!is.na(censor)), py = sum(as.numeric(censordt - month_15_date)/365)))

inf_15[, `:=`(inc_rate = hiv_inf/py * 100,
              se       = sqrt(hiv_inf/py/py) * 100,
              age_cat  = "18-45")]

setcolorder(inf_15, neworder = c(6, 1:5))

# Combine data together, including a duplicate row for incidence beyond 15 months to account for applying incidence in months 15+ to two time periods
inf_obs <- rbindlist(l = list(inf_obs, inf_15, inf_15))

save(inf_obs, file = "~/Documents/code/aspire/data-public/inf_obs.RData")

# Create table of placebo arm incidence over time to use in plotting
inc_obs_dt <- hiv[arm == "Placebo", .(ptid, enrolldt, censordt, censor, age)]

# Calculate number of visits expected given enrollment date and censor date
inc_obs_dt[, censordt := as.Date(censordt, format = "%d%b%Y")]
inc_obs_dt[, n_visits := round(as.numeric(censordt - enrolldt)/30)]

# Melt dataset
inc_obs_dt_long <- inc_obs_dt[rep(1:.N, n_visits)]
inc_obs_dt_long[, visit := 1:.N, by = ptid]

inc_obs_dt_long[, hiv := ifelse(censor == 1 & visit == n_visits, 1, 0)]

three_months <- lapply(seq(3, 32, 3), function(x) x:(x + 2))

# Create data table to store results
dt <- data.table(expand.grid(visit = 1:inc_obs_dt_long[, length(three_months)], inf = NA_real_, py = NA_real_))

for(i in 1:length(three_months)) {
  temp_dt <- inc_obs_dt_long[visit %in% three_months[[i]]]
  
  # Keep only last observation for each participant
  #temp_dt <- temp_dt[, .SD[.N], by = "ptid"]
  
  dt[i, `:=`(inf = temp_dt[, sum(hiv)],
             py  = temp_dt[, sum(nrow(temp_dt) * 30)/365])]
}

inc_ci_dt <- as.data.table(epi.conf(dat = as.matrix(dt[, .SD, .SDcols = 2:3]), ctype = "inc.rate", method = "exact") * 100)

inc_obs_time_dt <- as.data.table(cbind(inc_ci_dt, dt[, .(visit)]))

save(inc_obs_time_dt, file = "~/Documents/code/aspire/data-public/inc_obs_time.RData")

ggplot(data = inc_obs_time_dt) +
  geom_point(aes(x = visit, y = est))
