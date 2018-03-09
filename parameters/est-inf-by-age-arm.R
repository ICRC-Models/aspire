library(data.table)

load("~/Dropbox/UW Epi Program/RB/Effect of AI on ring efficacy/ai-ring/data/forKathrn.RData")
dt <- as.data.table(d1)

dt[, age_cat := c("18-21", "22-26", "27-45")[findInterval(x = dt[, Age], vec = c(18, 22, 27))]] # age categories (18-21, 22-26, 27-45)

dt[, arm := ifelse(test = arm == "Placebo",
                   yes  = 0,
                   no   = 1)]

inf_obs <- as.data.table(dt %>% group_by(arm, age_cat) %>% summarise(hiv_inf = sum(Infections)))

save(inf_obs, file = "~/Dropbox/UW Epi Program/RB/Effect of AI on ring efficacy/ai-ring/data/inf_obs.RData")
