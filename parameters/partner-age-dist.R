## Purpose: Derive age distribution of male partners from among participants in the VOICE trial.
## Author:  Kathryn Peebles
## Date:    5 April 2018

setwd("~/Documents/code/aspire")

library(data.table)
library(ggplot2)

## Load data from VOICE trial (from baseline)
dt <- as.data.table(read.csv(file = paste0(getwd(), "/data-private/risk_factors_base_clif.csv")))

## Keep participant and partner ages. Remove observations where either variable is missing
dt <- dt[!is.na(age) & !is.na(DEMagesp), .(age, DEMagesp)]

ggplot(data = dt, aes(x = age, y = DEMagesp)) + geom_jitter(width = 0.2, size = 0.75) + geom_smooth(method = "lm") + labs(x = "Female age", y = "Male partner age", title = "Age distribution in VOICE trial at baseline") + theme(plot.title = element_text(hjust = 0.5))

age_lm <- lm(DEMagesp ~ age, data = dt)

age_cats <- c("15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49", "50+")

dt[, f_age_cat := c("15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49", "50+")[findInterval(x = age, vec = c(15, 20, 25, 30, 35, 40, 45, 50))]]

dt[, m_age_cat := c("15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49", "50+")[findInterval(x = DEMagesp, vec = c(15, 20, 25, 30, 35, 40, 45, 50))]]

age_mix_mat <- sapply(dt[, sort(unique(f_age_cat))], function(x) {
  sapply(dt[, sort(unique(m_age_cat))], function(y) {
    dt[f_age_cat == x, sum(m_age_cat == y)/nrow(dt[f_age_cat == x])]
  })
})

## Max female age in available data from VOICE is 40. Assume that age distribution of male partners for those ages 40-44 is the same as for those ages 45-49
age_mix_mat <- cbind(age_mix_mat, age_mix_mat[, 6])
colnames(age_mix_mat)  <- age_cats[1:7]

## Check that probabilities of columns sum to one
apply(X = age_mix_mat, MARGIN = 2, FUN = sum)

save(age_mix_mat, file = paste0(getwd(), "/data-private/age_mix_mat.RDATA"))
