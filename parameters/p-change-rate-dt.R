## Purpose: Create table of partner change rates dependent on age, as reported in McGrath, 2017. There are a lot of assumptions in using this data. The rates are reported by age for men and women in aggregate, so do not reflect female-specific partner change rates by age. Additionally, these data were collected from a cohort of HIV-positive individuals (with the primary purpose of characterizing changes in relationship persistence following ART initiation), and so may not be reflective of partner change rates among the high-risk HIV-negative women included in the ASPIRE trial.
## Date:    15 March 2018
## Author:  Kathryn Peebles

setwd("~/Documents/code/aspire/data-public")

probs <- sapply(c(0.2806, 0.1574, 0.1046, 0.0424), function(x) {
  1 - exp(-x * 1/12)
})

age_cats <- list(18:21, 22:29, 30:39, 40:45)

p_rates_dt <- data.table(f_age  = unlist(age_cats), 
                         prob  = unlist(sapply(1:length(probs), function(x) rep(probs[x], length(age_cats[[x]])))))

save(p_rates_dt, file = paste0(getwd(), "/p_rates_dt.RDATA"))
