## Save table of condom probability by country. Per email from Elizabeth Brown dated 2017-11-17, values are proportion of women reporting condom use at last act.

library(data.table)

cond_dt <- data.table(country   = c("mal", "sa", "uga", "zim"),
                      cond_prob = c(0.384, 0.672, 0.324, 0.527))

save(cond_dt, file = "~/Dropbox/UW Epi Program/RB/Effect of AI on ring efficacy/ai-ring/data/cond_dt.RData")
