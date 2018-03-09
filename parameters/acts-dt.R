## Save table with mean, sd, min, and max number of sex acts in previous week by country, per email from Elizabeth Brown dated 2017-11-17.
## TO DO: Per Elizabeth's email, it looks like the number of reported acts is total acts, not per-partner acts. However, this data will be used to estimate per-partner acts, so will overestimate total acts.

library(data.table)

est_gamma_parms <- function(mu, var) {
  beta  <- mu/var
  alpha <- mu * beta
  return(list(alpha = alpha, beta = beta))
}

acts_dt <- data.table(country = c("mal", "sa", "uga", "zim"),
                      mean    = c(2.926199, 1.723164, 2.691700, 3.910979),
                      sd      = c(3.108247, 3.112891, 4.607599, 3.664751),
                      min     = c(rep(0, 4)),
                      max     = c(28, 81, 35, 44))

sapply(acts_dt[, country], function(x) {
  acts_dt[country == x, c("alpha", "beta") := est_gamma_parms(mu  = acts_dt[country == x, mean],
                                                              var = acts_dt[country == x, sd ^ 2])]
})

# ## Gamma dist plots by country
# sapply(acts_dt[, country], function(x) {
#   plot(dgamma(x = seq(0, 20, 0.01), shape = acts_dt[country == x, alpha], rate = acts_dt[country == x, beta]), type = "l", main = x)
# })

save(acts_dt, file = "~/Dropbox/UW Epi Program/RB/Effect of AI on ring efficacy/ai-ring/data/acts_dt.RData")
