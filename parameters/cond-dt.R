## Purpose: Save table of condom probability by country. Data are obtained from ASPIRE summaries in email from Elizabeth Brown dated 2017-11-17. Mean vaalues are proportion of women reporting condom use at last act.
## Author:  Kathryn Peebles
## Date:    19 March 2018

library(data.table)

load("~/Documents/code/aspire/data-private/forKathrn.RData")
dt <- as.data.table(d1)

cond_dt <- data.table(country   = c("mal", "sa", "uga", "zim"),
                      cond_prob = c(0.384, 0.672, 0.324, 0.527),
                      n         = dt[, sum(N), by = Country])

## Condom probability for each participant will be drawn from a beta distribution with mean equal to observed country-specific proportion condom use. Variance estimated as variance of binomial sampling distribution.
## TO DO: Ask EB: Using variance of binomial sampling distribution results in very narrow beta distributions, with consequence that no women will have very high or very low probability of condom use. Is that reasonable?
cond_dt[, var := (cond_prob * (1 - cond_prob))/n.V1]

## Estimate alpha and beta parameters of beta distribution given mean and variance
est_beta_params <- function(mu, var) {
  alpha <- ((1 - mu)/var - 1/mu) * mu ^ 2
  beta <- alpha * (1/mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

cond_dt[, c("alpha", "beta") := est_beta_params(mu = cond_prob, var = var)]

## Plot beta distribution of condom probability per country
sapply(1:nrow(cond_dt), function(x) { 
  cond_dt[x, plot(dbeta(seq(0, 1, 0.001), shape1 = alpha, shape2 = beta), type = "l")]
})

cond_dt[, `:=`(n.Country = NULL, n.V1 = NULL)]

save(cond_dt, file = "~/Documents/code/aspire/data-public/cond_dt.RData")
