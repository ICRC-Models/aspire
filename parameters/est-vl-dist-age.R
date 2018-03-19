## Purpose: Estimate age-specific viral load distribution among men. Huerga, et al., 2016 present data from a population-based survey employing cluster probability sampling in KZN, South Africa in 2013. Individuals ages 15-59 were eligible for the study. Viral load distributions are presented by gender and age separately. We assume the ratio of female-to-male engage in care is the same across age categories to infer the male-specific age-distribution of viral load.
## Date:    15 December 2017
## Author:  Kathryn Peebles

setwd("~/Documents/code/aspire")

## Specify viral load categories in Huerga
vl_cats <- c("<100", "100-999", "1000-9999", "10000-99999", "100000-999999")

## Viral load distribution by gender
g_vl <- matrix(data = c(551, 87, 154, 171, 101, # women
                        139, 19, 43, 70, 60),   # men
               nrow = length(vl_cats),
               dimnames = list(vl_cats, c("f", "m")))

## Viral load distribution by age
age_vl <- list("15-24" = c(66, 20, 52, 66, 38),
               "25-34" = c(202, 44, 77, 96, 68),
               "35-44" = c(235, 26, 40, 50, 30),
               "45-59" = c(187, 16, 28, 29, 25))

## Estimate proportion of each viral load category that is constituted by men
prop_m <- sapply(1:nrow(g_vl), function(x) { g_vl[x, "m"]/sum(g_vl[x, ]) })

## Assuming the proportion of men in each viral load category applies to all age groups, estimate the number of men in each age and viral load category.
m_vl <- lapply(names(age_vl), function(x) { age_vl[[x]] * prop_m })
m_vl_dt <- do.call(cbind, m_vl)
colnames(m_vl_dt) <- names(age_vl)
rownames(m_vl_dt) <- vl_cats

## Estimate proportion of males in each viral load category by age
m_vl_dist <- sapply(1:ncol(m_vl_dt), function(x) {
  sapply(1:nrow(m_vl_dt), function(y) { round(m_vl_dt[y, x]/sum(m_vl_dt[, x]), 2) })
})
colnames(m_vl_dist) <- names(age_vl)
rownames(m_vl_dist) <- vl_cats

## Check that column proportions sum to 1
sapply(1:ncol(m_vl_dist), function(x) sum(m_vl_dist[, x]))

## Age categories in Huerga do not match male categories in age-mixing matrix. Create new matrix with columns "min" and "max" and age categories corresponding to those used in model.
vl_dist <- cbind(c(1,  100, 1000, 10000, 100000), # minimum vl
                 c(99, 999, 9999, 99999, 999999), # maximum vl
                 m_vl_dist[, "15-24"], # age group 20-24
                 m_vl_dist[, "25-34"], # age group 25-29
                 m_vl_dist[, "25-34"], # age group 30-34
                 m_vl_dist[, "35-44"], # age group 35-39
                 m_vl_dist[, "35-44"], # age group 40-44,
                 m_vl_dist[, "45-59"], # age group 45-49,
                 m_vl_dist[, "45-59"]) # age group 50+
colnames(vl_dist) <- c("min", "max", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49", "50+")
rownames(vl_dist) <- NULL

## Save viral load distribution matrix
save(vl_dist, file = paste0(getwd(), "/data-public/vl_dist.RDATA"))
