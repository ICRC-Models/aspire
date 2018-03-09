## Estimate age-specific viral load distribution among men. Huerga, et al., 2016 present data from a population-based survey employing cluster probability sampling in KZN, South Africa in 2013. Individuals ages 15-59 were eligible for the study. Viral load distributions are presented by gender and age separately. We assume the ratio of female-to-male engage in care is the same across age categories to infer the male-specific age-distribution of viral load.

setwd("~/Dropbox/UW Epi Program/RB/ASPIRE")

## Specify viral load categories in Herga
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

col_func <- colorRampPalette(c("red", "blue"))
cols <- col_func(ncol(m_vl_dist))

pdf(file = paste0(getwd(), "/pres-mtn-2017-12-21/vl-dist.pdf"), height = 5, width = 8)
barplot(height    = t(m_vl_dist),
        beside    = T,
        names.arg = c("<100", "100-999", "1,000-9,999", "10,000-99,999", ">100,000"),
        main      = "Viral load distribution among HIV-positive male partners",
        xlab      = "Viral load (copies/mL)",
        ylab      = "Proportion",
        col       = cols)

legend("topright", legend = c("15-24", "25-34", "35-44", "45-59"), col = cols, fill = cols, bty = "n")

dev.off()
