#######################################################################################
# 
# Author: Kathryn Peebles
# Date:   15 March 2018
# assign_male_vl_trt_known and assign_male_vl_trt_unknown: Functions to assign viral load to male partners given male age and country. First function assigns viral load value to partners whose treatment status is known. Second function assigns viral load value to those with partners whose treatment status is unknown.
#
# input:  Subset of m_dt
# output: Vector of log10-transformed male viral load values.
# 
#######################################################################################

assign_male_vl_trt_known <- function(dt) {
  
  # At baseline, women who reported an HIV-positive partner also reported the treatment status of their partners (if known). For this group, use a viral load distribution that excludes viral loads < 1,000. 
  vl_dist_trt_known <- vl_dist[3:5, ]
  vl_dist_trt_known <- round(sapply(3:10, function(x) sapply(1:3, function(y) { vl_dist_trt_known[y, x]/sum(vl_dist_trt_known[, x])})), 2)
  
  vl_dist_trt_known <- cbind(vl_dist[3:5, 1:2], vl_dist_trt_known)
  colnames(vl_dist_trt_known) <- colnames(vl_dist)
  rownames(vl_dist_trt_known) <- NULL
  
  dt[hiv == 0, vl := 0]
  
  invisible(sapply(c(colnames(vl_dist_trt_known)[3:10]), function(x) {
    dt[m_age == x & hiv == 1, vl := log10(sapply(sample(x = 1:3, size = nrow(dt[m_age == x & hiv == 1]), replace = T, prob = vl_dist_trt_known[, x]), function(y) {
      runif(n = 1, min = vl_dist_trt_known[y, "min"], max = vl_dist_trt_known[y, "max"])
    }))]
  }))
  
  if(nrow(dt[hiv == 0 & vl > 0]) > 0) { browser() }

  return(dt[, vl])
}


assign_male_vl_trt_unknown <- function(dt) {

  dt[hiv == 0, vl := 0]
  
  invisible(sapply(c(colnames(vl_dist)[3:10]), function(x) {
    dt[m_age == x & hiv == 1, vl := log10(sapply(sample(x = 1:5, size = nrow(dt[m_age == x & hiv == 1]), replace = T, prob = vl_dist[, x]), function(y) {
      runif(n = 1, min = vl_dist[y, "min"], max = vl_dist[y, "max"])
    }))]
  }))
  
  if(nrow(dt[hiv == 0 & vl > 0]) > 0) { browser() }
  
  # ## Check distribution of viral load values. Prevalence is low among younger and older males, so more instability in proportions from those age groups.
  # dt <- data.table(vl = vl[hiv == 1], hiv = hiv[hiv == 1], age = age[hiv == 1])
  # dt[, vl_cat := c("<100", "100-999", "1,000-9,999", "10,000-99,999", ">100,000")[findInterval(x = 10 ^ vl, vec = c(0, 100, 1000, 10000, 100000))]]
  # dt[, age_cat := c("15-24", "25-34", "35-44", "45+")[findInterval(x = substr(x = dt$age, start = 1, stop = 2), vec = c(15, 25, 35, 45))]]
  # 
  # sapply(c("15-24", "25-34", "35-44", "45+"), function(x) {
  #   sapply(c("<100", "100-999", "1,000-9,999", "10,000-99,999", ">100,000"), function(y) {
  #     dt[age_cat == x, sum(vl_cat == y)/nrow(dt[age_cat == x])]
  #   })
  # })
  # 
  # ggplot(data = dt, aes(x = vl_cat, fill = age_cat)) + geom_bar(aes(group = age_cat, y = ..prop..), position = "dodge") + scale_x_discrete(limits = c("<100", "100-999", "1,000-9,999", "10,000-99,999", ">100,000"))
  
  return(dt[, vl])
}
