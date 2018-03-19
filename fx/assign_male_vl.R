#######################################################################################
# 
# Author: Kathryn Peebles
# Date:   15 March 2018
# assign_male_vl: Function to assign viral load to male partners given male age and country.
#
# input:  Vector of male 0/1 values indicating HIV status, vector of male age categories
# output: Vector of log10-transformed male viral load values.
# 
#######################################################################################

assign_male_vl <- function(hiv, age) {
  
  vl <- ifelse(test = hiv == 0,
               yes  = 0,
               no   = log10(unlist(sapply(unique(age), function(x) {
                 sapply(sample(x       = 1:5,
                               size    = sum(hiv[age == x] == 1),
                               replace = T,
                               prob    = vl_dist[, x]),
                        function(y) {
                          runif(n   = 1, 
                                            min = vl_dist[y, "min"],
                                            max = vl_dist[y, "max"])})}))))
  
  if(any(hiv == 0 & vl > 0)) { browser() }
  
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

  return(vl)
}