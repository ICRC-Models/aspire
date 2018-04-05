#######################################################################################
# 
# Author: Kathryn Peebles
# Date:   15 March 2018
# assign_male_age: Function to assign age to male partners given female partner age. Assumes that male-female partner age distribution is the same by country.
#
# input:  Vector of ids for male partners to whom age will be assigned, vector of ages of corresponding female partners.
# output: Vector of age categories for male partners.
# 
#######################################################################################

assign_male_age <- function(m_ids, f_age) {

  f_ages_cat <- colnames(age_mix_mat)[findInterval(x = f_age[m_ids], vec = c(15, 20, 25, 30, 35, 40, 45, 50))] # vec creates bins [15, 20), [20, 25), [25, 30), [30, 35), [35, 40), [40, 45), [45, 50)
  
  m_partner_age <- unname(sapply(f_ages_cat, function(age_cat) {
    sample(x = rownames(age_mix_mat), size = 1, replace = F, prob = age_mix_mat[, age_cat])
  }))
  
  # ## Check distribution of female ages by male partner age category
  # boxplot(formula = f_age[m_ids] ~ m_partner_age, ylab = "female age", xlab = "male partner age category")
  
  return(m_partner_age)
}