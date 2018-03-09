#######################################################################################
# 
# Kathryn Peebles
# 2017-09-17
# load_data: function to read and load parameters and data that will be used in ring model.
#
#######################################################################################

load_data <- function() {
  age_mix_mat <<- read_excel(paste0(getwd(), "/data/age-mixing-matrix.xlsx"),
                             range = "B1:I9",
                             col_names = T)
  rownames(age_mix_mat) <<- colnames(age_mix_mat)
  
  male_prev_mal <- t(read_excel(paste0(getwd(), "/data/male-hiv-prevalence-mal.xlsx"),
                                 col_names = T))
  
  male_prev_sa  <- t(read_excel(paste0(getwd(), "/data/male-hiv-prevalence-sa.xlsx"),
                                col_names = T))
  
  male_prev_uga <- t(read_excel(paste0(getwd(), "/data/male-hiv-prevalence-uga.xlsx"),
                                col_names = T))
  
  male_prev_zim <- t(read_excel(paste0(getwd(), "/data/male-hiv-prevalence-zim.xlsx"),
                                 col_names = T))
  
  male_prev     <<- data.frame(mal = male_prev_mal,
                               sa  = round(male_prev_sa, 3),
                               uga = round(male_prev_uga, 3),
                               zim = male_prev_zim)

  vl_dist <<- as.matrix(read_excel(paste0(getwd(), "/data/viral-load-distribution.xlsx"),
                                   col_names = T))
  ## TO DO: Consider if this VL dist is the best one to use. Review Pronin, 2012.
  
}
