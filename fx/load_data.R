#######################################################################################
# 
# Author: Kathryn Peebles
# Date:   9 March 2018
# load_data: function to read and load parameters and data that will be used in ring model.
#
#######################################################################################

load_data <- function(working_directory = wd) {
  
  male_prev_mal <- t(read_excel(paste0(wd, "/data-public/male-hiv-prevalence-mal.xlsx"),
                                 col_names = T))
  
  male_prev_sa  <- t(read_excel(paste0(wd, "/data-public/male-hiv-prevalence-sa.xlsx"),
                                col_names = T))
  
  male_prev_uga <- t(read_excel(paste0(wd, "/data-public/male-hiv-prevalence-uga.xlsx"),
                                col_names = T))
  
  male_prev_zim <- t(read_excel(paste0(wd, "/data-public/male-hiv-prevalence-zim.xlsx"),
                                 col_names = T))
  
  male_prev     <<- data.frame(mal = male_prev_mal,
                               sa  = round(male_prev_sa, 3),
                               uga = round(male_prev_uga, 3),
                               zim = male_prev_zim)

  load(file = paste0(wd, "/data-public/vl_dist.RDATA"),      envir = .GlobalEnv)
  load(file = paste0(wd, "/data-public/inf_obs.RData"),      envir = .GlobalEnv)
  load(file = paste0(wd, "/data-public/inf_obs_arm.RData"),  envir = .GlobalEnv)
  load(file = paste0(wd, "/data-private/f_dt.RData"),        envir = .GlobalEnv)
  load(file = paste0(wd, "/data-public/p_rates_dt.RDATA"),   envir = .GlobalEnv)
  load(file = paste0(wd, "/data-private/age_mix_mat.RDATA"), envir = .GlobalEnv)
  load(file = paste0(wd, "/data-public/cond_dt.RData"),      envir = .GlobalEnv)
  load(file = paste0(wd, "/data-public/ai_dt.RData"),        envir = .GlobalEnv)

}
