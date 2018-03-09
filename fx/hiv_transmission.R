hiv_transmission <- function(dt, t, l) {
  ## TO DO: Add HSV2 infection
  ## TO DO: Make all of the relative risks parameters that are loaded from elsewhere.
  
  # Partner 1
  # p10: cumulative probability of not being infected by unprotected vaginal acts
  dt[!is.na(vl_p1) & vl_p1 > 0 & arm == 1, p10 := (1 - l) ^ (n_acts_vi_no_condom * exp(log(2.89) * (vl_p1 - 4.0) + log(0.96) * (age - 35) + log(0.25) * adh))]
  
  dt[!is.na(vl_p1) & vl_p1 > 0 & arm == 0, p10 := (1 - l) ^ (n_acts_vi_no_condom * exp(log(2.89) * (vl_p1 - 4.0) + log(0.96) * (age - 35)))]
  
  # p11: cumulative probability of not being infected by protected vaginal acts
  dt[!is.na(vl_p1) & vl_p1 > 0 & arm == 1, p11 := (1 - l) ^ (n_acts_vi_condom * exp(log(2.89) * (vl_p1 - 4.0) + log(0.96) * (age - 35) + log(0.25) * adh + log(0.22)))]
  
  dt[!is.na(vl_p1) & vl_p1 > 0 & arm == 0, p11 := (1 - l) ^ (n_acts_vi_condom * exp(log(2.89) * (vl_p1 - 4.0) + log(0.96) * (age - 35) + log(0.22)))]
  
  # p12: cumulative probability of not being infected by unprotected anal acts
  dt[!is.na(vl_p1) & vl_p1 > 0, p12 := (1 - l) ^ (n_acts_ai_no_condom * exp(log(2.89) * (vl_p1 - 4.0) + log(0.96) * (age - 35) + log(17.25)))]
  
  # p13: cumulative probability of not being infected by protected anal acts
  dt[!is.na(vl_p1) & vl_p1 > 0, p13 := (1 - l) ^ (n_acts_ai_condom * exp(log(2.89) * (vl_p1 - 4.0) + log(0.96) * (age - 35) + log(17.25) + log(0.22)))]
  
  dt[is.na(vl_p1) | vl_p1 == 0, 13:16 := 1]
  
  # Partner 2
  # p20: cumulative probability of not being infected by unprotected vaginal acts
  dt[!is.na(vl_p2) & vl_p2 > 0 & arm == 1, p20 := (1 - l) ^ (n_acts_vi_no_condom * exp(log(2.89) * (vl_p1 - 4.0) + log(0.96) * (age - 35) + log(0.25) * adh))]
  
  dt[!is.na(vl_p2) & vl_p2 > 0 & arm == 0, p20 := (1 - l) ^ (n_acts_vi_no_condom * exp(log(2.89) * (vl_p1 - 4.0) + log(0.96) * (age - 35)))]
  
  # p21: cumulative probability of not being infected by protected vaginal acts
  dt[!is.na(vl_p2) & vl_p2 > 0 & arm == 1, p21 := (1 - l) ^ (n_acts_vi_condom * exp(log(2.89) * (vl_p1 - 4.0) + log(0.96) * (age - 35) + log(0.25) * adh + log(0.22)))]
  
  dt[!is.na(vl_p2) & vl_p2 > 0 & arm == 0, p21 := (1 - l) ^ (n_acts_vi_condom * exp(log(2.89) * (vl_p1 - 4.0) + log(0.96) * (age - 35) + log(0.22)))]
  
  # p22: cumulative probability of not being infected by unprotected anal acts
  dt[!is.na(vl_p2) & vl_p2 > 0, p22 := (1 - l) ^ (n_acts_ai_no_condom * exp(log(2.89) * (vl_p1 - 4.0) + log(0.96) * (age - 35) + log(17.25)))]
  
  # p23: cumulative probability of not being infected by protected anal acts
  dt[!is.na(vl_p2) & vl_p2 > 0, p23 := (1 - l) ^ (n_acts_ai_condom * exp(log(2.89) * (vl_p1 - 4.0) + log(0.96) * (age - 35) + log(17.25) + log(0.22)))]
  
  dt[is.na(vl_p2) | vl_p2 == 0, 17:20 := 1]
  
  # Total monthly risk
  dt[, risk := 1 - (Reduce(`*`, .SD)), .SDcols = 13:20, by = id]
  
  #  If HIV-positive at previous timestep, assign HIV status of 1 at all subsequent time steps
  dt[hiv == 1, risk := 1]

  # Assign HIV transmission as a binomial draw with probability equal to total monthly risk
  dt[, hiv := rbinom(n = nrow(dt), size = 1, prob = risk)]
  
  return(dt$hiv)
}