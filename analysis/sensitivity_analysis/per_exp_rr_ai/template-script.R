
# Source functions
source(paste0(getwd(), "/set-up.R"))

# Run set of 250 simulations
results_list <- lapply(X   = 1:250, FUN = function(x) { run_sim(lambda    = val_lambda,
                                                                cond_rr   = val_cond_rr,
                                                                c         = val_c,
                                                                s         = val_s,
                                                                rr_ai     = val_rr_ai,
                                                                p_rate_rr = val_p_rate_rr,
                                                                base_male_hiv_incidence = val_base_male_hiv_incidence,
                                                                prev_ai   = "100",
                                                                prop_ai   = sim_prop_ai,
                                                                prop_adh  = 0.90,
                                                                rr_ring   = sim_rr_ring,
                                                                re_randomize  = sim_re_randomize,
                                                                reduce_output = T) })

save(results_list, file = paste0("/gscratch/csde/kpeebles/results_prop_ai_", sim_prop_ai, "_rr_ring_", sim_rr_ring, "_rerandomize_", sim_re_randomize, "_", rr_ai_round, ".RData"))