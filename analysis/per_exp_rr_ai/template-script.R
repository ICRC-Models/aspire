
# Source functions
source(paste0(getwd(), "/set-up.R"))

# Run set of 250 simulations
results_list <- lapply(X   = 1:250, FUN = function(x) { run_sim(lambda    = params$lambda,
                                                                cond_rr   = params$cond_rr,
                                                                c         = params$c,
                                                                s         = params$s,
                                                                rr_ai     = params$rr_ai,
                                                                p_rate_rr = params$p_rate_rr,
                                                                base_male_hiv_incidence = params$base_male_hiv_incidence,
                                                                prev_ai   = "100",
                                                                prop_ai   = sim_prop_ai,
                                                                prop_adh  = 0.90,
                                                                rr_ring   = sim_rr_ring,
                                                                re_randomize  = sim_re_randomize,
                                                                reduce_output = T) })

save(results_list, file = paste0("/gscratch/csde/kpeebles/results_prop_ai_", sim_prop_ai, "_rr_ring_", sim_rr_ring, "_rerandomize_", sim_re_randomize, ".RData"))