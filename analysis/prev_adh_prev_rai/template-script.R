
# Source functions
source(paste0(getwd(), "/set-up.R"))

# Run set of 300 simulations
results_list <- lapply(X   = 1:300, FUN = function(x) { run_sim(lambda    = params$lambda,
                                                                cond_rr   = params$cond_rr,
                                                                c         = params$c,
                                                                s         = params$s,
                                                                rr_ai     = params$rr_ai,
                                                                p_rate_rr = params$p_rate_rr,
                                                                base_male_hiv_incidence = params$base_male_hiv_incidence,
                                                                prev_ai   = sim_prev_ai,
                                                                prop_ai   = 0.078,
                                                                prop_adh  = sim_prop_adh,
                                                                rr_ring   = sim_rr_ring,
                                                                re_randomize  = F,
                                                                reduce_output = T) })

save(results_list, file = paste0("/gscratch/csde/kpeebles/results_ai_", sim_prev_ai, "_adh_", sim_prop_adh, "_rr_ring_", sim_rr_ring, ".RData"))