## Purpose: Call intermediate ABC processing functions for iteration t and output data and scripts for iteration t + 1.
## Author:  Kathryn Peebles
## Date:    20 October 2019

source("~/Documents/code/aspire/abc/abc_samples_process.R")
# system.time({ abc_samples_process(t = 0, N = 100000, alpha = 0.6, p_acc_min = 0.10, n_nodes = 14) }) # user: 2057.639, system: 27.475, elapsed: 2112.138

# system.time({ abc_samples_process(t = 1, N = 100000, alpha = 0.6, p_acc_min = 0.10, n_nodes = 14) }) # user: 13310.782, system: 45.025, elapsed: 13407.428

# system.time({ abc_samples_process(t = 2, N = 100000, alpha = 0.6, p_acc_min = 0.10, n_nodes = 14) }) # user: 13418.167, system: 44.947, elapsed: 13526.207

# system.time({ abc_samples_process(t = 3, N = 100000, alpha = 0.6, p_acc_min = 0.10, n_nodes = 14) }) # user: 13001.403, system: 28.989, elapsed: 13051.543

# system.time({ abc_samples_process(t = 4, N = 100000, alpha = 0.6, p_acc_min = 0.10, n_nodes = 14) }) # user: 13619.311, system: 45.094, elapsed: 15838.946

# system.time({ abc_samples_process(t = 5, N = 100000, alpha = 0.6, p_acc_min = 0.10, n_nodes = 14) }) # user: 13171.854, system: 31.021, elapsed: 13229.304

# system.time({ abc_samples_process(t = 6, N = 100000, alpha = 0.6, p_acc_min = 0.10, n_nodes = 14) }) # user: 13188.733, system: 29.716, elapsed: 13228.230

# system.time({ abc_samples_process(t = 7, N = 100000, alpha = 0.6, p_acc_min = 0.10, n_nodes = 14) }) # user: 12907.978, system: 29.608 , elapsed: 12958.432

# system.time({ abc_samples_process(t = 8, N = 100000, alpha = 0.6, p_acc_min = 0.10, n_nodes = 14) }) # user: 13752.304, system: 47.761, elapsed: 14658.491

# system.time({ abc_samples_process(t = 9, N = 100000, alpha = 0.6, p_acc_min = 0.10, n_nodes = 14) }) # user: 13094.246, system: 38.971, elapsed: 13167.839

# system.time({ abc_samples_process(t = 10, N = 100000, alpha = 0.6, p_acc_min = 0.10, n_nodes = 28) }) # user: 13505.876, system: 44.533, elapsed: 13597.481

# system.time({ abc_samples_process(t = 11, N = 100000, alpha = 0.6, p_acc_min = 0.10, n_nodes = 40) }) # user: 13211.326, system: 38.505, elapsed: 13280.816

# system.time({ abc_samples_process(t = 12, N = 100000, alpha = 0.6, p_acc_min = 0.10, n_nodes = 60) })
