## Purpose: Call intermediate ABC processing functions for iteration t and output data and scripts for iteration t + 1.
## Author:  Kathryn Peebles
## Date:    28 January 2020

source("~/Documents/code/aspire/abc/abc_samples_process.R")
# system.time({ abc_samples_process(t = 0, N = 100000, alpha = 0.6, p_acc_min = 0.10, n_nodes = 28) }) # user: 2073.784, system: 31.224, elapsed: 2139.455

# system.time({ abc_samples_process(t = 1, N = 100000, alpha = 0.6, p_acc_min = 0.10, n_nodes = 42) }) # user: 13381.479, system: 49.089, elapsed: 13480.970

# system.time({ abc_samples_process(t = 2, N = 100000, alpha = 0.6, p_acc_min = 0.10, n_nodes = 42) }) # user: 13267.110, system: 41.764, elapsed: 13340.752

# system.time({ abc_samples_process(t = 3, N = 100000, alpha = 0.6, p_acc_min = 0.10, n_nodes = 42) }) # user: 13215.743, system: 36.518, elapsed: 26842.477

# system.time({ abc_samples_process(t = 4, N = 100000, alpha = 0.6, p_acc_min = 0.10, n_nodes = 42) }) # user: 13488.428, system: 44.985, elapsed: 13582.059

# system.time({ abc_samples_process(t = 5, N = 100000, alpha = 0.6, p_acc_min = 0.10, n_nodes = 42) }) # user: 13402.396, system: 42.242, elapsed: 15206.369

# system.time({ abc_samples_process(t = 6, N = 100000, alpha = 0.6, p_acc_min = 0.10, n_nodes = 42) }) # user: 13753.949, system: 73.761, elapsed: 13911.395

# system.time({ abc_samples_process(t = 7, N = 100000, alpha = 0.6, p_acc_min = 0.10, n_nodes = 42) }) # user: 13432.61, system:  53.86, elapsed: 13541.15

# system.time({ abc_samples_process(t = 8, N = 100000, alpha = 0.6, p_acc_min = 0.10, n_nodes = 42) }) # user: 13025.368, system: 60.849, elapsed: 13375.711

# system.time({ abc_samples_process(t = 9, N = 100000, alpha = 0.6, p_acc_min = 0.10, n_nodes = 42) }) # user: 13115.422, system: 34.037, elapsed: 13171.424

# system.time({ abc_samples_process(t = 10, N = 100000, alpha = 0.6, p_acc_min = 0.10, n_nodes = 42) }) # user: 13443.894, system: 49.768, elapsed: 13548.390

# system.time({ abc_samples_process(t = 11, N = 100000, alpha = 0.6, p_acc_min = 0.10, n_nodes = 42) }) # user: 13609.458, system: 44.755, elapsed: 13695.936

# system.time({ abc_samples_process(t = 12, N = 100000, alpha = 0.6, p_acc_min = 0.10, n_nodes = 42) }) # user: 13215.162, system: 39.117, elapsed: 14093.067

# system.time({ abc_samples_process(t = 13, N = 100000, alpha = 0.6, p_acc_min = 0.10, n_nodes = 42) }) # user: 11640.658, system: 38.067, elapsed: 11693.687
# Proportion accepted particles less than minimum acceptance criterion. Approximate convergence achieved - yay!
