#!/bin/bash
for i in {1..231}; do
	export SIMNO=$i
	sbatch runsim_aspire_prev_adh_ai.sh
done
