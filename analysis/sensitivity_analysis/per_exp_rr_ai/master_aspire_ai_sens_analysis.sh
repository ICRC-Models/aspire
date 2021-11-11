#!/bin/bash
for i in {1..70}; do
	export SIMNO=$i
	sbatch runsim_aspire_ai_sens_analysis.sh
done
