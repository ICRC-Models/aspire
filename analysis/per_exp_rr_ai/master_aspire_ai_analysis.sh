#!/bin/bash
for i in {1..51}; do
	export SIMNO=$i
	sbatch runsim_aspire_ai_analysis.sh
done
