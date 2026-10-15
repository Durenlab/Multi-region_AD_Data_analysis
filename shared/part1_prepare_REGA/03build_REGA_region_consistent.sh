#!/bin/bash
#SBATCH --job-name=input1
#SBATCH --cpus-per-task=4
#SBATCH --mem=500G
#SBATCH --time=48:00:00
#SBATCH --output=logs/input1.out
#SBATCH --error=logs/input1.err

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"
export REGA_SHARED_ROOT="${REGA_SHARED_ROOT:-$(dirname "${SCRIPT_DIR}")}"

ml anaconda3
source activate bench_r

Rscript 03build_REGA_region_consistent.r