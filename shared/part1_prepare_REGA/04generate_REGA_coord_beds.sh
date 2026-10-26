#!/bin/bash
#SBATCH --job-name=input2
#SBATCH --cpus-per-task=4
#SBATCH --mem=500G
#SBATCH --time=48:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"
export REGA_SHARED_ROOT="${REGA_SHARED_ROOT:-$(dirname "${SCRIPT_DIR}")}"

ml anaconda3
source activate bench_r

Rscript 04generate_REGA_coord_beds.r