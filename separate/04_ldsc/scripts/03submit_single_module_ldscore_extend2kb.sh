#!/bin/bash
#SBATCH --job-name=ld_2kb
#SBATCH --cpus-per-task=4
#SBATCH --mem=256GB
#SBATCH --time=48:00:00
#SBATCH -o logs/ld_2kb_%A_%a.out
#SBATCH -e logs/ld_2kb_%A_%a.err

set -euo pipefail

mkdir -p ldscore_phase3_single_module_extend2kb

source activate ldsc_py2

bash scripts/03run_single_module_ldscore_extend2kb_array.sh
