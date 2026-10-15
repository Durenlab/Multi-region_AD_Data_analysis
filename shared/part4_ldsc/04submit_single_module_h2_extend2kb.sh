#!/bin/bash
#SBATCH --job-name=h2_2kb_Bell
#SBATCH --cpus-per-task=2
#SBATCH --mem=128GB
#SBATCH --time=24:00:00
#SBATCH -o logs/h2_2kb_%A_%a.out
#SBATCH -e logs/h2_2kb_%A_%a.err

set -euo pipefail

mkdir -p h2_result_phase3_single_module_extend2kb_FC

source activate ldsc_py2

bash scripts/04run_single_module_h2_extend2kb_array.sh