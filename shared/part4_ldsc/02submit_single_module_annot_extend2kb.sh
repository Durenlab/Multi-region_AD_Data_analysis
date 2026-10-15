#!/bin/bash
#SBATCH --job-name=anno_2kb
#SBATCH --cpus-per-task=2
#SBATCH --mem=128GB
#SBATCH --time=24:00:00
#SBATCH -o logs/anno_2kb_%A_%a.out
#SBATCH -e logs/anno_2kb_%A_%a.err

set -euo pipefail

bash scripts/02run_single_module_annot_extend2kb_array.sh