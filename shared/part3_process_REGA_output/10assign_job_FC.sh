#!/bin/bash
#SBATCH --job-name=assignnV
#SBATCH --output=logs/assignRE_%j_nV.out
#SBATCH --error=logs/assignRE_%j_nV.err
#SBATCH --time=72:00:00
#SBATCH --mem=600G
#SBATCH --cpus-per-task=4

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"
export REGA_SHARED_ROOT="${REGA_SHARED_ROOT:-$(dirname "${SCRIPT_DIR}")}"

mkdir -p logs
echo "============================================================"
echo "Job started"
echo "Date: $(date)"
echo "Host: $(hostname)"
echo "Working directory: $(pwd)"
echo "============================================================"

# ============================================================
# Environment
# ============================================================

ml anaconda3
source activate bench_r

Rscript 10-1assign_gene_FC.r
Rscript 10-2assign_RE_FC.r
echo "============================================================"
echo "Job finished"
echo "Date: $(date)"
echo "============================================================"