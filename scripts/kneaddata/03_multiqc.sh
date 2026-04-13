#!/bin/bash
#SBATCH --job-name=multiqc_post
#SBATCH --time=1:00:00
#SBATCH --account=st-ctropini-1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=/scratch/st-ctropini-1/ahong/bfrag4/output/logs/multiqc_post_%j.out
#SBATCH --error=/scratch/st-ctropini-1/ahong/bfrag4/output/logs/multiqc_post_%j.err
#SBATCH --mail-user=alice.hong@ubc.ca
#SBATCH --mail-type=ALL

###############################################################################

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate base

set -euo pipefail

PROJ=/scratch/st-ctropini-1/ahong/bfrag4

mkdir -p "$PROJ/output/qc/multiqc_post" "$PROJ/output/qc/multiqc_comparison"

echo "[$(date)] Running MultiQC on post-KneadData FastQC output"

# Post-KneadData report
multiqc "$PROJ/output/qc/fastqc_post" \
    --outdir "$PROJ/output/qc/multiqc_post" \
    --filename multiqc_report \
    2>&1 | tee "$PROJ/output/qc/multiqc_post/multiqc.log"

echo "[$(date)] Running combined pre/post MultiQC for comparison"

# create combined pre and post report for direct comparisons
multiqc \
    "$PROJ/output/qc/fastqc_pre" \
    "$PROJ/output/qc/fastqc_post" \
    --outdir "$PROJ/output/qc/multiqc_comparison" \
    --filename multiqc_comparison \
    2>&1 | tee "$PROJ/output/qc/multiqc_comparison/multiqc.log"

echo "[$(date)] Done."
echo "  Post-only report:    $PROJ/output/qc/multiqc_post/multiqc_report.html"
echo "  Comparison report:   $PROJ/output/qc/multiqc_comparison/multiqc_comparison.html"
