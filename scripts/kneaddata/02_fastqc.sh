#!/bin/bash
#SBATCH --job-name=fastqc_post
#SBATCH --time=12:00:00
#SBATCH --account=st-ctropini-1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --output=/scratch/st-ctropini-1/ahong/bfrag4/output/logs/fastqc_post_%j.out
#SBATCH --error=/scratch/st-ctropini-1/ahong/bfrag4/output/logs/fastqc_post_%j.err
#SBATCH --mail-user=alice.hong@ubc.ca
#SBATCH --mail-type=ALL

###############################################################################

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate bfrag4

set -euo pipefail

PROJ=/scratch/st-ctropini-1/ahong/bfrag4
INDIR=$PROJ/kneaddata_output
OUTDIR=$PROJ/output/qc/fastqc_post

mkdir -p "$OUTDIR"

echo "[$(date)] Running FastQC on KneadData paired output"

# Run on paired output files only (excludes unmatched/contaminated)
fastqc "$INDIR"/*/*_kneaddata_paired_*.fastq \
    --outdir "$OUTDIR" \
    --threads "$SLURM_CPUS_PER_TASK" \
    2>&1 | tee "$OUTDIR/fastqc_post.log"

echo "[$(date)] Done. Output in $OUTDIR"
