#!/bin/bash
#SBATCH --job-name=flagstat
#SBATCH --time=1:00:00
#SBATCH --account=st-ctropini-1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --array=0-118
#SBATCH --output=output/logs/flagstat_%A_%a.out
#SBATCH --error=output/logs/flagstat_%A_%a.err
#SBATCH --mail-user=alice.hong@ubc.ca
#SBATCH --mail-type=ALL

###############################################################################

set -euo pipefail

PROJ=/scratch/st-ctropini-1/ahong/bfrag3
BAMDIR=$PROJ/output/bam
OUTDIR=$PROJ/output/qc/flagstat
SAMPLES=$PROJ/samples.txt

mkdir -p "$OUTDIR"

SAMPLE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$SAMPLES")

echo "[$(date)] flagstat $SAMPLE"

samtools flagstat \
    -@ "$SLURM_CPUS_PER_TASK" \
    "$BAMDIR/${SAMPLE}.bam" \
    > "$OUTDIR/${SAMPLE}.flagstat"

echo "[$(date)] Done"
