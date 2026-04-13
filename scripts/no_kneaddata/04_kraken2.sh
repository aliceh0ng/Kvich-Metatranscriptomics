#!/bin/bash
#SBATCH --job-name=kraken2_conf0.05
#SBATCH --time=24:00:00
#SBATCH --account=st-ctropini-1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --array=0-118
#SBATCH --output=output/logs/kraken2_conf0.05_%A_%a.out
#SBATCH --error=output/logs/kraken2_conf0.05_%A_%a.err
#SBATCH --mail-user=alice.hong@ubc.ca
#SBATCH --mail-type=ALL

###############################################################################

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate kraken_env

set -euo pipefail

PROJ=/scratch/st-ctropini-1/ahong/bfrag3
DB=/scratch/st-ctropini-1/ahong/kraken2_standard_db
INPUT=$PROJ/input
OUTDIR=$PROJ/output/kraken2_conf0.05
SAMPLES=$PROJ/samples.txt

mkdir -p "$OUTDIR"

SAMPLE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$SAMPLES")

R1=$INPUT/${SAMPLE}_unmapped_R1.fq
R2=$INPUT/${SAMPLE}_unmapped_R2.fq

echo "[$(date)] Running Kraken2 (--confidence 0.05) on $SAMPLE"

kraken2 \
    --db "$DB" \
    --threads "$SLURM_CPUS_PER_TASK" \
    --paired \
    --confidence 0.05 \
    --report "$OUTDIR/${SAMPLE}.report" \
    --output /dev/null \
    "$R1" "$R2"

echo "[$(date)] Done: $OUTDIR/${SAMPLE}.report"
