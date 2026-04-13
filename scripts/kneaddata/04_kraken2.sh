#!/bin/bash
#SBATCH --job-name=kraken2
#SBATCH --time=24:00:00
#SBATCH --account=st-ctropini-1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --array=0-85
#SBATCH --output=/scratch/st-ctropini-1/ahong/bfrag4/output/logs/kraken2_%A_%a.out
#SBATCH --error=/scratch/st-ctropini-1/ahong/bfrag4/output/logs/kraken2_%A_%a.err
#SBATCH --mail-user=alice.hong@ubc.ca
#SBATCH --mail-type=ALL

###############################################################################

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate kraken_env

set -euo pipefail

PROJ=/scratch/st-ctropini-1/ahong/bfrag4
DB=/scratch/st-ctropini-1/ahong/kraken2_standard_db
INPUT=$PROJ/kneaddata_output
OUTDIR=$PROJ/output/kraken2
SAMPLES=$PROJ/samples_filtered.txt

mkdir -p "$OUTDIR"

SAMPLE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$SAMPLES")

R1=$INPUT/${SAMPLE}/${SAMPLE}_unmapped_R1_kneaddata_paired_1.fastq
R2=$INPUT/${SAMPLE}/${SAMPLE}_unmapped_R1_kneaddata_paired_2.fastq

echo "[$(date)] Running Kraken2 on $SAMPLE"

kraken2 \
    --db "$DB" \
    --threads "$SLURM_CPUS_PER_TASK" \
    --confidence 0.05 \
    --paired \
    --report "$OUTDIR/${SAMPLE}.report" \
    --output /dev/null \
    "$R1" "$R2"

echo "[$(date)] Done: $OUTDIR/${SAMPLE}.report"
