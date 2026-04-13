#!/bin/bash
#SBATCH --job-name=kneaddata
#SBATCH --time=12:00:00
#SBATCH --account=st-ctropini-1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --array=0-118
#SBATCH --output=/scratch/st-ctropini-1/ahong/bfrag4/output/logs/kneaddata_%A_%a.out
#SBATCH --error=/scratch/st-ctropini-1/ahong/bfrag4/output/logs/kneaddata_%A_%a.err
#SBATCH --mail-user=alice.hong@ubc.ca
#SBATCH --mail-type=ALL

###############################################################################

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate bfrag4

set -euo pipefail

PROJ=/scratch/st-ctropini-1/ahong/bfrag4
DB=/arc/project/st-ctropini-1/ahong/bfrag4/dbs/human_genome/hg_39
INPUT=$PROJ/input
OUTDIR=$PROJ/kneaddata_output
SAMPLES=$PROJ/samples.txt
TRIMMOMATIC=/home/alicesh/miniconda3/envs/bfrag4/share/trimmomatic-0.40-0

SAMPLE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$SAMPLES")

R1=$INPUT/${SAMPLE}_unmapped_R1.fq
R2=$INPUT/${SAMPLE}_unmapped_R2.fq

mkdir -p "$OUTDIR/$SAMPLE"

echo "[$(date)] Running KneadData on $SAMPLE"

kneaddata \
    --input1 "$R1" \
    --input2 "$R2" \
    --reference-db "$DB" \
    --output "$OUTDIR/$SAMPLE" \
    --threads "$SLURM_CPUS_PER_TASK" \
    --trimmomatic "$TRIMMOMATIC" \
    --trimmomatic-options "SLIDINGWINDOW:4:20 MINLEN:50 LEADING:3 TRAILING:3" \
    --bypass-trf \
    --remove-intermediate-output \
    --log "$OUTDIR/$SAMPLE/${SAMPLE}.log"

echo "[$(date)] Done: $OUTDIR/$SAMPLE"

# check for when only one mate mapps to human from bowtie2 stats in KneaData's log
DISCORDANT=$(grep -oP '\d+(?= pairs aligned discordantly)' "$OUTDIR/$SAMPLE/${SAMPLE}.log" | awk '{s+=$1} END{print s+0}' || true)
echo "[$(date)] Pairs with only 1 mate mapping to human (discordant): $DISCORDANT"
