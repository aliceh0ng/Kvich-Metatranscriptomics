#!/bin/bash
#SBATCH --job-name=bracken_conf0.05
#SBATCH --time=4:00:00
#SBATCH --account=st-ctropini-1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --array=0-118
#SBATCH --output=output/logs/bracken_conf0.05_%A_%a.out
#SBATCH --error=output/logs/bracken_conf0.05_%A_%a.err
#SBATCH --mail-user=alice.hong@ubc.ca
#SBATCH --mail-type=ALL

###############################################################################

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate kraken_env

set -euo pipefail

PROJ=/scratch/st-ctropini-1/ahong/bfrag3
DB=/scratch/st-ctropini-1/ahong/kraken2_standard_db
KRAKEN_DIR=$PROJ/output/kraken2_conf0.05
OUTDIR=$PROJ/output/bracken_conf0.05
SAMPLES=$PROJ/samples.txt
READ_LEN=150    # nearest pre-built Bracken length to 151 bp reads
LEVEL=S         # species-level abundance estimation

mkdir -p "$OUTDIR"

SAMPLE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$SAMPLES")

echo "[$(date)] Running Bracken (conf0.05 input) on $SAMPLE"

bracken \
    -d "$DB" \
    -i "$KRAKEN_DIR/${SAMPLE}.report" \
    -o "$OUTDIR/${SAMPLE}.bracken" \
    -w "$OUTDIR/${SAMPLE}_bracken.report" \
    -r "$READ_LEN" \
    -l "$LEVEL"

echo "[$(date)] Done: $OUTDIR/${SAMPLE}.bracken"
