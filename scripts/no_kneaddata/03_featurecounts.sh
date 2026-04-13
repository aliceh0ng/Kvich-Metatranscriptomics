#!/bin/bash
#SBATCH --job-name=featurecounts
#SBATCH --time=24:00:00
#SBATCH --account=st-ctropini-1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=output/logs/featurecounts_%j.out
#SBATCH --error=output/logs/featurecounts_%j.err
#SBATCH --mail-user=alice.hong@ubc.ca
#SBATCH --mail-type=ALL

###############################################################################

set -euo pipefail

PROJ=/scratch/st-ctropini-1/ahong/bfrag3
GTF=$PROJ/refs/combined/combined_bfrag_fnuc.gtf
BAMDIR=$PROJ/output/bam
OUTDIR=$PROJ/output/counts

mkdir -p "$OUTDIR"

echo "[$(date)] Running featureCounts"

featureCounts \
    -p \
    -B \
    -C \
    -T "$SLURM_CPUS_PER_TASK" \
    -a "$GTF" \
    -F GTF \
    -t CDS \
    -g gene_id \
    -o "$OUTDIR/counts.txt" \
    "$BAMDIR"/*.bam

echo "[$(date)] Done: $OUTDIR/counts.txt"
