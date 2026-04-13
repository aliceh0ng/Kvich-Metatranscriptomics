#!/bin/bash
#SBATCH --job-name=bwa_align
#SBATCH --time=24:00:00            
#SBATCH --account=st-ctropini-1   
#SBATCH --nodes=1                 
#SBATCH --ntasks=1                
#SBATCH --cpus-per-task=8         
#SBATCH --mem=64G               
#SBATCH --output=output/logs/align_%A_%a.out   
#SBATCH --error=output/logs/align_%A_%a.err    
#SBATCH --mail-user=alice.hong@ubc.ca    
#SBATCH --mail-type=ALL           
#SBATCH --array=0-118             

###############################################################################

set -euo pipefail

# Paths
PROJ=/scratch/st-ctropini-1/ahong/bfrag3
REF=$PROJ/refs/combined/combined_bfrag_fnuc.fna
INPUT=$PROJ/input
OUTDIR=$PROJ/output/bam
SAMPLES=$PROJ/samples.txt
THREADS=$SLURM_CPUS_PER_TASK

mkdir -p "$OUTDIR" "$PROJ/output/logs"

# Get sample name for this array task
SAMPLE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$SAMPLES")

R1=$INPUT/${SAMPLE}_unmapped_R1.fq
R2=$INPUT/${SAMPLE}_unmapped_R2.fq

echo "[$(date)] Aligning $SAMPLE"

bwa mem \
    -t "$THREADS" \
    -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:${SAMPLE}\tPL:ILLUMINA" \
    "$REF" "$R1" "$R2" \
  | samtools sort \
    -@ "$THREADS" \
    -o "$OUTDIR/${SAMPLE}.bam"

samtools index "$OUTDIR/${SAMPLE}.bam"

echo "[$(date)] Done: $OUTDIR/${SAMPLE}.bam"
