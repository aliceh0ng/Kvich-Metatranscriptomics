#!/bin/bash
#SBATCH --job-name=humann3
#SBATCH --time=24:00:00
#SBATCH --account=st-ctropini-1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --array=0-118
#SBATCH --output=output/logs/humann3_%A_%a.out
#SBATCH --error=output/logs/humann3_%A_%a.err
#SBATCH --mail-user=alice.hong@ubc.ca
#SBATCH --mail-type=ALL

# HUMAnN3 functional profiling using Bracken-derived taxonomic profiles.
#
# Workflow per sample:
#   1. Convert _bracken.report -> MetaPhlAn3 format (.mpa) with bracken_to_mpa.py
#   2. Concatenate R1 + R2 reads into a single FASTQ (HUMAnN3 is single-file)
#   3. Run humann --taxonomic-profile to skip internal MetaPhlAn and use our profile
#   4. Remove the concatenated reads temp file to save space
#
# Outputs (per sample, in output/humann3/{sample}/):
#   {sample}_genefamilies.tsv    — UniRef90 gene families (RPK)
#   {sample}_pathabundance.tsv   — MetaCyc pathway abundances (RPK)
#   {sample}_pathcoverage.tsv    — MetaCyc pathway coverage (0–1)
#   {sample}_humann_temp/        — intermediate files (bowtie2 SAM, diamond hits);
#                                   can be deleted after confirming outputs look good
#
# After all jobs finish, merge and normalize:
#   humann_join_tables -i output/humann3 -t genefamilies -o output/humann3_merged/all_genefamilies.tsv --file_name genefamilies
#   humann_join_tables -i output/humann3 -t pathabundance -o output/humann3_merged/all_pathabundance.tsv --file_name pathabundance
#   humann_renorm_table --input output/humann3_merged/all_genefamilies.tsv --units cpm --output output/humann3_merged/all_genefamilies_cpm.tsv
#   humann_split_stratified_table --input output/humann3_merged/all_pathabundance.tsv --output output/humann3_merged/
#
#.Credit: This file was partially generated using using Claude 3.5 Sonnet by Anthropic - AH

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate humann3

SAMPLES=($(cat samples.txt))
SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"

R1="input/${SAMPLE}_unmapped_R1.fq"
R2="input/${SAMPLE}_unmapped_R2.fq"
BRACKEN_REPORT="output/bracken/${SAMPLE}_bracken.report"

MPA_DIR="output/humann3_mpa"
OUT_DIR="output/humann3/${SAMPLE}"
READS_TMP="output/humann3_reads/${SAMPLE}_cat.fq"

mkdir -p "$MPA_DIR" "$OUT_DIR" output/humann3_reads output/logs

# ── Step 1: Convert Bracken report to MetaPhlAn3 format ──────────────────────
MPA="${MPA_DIR}/${SAMPLE}.mpa"
python scripts/bracken_to_mpa.py "$BRACKEN_REPORT" "$MPA"

# ── Step 2: Concatenate R1 + R2 ───────────────────────────────────────────────
# HUMAnN3 takes a single file; concatenating R1+R2 is the standard approach
# for paired-end data. Both mates contribute independently to alignment/search.
cat "$R1" "$R2" > "$READS_TMP"

# ── Step 3: Run HUMAnN3 ───────────────────────────────────────────────────────
humann \
    --input "$READS_TMP" \
    --taxonomic-profile "$MPA" \
    --output "$OUT_DIR" \
    --threads "$SLURM_CPUS_PER_TASK" \
    --nucleotide-database /home/alicesh/project/ahong/databases/humann_databases/chocophlan \
    --protein-database /home/alicesh/project/ahong/databases/humann_databases/uniref \
    --memory-use maximum

# ── Step 4: Cleanup concatenated reads ───────────────────────────────────────
rm -f "$READS_TMP"
