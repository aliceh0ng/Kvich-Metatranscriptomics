## Outputs from samples directly from Kvich et al. (no further pre-processing)

No additional human read removal applied. Reads are already the unmapped fraction from the upstream human alignment pipeline (Kvich et al.: Cutadapt → SortMeRNA → bwa-mem to GRCh38).

**Samples:** 119 (all)  
**Input:** `../input/{sample}_unmapped_R1/R2.fq`

### Pipeline summary

| Step | Script | Output |
|---|---|---|
| 1. Align to B. frag + F. nuc | `01_align_bwa.sh` | `alignment/bam/` (not copied) |
| 2. Alignment QC | `02_flagstat.sh` | `qc/flagstat/` |
| 3. Quantify gene expression | `03_featurecounts.sh` | `alignment/counts/counts.txt` |
| 3b. Parse gene annotations | `03b_make_gene_annotations.py` | `alignment/gene_annotations.tsv` |
| 4. Taxonomic classification | `04_kraken2.sh` | `community_composition/kraken2/` |
| 5. Abundance estimation | `05_bracken.sh` | `community_composition/bracken/` |
| 6. Normalize Bracken counts | `07_normalize_bracken.py` | `community_composition/normalized/` |

HUMANn3 was not analysed for this track (kneaddata track only).

---

### qc/flagstat/

**Tool:** samtools flagstat  
**Script:** `scripts/no_kneaddata/02_flagstat.sh`

119 files, one per sample: `{sample}.flagstat` — alignment summary statistics from `samtools flagstat` on the bwa-mem BAM (total reads, mapped reads, paired reads, properly paired, etc.).

---

### qc/multiqc/

**Tool:** MultiQC v1.33

| File | Description |
|---|---|
| `bwa_alignment_qc.html` | MultiQC report aggregating samtools flagstat for all 119 samples |

Per-sample FastQC on raw input reads was not run for this track (reads came pre-trimmed from the upstream human pipeline).

---

### alignment/

**Tool:** bwa-mem v0.7.19  
**Reference:** `refs/combined/combined_bfrag_fnuc.fna` (B. fragilis GCF_003019295.1 + F. nucleatum GCF_016889925.1, concatenated)  
**Scripts:** `scripts/no_kneaddata/01_align_bwa.sh`, `scripts/no_kneaddata/03_featurecounts.sh`  
**Samples:** 119

#### `counts/counts.txt`

featureCounts output matrix — gene × sample raw read counts.  
Parameters: `-p` (paired-end), `-B` (both mates must map), `-C` (discard chimeric pairs), `-t CDS`, `-g gene_id`.  
Annotation: `refs/combined/combined_bfrag_fnuc.gtf`

#### `gene_annotations.tsv`

Parsed from `refs/combined/combined_bfrag_fnuc.gtf` by `scripts/no_kneaddata/03b_make_gene_annotations.py`.  
Columns: `gene_id`, `gene` (name), `product` (functional description from CDS lines).

To regenerate (run from `bmeg524/` root):
```bash
python scripts/no_kneaddata/03b_make_gene_annotations.py \
    refs/combined/combined_bfrag_fnuc.gtf \
    > no_kneaddata/alignment/gene_annotations.tsv
```

BAMs are not copied (too large). Find them at:
```
/scratch/st-ctropini-1/ahong/bfrag3/output/bam/{sample}.bam
/scratch/st-ctropini-1/ahong/bfrag3/output/bam/{sample}.bam.bai
```

---

### community_composition/kraken2/

**Source:** bfrag3 `output/kraken2_conf0.05/`  
**Tool:** Kraken2 2.17.1  
**Script:** `scripts/no_kneaddata/04_kraken2.sh`  
**Database:** Kraken2 standard (k2_standard_20260226)  
**Parameters:** `--confidence 0.05`, `--paired`; classified reads discarded (`--output /dev/null`)

119 files: `{sample}.report` — Kraken2 report format (% reads, clade counts, rank, taxon name)

---

### community_composition/bracken/

**Source:** bfrag3 `output/bracken_conf0.05/`  
**Tool:** Bracken 3.1  
**Script:** `scripts/no_kneaddata/05_bracken.sh`  
**Parameters:** `-r 150` (read length), `-l S` (species level), `-t 10`

238 files, two per sample:

| File | Description |
|---|---|
| `{sample}.bracken` | Tab-separated: name, taxonomy_id, level, kraken reads, bracken reads, fraction |
| `{sample}_bracken.report` | Bracken-corrected Kraken2 report format |

---

### community_composition/normalized/

**Script:** `scripts/no_kneaddata/07_normalize_bracken.py`  
**Method:** Kvich et al. depth normalization

1. Read counts from raw R1 FASTQ (`wc -l / 4`)
2. Scaling factor = `min_reads_across_samples / reads_in_sample`
3. Scaled counts = `bracken new_est_reads × scaling_factor`
4. Remove taxa where `log10(scaled_count + 1) < 0.9` in all samples (~8 reads)

| File | Description | Use for |
|---|---|---|
| `counts_unscaled.tsv` | Raw Bracken `new_est_reads` | Within-sample comparisons |
| `counts_scaled.tsv` | Depth-normalized counts | Across-sample comparisons |
| `counts_filtered.tsv` | Scaled + low-abundance taxa removed | Downstream stats |
| `scaling_factors.tsv` | Per-sample read counts and scaling factors | QC / reporting |

All 119 samples included (no exclusions in this track). Run from `bmeg524/` root to reproduce:
```bash
python scripts/no_kneaddata/07_normalize_bracken.py
```
