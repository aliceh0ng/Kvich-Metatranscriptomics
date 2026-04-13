# bmeg524 — Metatranscriptomics Re-analysis

Re-analysis of Kvich et al. unmapped bacterial RNA-seq reads using two parallel preprocessing tracks.

**Date created:** 2026-04-12

---

## Dataset

119 paired-end RNA-seq samples (unmapped reads from human alignment).

- La samples: La34–La125 (with gaps) — ~92 samples
- Lasse samples: Lasse4–Lasse33 — 30 samples (`Lasse{N}` in filenames = `La{N}` in metadata)
- **La125:** sequencing control, excluded from all analyses
- Read length: ~150 bp (Illumina A00962)

Sample lists: `samples_all.txt` (119), `samples_filtered.txt` (86), `excluded_samples.txt` (33 excluded)

---

## Two-Track Structure

### Track 1: `no_kneaddata/` — all 119 samples
No additional human read removal (reads already from upstream unmapped fraction).

bwa-mem → featureCounts (B. frag + F. nuc reference) → Kraken2 (--confidence 0.05) → Bracken → normalize

HUMANn3 not analysed for this track.

### Track 2: `kneaddata/` — 86 samples
33 excluded (Lasse4–33 except Lasse9/20/25; La50, La54, La68, La73, La77, La89).

KneadData v0.12.4 (hg39 T2T, SLIDINGWINDOW:4:20 MINLEN:50 LEADING:3 TRAILING:3) → Kraken2 (--confidence 0.05) → Bracken → normalize → HUMANn3 (84 samples; La98 + La125 excluded)

See `no_kneaddata/README.md` and `kneaddata/README.md` for per-track details.

---

## Directory Structure

```
bmeg524/
├── input/                  → symlink → bfrag4/input (238 FASTQ files)
├── refs/                   → symlink → bfrag3/refs (B. frag + F. nuc genomes + bwa index)
├── dbs/                    → symlink → bfrag4/dbs (hg39 T2T KneadData reference)
├── samples_all.txt
├── samples_filtered.txt
├── excluded_samples.txt
├── scripts/
│   ├── no_kneaddata/       — pipeline scripts for Track 1
│   └── kneaddata/          — pipeline scripts for Track 2
├── no_kneaddata/           — Track 1 outputs
└── kneaddata/              — Track 2 outputs
```

---

## Databases

| Database | Location |
|---|---|
| Kraken2 standard (k2_standard_20260226) | `/arc/project/st-ctropini-1/ahong/databases/kraken2_standard_db/` |
| HUMANn3 ChocoPhlAn | `/home/alicesh/project/ahong/databases/humann_databases/chocophlan` |
| HUMANn3 UniRef90 | `/home/alicesh/project/ahong/databases/humann_databases/uniref` |
| KneadData hg39 T2T | `dbs/human_genome/hg_39` |

## Conda Environments

| Environment | Used for |
|---|---|
| `bfrag4` | KneadData v0.12.4 |
| `kraken_env` | Kraken2, Bracken, normalization scripts |
| `humann3` | HUMANn3 v3.6 |
| `base` | MultiQC |

---

## Notes

- BAMs not copied (too large) — at `/scratch/st-ctropini-1/ahong/bfrag3/output/bam/`
- `kneaddata/kneaddata_output/` is a symlink — copy manually (`cp -rL`) if moving off scratch
- `input/`, `refs/`, `dbs/` are symlinks to large files — exclude when rsyncing
- HUMANn3 output filenames contain `_cat_` (input was concatenated R1+R2: `{sample}_cat.fq`)
- High duplication in La samples is expected for RNA-seq metatranscriptomics
