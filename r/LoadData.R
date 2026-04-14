source("/Users/alicehong/kvich_meta/r/LoadLibs.R")
source("/Users/alicehong/kvich_meta/r/HelperFunctions.R")

# Constants ──────────────────────────────────────────────────────────────────
BFRAG_CHR   <- "NZ_CP069563.1"
FNUC_CHR    <- "NZ_CP028101.1"
type_levels <- c("CRC", "Healthy", "Paired")

# Metadata ───────────────────────────────────────────────────────────────────
meta <- read.delim("/Users/alicehong/kvich_meta/meta.tsv") %>%
  mutate(type = factor(type, levels = c("CRC", "Paired", "Healthy")))

# FastQC general stats (pre-KneadData): total sequences, GC, length, duplicates
multiqc <- read.delim("/Users/alicehong/kvich_meta/kneaddata/qc/multiqc_pre/multiqc_report_data/multiqc_general_stats.txt") %>%
  mutate(Sample = clean_sample_ids(Sample)) %>%
  extract(Sample, into = c("sample", "read"), regex = "^(La\\d+)_(.*)$") %>%
  pivot_wider(id_cols = sample, names_from = read, values_from = -c(sample, read))

# Per-base sequence quality (pre-KneadData): mean Phred per sample
qual <- read.delim("/Users/alicehong/kvich_meta/kneaddata/qc/multiqc_pre/multiqc_report_data/fastqc_per_base_sequence_quality_plot.txt",
                   check.names = FALSE) %>%
  pivot_longer(-Sample, names_to = "pos_idx", values_to = "tuple") %>%
  filter(!is.na(tuple) & tuple != "") %>%
  mutate(Sample  = clean_sample_ids(Sample),
         quality = as.numeric(str_extract(tuple, "(?<=,\\s)[0-9.]+"))) %>%
  extract(Sample, into = c("sample", "read"), regex = "^(La\\d+)_(.*)$") %>%
  group_by(sample) %>%
  summarise(mean_quality = mean(quality, na.rm = TRUE))

# Enrich meta with FastQC metrics
meta <- meta %>%
  left_join(multiqc, by = "sample") %>%
  mutate(
    total_sequences = rowMeans(dplyr::select(., contains("fastqc.total_sequences")), na.rm = TRUE),
    percent_gc      = rowMeans(dplyr::select(., contains("fastqc.percent_gc")),      na.rm = TRUE),
    avg_seq_length  = rowMeans(dplyr::select(., contains("fastqc.avg_sequence_length")), na.rm = TRUE),
    percent_dups    = rowMeans(dplyr::select(., contains("fastqc.percent_duplicates")), na.rm = TRUE)
  ) %>%
  left_join(qual, by = "sample")

# featureCounts ──────────────────────────────────────────────────────────────
counts_raw    <- load_featurecounts("/Users/alicehong/kvich_meta/no_kneaddata/alignment/counts/counts.txt")
kd_counts_raw <- load_featurecounts("/Users/alicehong/kvich_meta/kneaddata/alignment/counts/counts.txt")
sample_cols   <- setdiff(colnames(counts_raw), c("Geneid", "Chr"))
gene_annot    <- read.delim("/Users/alicehong/kvich_meta/no_kneaddata/alignment/gene_annotations.tsv") %>%
  filter(product != "")

# Flagstat ───────────────────────────────────────────────────────────────────
# Per-sample .flagstat files → total read pairs per sample (used in bfrag_fnuc)
flagstat_df <- list.files("/Users/alicehong/kvich_meta/no_kneaddata/qc/flagstat",
                          pattern = "\\.flagstat$", full.names = TRUE) %>%
  map_dfr(parse_flagstat)

kd_flagstat_df <- read.delim("/Users/alicehong/kvich_meta/kneaddata/qc/multiqc_bwa/bwa_alignment_qc_data/multiqc_samtools_flagstat.txt") %>%
  transmute(sample = clean_sample_ids(Sample), total_fragments = read1_passed)

# Full multiqc flagstat tables (used in Flagstat.R)
flagstat_kd   <- read.delim("/Users/alicehong/kvich_meta/kneaddata/qc/multiqc_bwa/bwa_alignment_qc_data/multiqc_samtools_flagstat.txt") %>%
  mutate(pipeline = "KneadData")
flagstat_nokd <- read.delim("/Users/alicehong/kvich_meta/no_kneaddata/qc/multiqc/bwa_alignment_qc_data/multiqc_samtools_flagstat.txt") %>%
  mutate(pipeline = "No KneadData")

flagstat <- bind_rows(flagstat_kd, flagstat_nokd) %>%
  mutate(sample   = clean_sample_ids(Sample),
         pipeline = factor(pipeline, levels = c("No KneadData", "KneadData"))) %>%
  left_join(meta, by = "sample") %>%
  filter(!is.na(type)) %>%
  dplyr::select(sample, type, pipeline,
                total_reads         = total_passed,
                mapped              = mapped_passed,
                pct_mapped          = mapped_passed_pct,
                properly_paired     = properly.paired_passed,
                pct_properly_paired = properly.paired_passed_pct,
                singletons          = singletons_passed,
                pct_singletons      = singletons_passed_pct)

# KneadData read counts ──────────────────────────────────────────────────────
kd <- read.delim("/Users/alicehong/kvich_meta/kneaddata/kneaddata_read_counts.tsv",
                 check.names = FALSE) %>%
  mutate(reads_lost_trim  = `raw pair1` - `trimmed pair1`,
         reads_lost_decon = `trimmed pair1` - `decontaminated hg_39 pair1`,
         final_M  = `final pair1` / 1e6,
         flagged  = `final pair1` < 1e5) %>%
  arrange(`final pair1`) %>%
  mutate(Sample = factor(Sample, levels = Sample))

# Bracken normalized counts ──────────────────────────────────────────────────
rename_lasse <- function(df) { colnames(df) <- clean_sample_ids(colnames(df)); df }

counts_scaled   <- read.delim("no_kneaddata/community_composition/normalized/counts_scaled.tsv",
                               row.names = 1, check.names = FALSE) %>% rename_lasse()
counts_unscaled <- read.delim("no_kneaddata/community_composition/normalized/counts_unscaled.tsv",
                               row.names = 1, check.names = FALSE) %>% rename_lasse()
keep            <- intersect(colnames(counts_scaled), meta$sample)
counts_scaled   <- counts_scaled[, keep]
counts_unscaled <- counts_unscaled[, keep]

sf_df <- read.delim("no_kneaddata/community_composition/normalized/scaling_factors.tsv") %>%
  mutate(sample = clean_sample_ids(sample))

counts_scaled_kd   <- read.delim("kneaddata/community_composition/normalized/counts_scaled.tsv",
                                  row.names = 1, check.names = FALSE)
counts_unscaled_kd <- read.delim("kneaddata/community_composition/normalized/counts_unscaled.tsv",
                                  row.names = 1, check.names = FALSE)
keep_kd            <- intersect(colnames(counts_scaled_kd), meta$sample)
counts_scaled_kd   <- counts_scaled_kd[, keep_kd]
counts_unscaled_kd <- counts_unscaled_kd[, keep_kd]

sf_df_kd <- read.delim("kneaddata/community_composition/normalized/scaling_factors.tsv") %>%
  mutate(sample = clean_sample_ids(sample))

# Kraken2 / Bracken classification output ────────────────────────────────────
kraken_files <- list.files("/Users/alicehong/kvich_meta/kneaddata/community_composition/kraken2",
                           pattern = "\\.report$", full.names = TRUE)
kraken <- setNames(kraken_files,
                   str_extract(clean_sample_ids(basename(kraken_files)), "La\\d+")) %>%
  map_dfr(~ read.delim(.x, header = FALSE,
                       col.names = c("pct", "clade_reads", "direct_reads", "rank", "taxid", "name")),
          .id = "sample") %>%
  left_join(meta, by = "sample")

bracken_files <- list.files("/Users/alicehong/kvich_meta/kneaddata/community_composition/bracken",
                            pattern = "\\.bracken$", full.names = TRUE)
bracken <- setNames(bracken_files,
                    str_extract(clean_sample_ids(basename(bracken_files)), "La\\d+")) %>%
  map_dfr(~ read.delim(.x), .id = "sample") %>%
  left_join(meta, by = "sample")

# ── Pre-computed DE results (from bfrag_fnuc_gene_expression.R) ───────────────

OUT <- "/Users/alicehong/kvich_meta/r/output"
if (file.exists(file.path(OUT, "deseq2_crc_vs_healthy.tsv"))) {
  res_crc_healthy    <- read.delim(file.path(OUT, "deseq2_crc_vs_healthy.tsv"))
  res_crc_paired     <- read.delim(file.path(OUT, "deseq2_crc_vs_paired.tsv"))
  res_1M_crc_healthy <- read.delim(file.path(OUT, "deseq2_1M_crc_vs_healthy.tsv"))
  res_1M_crc_paired  <- read.delim(file.path(OUT, "deseq2_1M_crc_vs_paired.tsv"))
  aldex_crc_healthy  <- read.delim(file.path(OUT, "aldex2_crc_vs_healthy.tsv"))
  aldex_crc_paired   <- read.delim(file.path(OUT, "aldex2_crc_vs_paired.tsv"))}
