# Script to explore results of bwa-mem to b. fragilis and f. nucleatum reference genomes
# All data (meta, counts_raw, kd_counts_raw, gene_annot, flagstat_df, kd_flagstat_df) loaded in LoadData.R

source("/Users/alicehong/bmeg524/r/LoadData.R")

# Similar colors to as in supplemental of paper!
fill_scale <- scale_fill_manual(name = "Mapping",values = c("B. fragilis" = "#A0522D", "F. nucleatum" = "#1B6CA8", "Unmapped" = "#BEBEBE"))

# Local helpers

build_reads_df <- function(flagstat_df, counts_raw) {
  sc <- setdiff(colnames(counts_raw), c("Geneid", "Chr"))
  bfrag <- base::colSums(counts_raw %>% filter(Chr == BFRAG_CHR) %>% dplyr::select(all_of(sc)))
  fnuc  <- base::colSums(counts_raw %>% filter(Chr == FNUC_CHR)  %>% dplyr::select(all_of(sc)))
  flagstat_df %>%
    left_join(tibble(sample = names(bfrag), bfrag = bfrag), by = "sample") %>%
    left_join(tibble(sample = names(fnuc),  fnuc  = fnuc),  by = "sample") %>%
    mutate(across(c(bfrag, fnuc), ~ replace_na(.x, 0)),
           unmapped = pmax(total_fragments - bfrag - fnuc, 0)) %>%
    left_join(meta %>% dplyr::select(sample, type), by = "sample") %>%
    filter(!is.na(type))
}

############################################################
# Plot to recreate form paper
plot_s6 <- function(reads_df, title_suffix = "") {
  rl <- reads_df %>%
    pivot_longer(c(bfrag, fnuc, unmapped), names_to = "mapping", values_to = "reads") %>%
    mutate(reads_M = reads / 1e6,
           mapping = factor(mapping, levels = c("bfrag", "fnuc", "unmapped"),
                            labels = c("B. fragilis", "F. nucleatum", "Unmapped")))
  base <- list(fill_scale, scale_y_continuous(expand = expansion(mult = c(0, 0.05))),
               theme_classic())
  pa <- rl %>% sort_by_total() %>%
    ggplot(aes(x = sample, y = reads_M, fill = mapping)) +
    geom_col(width = 1) + geom_hline(yintercept = 1, linetype = "dashed") + base +
    labs(x = "Sample", y = "Non-Human Reads (M)", title = paste0("A  All samples", title_suffix)) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + text_size
  pb <- rl %>% group_by(sample) %>% filter(sum(reads_M) > 1) %>% ungroup() %>%
    sort_by_total() %>%
    ggplot(aes(x = sample, y = reads_M, fill = mapping)) +
    geom_col() + base +
    labs(x = "Sample", y = "Non-Human Reads (M)", title = paste0("B  Samples with >1M reads", title_suffix)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + text_size
  pa + pb + plot_layout(guides = "collect")
}

make_table_s2 <- function(counts_raw, reads_df) {
  crc_1M <- reads_df %>% filter(type == "CRC", total_fragments >= 1e6) %>% pull(sample)
  counts_mat <- counts_raw %>% column_to_rownames("Geneid") %>% dplyr::select(-Chr) %>%
    dplyr::select(any_of(crc_1M)) %>% mutate(across(everything(), as.integer))
  run_vst <- function(gene_ids) {
    mat <- counts_mat[rownames(counts_mat) %in% gene_ids, ]
    mat <- mat[rowSums(mat) > 0, ]
    coldata <- data.frame(row.names = colnames(mat), sample = colnames(mat))
    dds <- DESeqDataSetFromMatrix(mat, colData = coldata, design = ~ 1)
    rowMeans(assay(varianceStabilizingTransformation(dds, blind = TRUE)))
  }
  make_top20 <- function(avg_expr, species_label) {
    data.frame(gene_id = names(avg_expr), avg_vst = avg_expr) %>%
      inner_join(gene_annot, by = "gene_id") %>%
      arrange(desc(avg_vst)) %>% slice_head(n = 5) %>%
      transmute(Species = species_label, Gene = gene, Product = product, Avg_VST = round(avg_vst, 2))
  }
  bind_rows(
    make_top20(run_vst(counts_raw %>% filter(Chr == BFRAG_CHR) %>% pull(Geneid)), "B. fragilis"),
    make_top20(run_vst(counts_raw %>% filter(Chr == FNUC_CHR)  %>% pull(Geneid)), "F. nucleatum")
  )
}

# No KneadData plot/table
reads_df <- build_reads_df(flagstat_df, counts_raw)

plot_s6(reads_df)
table_s2 <- make_table_s2(counts_raw, reads_df)
table_s2

# KneadData plot/table
kd_reads_df    <- build_reads_df(kd_flagstat_df, kd_counts_raw)

plot_s6(kd_reads_df, " (KneadData)")
kd_table_s2    <- make_table_s2(kd_counts_raw, kd_reads_df)
kd_table_s2

# ############################################################
# DESeq2: CRC vs Healthy and CRC vs paired
# Alread ran and written out

# All three conditions — fit one model, extract both contrasts
# dds_meta <- meta %>%
#   filter(sample %in% colnames(counts_raw), !is.na(type)) %>%
#   mutate(type = droplevels(type))
# 
# dds_mat <- counts_raw %>%
#   column_to_rownames("Geneid") %>%
#   dplyr::select(-Chr) %>%
#   dplyr::select(all_of(dds_meta$sample)) %>%
#   mutate(across(everything(), as.integer))
# 
# dds <- DESeqDataSetFromMatrix(dds_mat,
#                               colData = column_to_rownames(dds_meta, "sample"),
#                               design  = ~ type) %>% DESeq()

## Helper to extract results and annotate with organism
# get_res <- function(dds, contrast_label, contrast) {
#   cnt <- counts(dds)
#   gene_stats <- tibble(Geneid = rownames(cnt),
#                        n_samples_with_counts = rowSums(cnt > 0),
#                        pct_total_reads       = rowSums(cnt) / sum(cnt) * 100)
#   results(dds, contrast = contrast, alpha = 0.05) %>%
#     as.data.frame() %>%
#     rownames_to_column("Geneid") %>%
#     left_join(gene_stats, by = "Geneid") %>%
#     left_join(counts_raw %>% dplyr::select(Geneid, Chr), by = "Geneid") %>%
#     left_join(gene_annot, by = c("Geneid" = "gene_id")) %>%
#     mutate(organism = case_when(Chr == BFRAG_CHR ~ "B. fragilis",
#                                 Chr == FNUC_CHR  ~ "F. nucleatum"),
#            contrast = contrast_label) %>%
#     arrange(padj)
# }
# 
# res_crc_healthy <- get_res(dds, "CRC vs Healthy", c("type", "CRC", "Healthy"))
# res_crc_paired  <- get_res(dds, "CRC vs Paired",  c("type", "CRC", "Paired"))

# Volcano plots side by side
bind_rows(res_crc_healthy, res_crc_paired) %>%
  filter(!is.na(padj)) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj), color = organism)) +
  geom_point(size = 0.8, alpha = 0.6) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  scale_color_manual(values = c("B. fragilis" = "#A0522D", "F. nucleatum" = "#1B6CA8"),
                     na.value = "grey70", name = NULL) +
  facet_wrap(~ contrast) +
  labs(title = "DESeq2 CRC comparisons", x = "Log2 Fold Change", y = "-log10(adj. p-value)") +
  theme_classic() + text_size

############################################################
# DESeq2 but motr htan 1M reads only
## Already ran and written out

# dds_meta_1M <- meta %>%
#   filter(sample %in% colnames(counts_raw), !is.na(type),
#          sample %in% (reads_df %>% filter(total_fragments >= 1e6) %>% pull(sample))) %>%
#   mutate(type = droplevels(type))
# 
# dds_mat_1M <- counts_raw %>%
#   column_to_rownames("Geneid") %>%
#   dplyr::select(-Chr) %>%
#   dplyr::select(all_of(dds_meta_1M$sample)) %>%
#   mutate(across(everything(), as.integer))
# 
# dds_1M <- DESeqDataSetFromMatrix(dds_mat_1M,
#                                  colData = column_to_rownames(dds_meta_1M, "sample"),
#                                  design  = ~ type) %>%
#   DESeq()

# res_1M_crc_healthy <- get_res(dds_1M, "CRC vs Healthy", c("type", "CRC", "Healthy"))
# res_1M_crc_paired  <- get_res(dds_1M, "CRC vs Paired",  c("type", "CRC", "Paired"))

bind_rows(res_1M_crc_healthy, res_1M_crc_paired) %>%
  filter(!is.na(padj)) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj), color = organism)) +
  geom_point(size = 0.8, alpha = 0.6) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  scale_color_manual(values = c("B. fragilis" = "#A0522D", "F. nucleatum" = "#1B6CA8"),
                     na.value = "grey70", name = NULL) +
  facet_wrap(~ contrast) +
  labs(title = "DESeq2 CRC comparisons (>=1M reads)", x = "Log2 Fold Change", y = "-log10(adj. p-value)") +
  theme_classic() + text_size

############################################################
# ALDEx2: CRC vs Healthy and CRC vs Paired
# No depth filter — ALDEx2 handles variable depth via CLR; keep gene filter only
aldex_mat <- counts_raw %>%
  column_to_rownames("Geneid") %>%
  dplyr::select(-Chr) %>%
  dplyr::select(any_of(reads_df$sample)) %>%
  mutate(across(everything(), as.integer)) %>%
  .[rowSums(. >= 10) >= 2, ]


# Already ran and written out!
# aldex_crc_healthy <- run_aldex2("CRC", "Healthy", "CRC vs Healthy")
# aldex_crc_paired  <- run_aldex2("CRC", "Paired",  "CRC vs Paired")

# Effect plot: effect size vs median CLR abundance (ALDEx2 standard visualisation)
bind_rows(aldex_crc_healthy, aldex_crc_paired) %>%
  ggplot(aes(x = rab.all, y = effect, color = organism)) +
  geom_point(size = 0.8, alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("B. fragilis" = "#A0522D", "F. nucleatum" = "#1B6CA8"),
                     na.value = "grey70", name = NULL) +
  facet_wrap(~ contrast) +
  labs(title = "ALDEx2 CRC comparisons",
       x = "Median CLR Abundance (all samples)", y = "Effect Size") +
  theme_classic() + text_size

# ── Write results to disk ─────────────────────────────────────────────────────
# OUT <- "/Users/alicehong/bmeg524/r/output"
# write.table(res_crc_healthy,    file.path(OUT, "deseq2_crc_vs_healthy.tsv"),    sep = "\t", row.names = FALSE, quote = FALSE)
# write.table(res_crc_paired,     file.path(OUT, "deseq2_crc_vs_paired.tsv"),     sep = "\t", row.names = FALSE, quote = FALSE)
# write.table(res_1M_crc_healthy, file.path(OUT, "deseq2_1M_crc_vs_healthy.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
# write.table(res_1M_crc_paired,  file.path(OUT, "deseq2_1M_crc_vs_paired.tsv"),  sep = "\t", row.names = FALSE, quote = FALSE)
# write.table(aldex_crc_healthy,  file.path(OUT, "aldex2_crc_vs_healthy.tsv"),    sep = "\t", row.names = FALSE, quote = FALSE)
# write.table(aldex_crc_paired,   file.path(OUT, "aldex2_crc_vs_paired.tsv"),     sep = "\t", row.names = FALSE, quote = FALSE)
