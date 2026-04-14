source("/Users/alicehong/kvich_meta/r/LoadData.R")

# res_crc_healthy / res_crc_paired loaded from r/output via LoadData.R

SIG <- 0.05
LFC <- 1  # abs(log2FoldChange) cutoff

org_colors <- c("B. fragilis" = "#A0522D", "F. nucleatum" = "#1B6CA8")

# ── Filter significant genes ──────────────────────────────────────────────────
sig_healthy <- res_crc_healthy %>%
  filter(!is.na(padj), padj < SIG, abs(log2FoldChange) > LFC)

sig_paired <- res_crc_paired %>%
  filter(!is.na(padj), padj < SIG, abs(log2FoldChange) > LFC)

cat(sprintf("CRC vs Healthy: %d significant genes\n", nrow(sig_healthy)))
cat(sprintf("CRC vs Paired:  %d significant genes\n", nrow(sig_paired)))
cat(sprintf("Shared across both contrasts: %d genes\n",
            length(intersect(sig_healthy$Geneid, sig_paired$Geneid))))

# ── 1. Volcano plots ──────────────────────────────────────────────────────────
make_volcano <- function(res, title, n_label = 15) {
  df <- res %>%
    filter(!is.na(padj)) %>%
    mutate(
      sig = padj < SIG & abs(log2FoldChange) > LFC,
      direction = case_when(
        sig & log2FoldChange > 0 ~ "Up in CRC",
        sig & log2FoldChange < 0 ~ "Down in CRC",
        TRUE ~ "NS"
      ),
      label = if_else(sig, coalesce(gene, product, Geneid), NA_character_)
    )

  top_genes <- df %>% filter(sig) %>% slice_min(padj, n = n_label, with_ties = FALSE)

  ggplot(df, aes(x = log2FoldChange, y = -log10(padj), color = organism)) +
    geom_point(data = filter(df, direction == "NS"),
               color = "grey80", size = 0.7, alpha = 0.5) +
    geom_point(data = filter(df, direction != "NS"),
               size = 1.2, alpha = 0.8) +
    geom_vline(xintercept = c(-LFC, LFC), linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = -log10(SIG), linetype = "dashed", color = "grey50") +
    ggrepel::geom_text_repel(
      data = top_genes,
      aes(label = label),
      size = 2.8, max.overlaps = 20, show.legend = FALSE
    ) +
    scale_color_manual(values = org_colors, na.value = "grey60", na.translate = FALSE, name = NULL) +
    labs(title = title, x = "Log2 Fold Change", y = "-log10(adj. p-value)") +
    theme_classic() +
    theme(legend.position = "bottom")
}

p_volcano_healthy <- make_volcano(res_crc_healthy, "CRC vs Healthy")
p_volcano_paired  <- make_volcano(res_crc_paired,  "CRC vs Paired")

p_volcano_healthy + p_volcano_paired + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# ── 2. Venn diagram: overlap between the two contrasts ────────────────────────
venn_list <- list(
  "CRC vs Healthy" = sig_healthy$Geneid,
  "CRC vs Paired"  = sig_paired$Geneid
)

grid.newpage()
grid.draw(venn.diagram(
  venn_list,
  filename  = NULL,
  fill      = c("#e6194b", "#4363d8"),
  alpha     = 0.45,
  cex       = 1.1,
  cat.cex   = 0.9,
  cat.dist  = 0.05,
  main      = "Significant DE genes (DESeq2)",
  main.cex  = 1.1
))

# ── 3. Heatmap of top significant genes ───────────────────────────────────────
# Pick top N genes by padj from either contrast, then VST-normalise counts
TOP_N <- 25

top_genes_heatmap <- bind_rows(sig_healthy, sig_paired) %>%
  arrange(padj) %>%
  distinct(Geneid, .keep_all = TRUE) %>%
  slice_head(n = TOP_N)

# Build VST count matrix for these genes across all samples with a type
hm_meta <- meta %>% filter(sample %in% colnames(counts_raw), !is.na(type)) %>%
  arrange(type)

hm_mat <- counts_raw %>%
  filter(Geneid %in% top_genes_heatmap$Geneid) %>%
  column_to_rownames("Geneid") %>%
  dplyr::select(-Chr) %>%
  dplyr::select(all_of(hm_meta$sample)) %>%
  mutate(across(everything(), as.integer))

# log2 CPM — use full library sizes so sparse heatmap genes don't cause /0
full_lib <- counts_raw %>%
  column_to_rownames("Geneid") %>%
  dplyr::select(-Chr) %>%
  dplyr::select(all_of(hm_meta$sample)) %>%
  mutate(across(everything(), as.integer)) %>%
  base::colSums()

# drop samples with zero library size
keep_samp <- names(full_lib)[full_lib > 0]
hm_mat    <- hm_mat[, keep_samp, drop = FALSE]
full_lib  <- full_lib[keep_samp]

hm_mat  <- as.matrix(hm_mat)
cpm_mat <- sweep(hm_mat, 2, full_lib / 1e6, "/")
vst_mat <- log2(cpm_mat + 1)
vst_mat[!is.finite(vst_mat)] <- 0   # guard against any residual Inf/NaN

# Row labels: gene name or product (truncated)
row_labels <- top_genes_heatmap %>%
  dplyr::select(Geneid, gene, product, organism) %>%
  mutate(label = coalesce(
    na_if(gene, ""),
    str_trunc(product, 50),
    Geneid
  )) %>%
  dplyr::select(Geneid, label, organism)

rownames(vst_mat) <- row_labels$label[match(rownames(vst_mat), row_labels$Geneid)]

# Column annotation: sample type (aligned to filtered sample set)
col_ann <- data.frame(
  Type = hm_meta$type[hm_meta$sample %in% keep_samp],
  row.names = hm_meta$sample[hm_meta$sample %in% keep_samp]
)
ann_colors <- list(
  Type     = c(CRC = "#d73027", Healthy = "#4575b4", Paired = "#fdae61"),
  Organism = org_colors
)

# Row annotation: organism
row_ann <- data.frame(
  Organism = row_labels$organism[match(
    rownames(vst_mat),
    row_labels$label
  )],
  row.names = rownames(vst_mat)
)

ph_gt <- pheatmap::pheatmap(
  vst_mat,
  annotation_col    = col_ann,
  annotation_row    = row_ann,
  annotation_colors = ann_colors,
  scale             = "row",
  show_colnames     = FALSE,
  fontsize_row      = 7,
  cluster_rows      = FALSE,
  clustering_distance_cols = "euclidean",
  angle_row         = 0,
  main              = paste0("Top ", TOP_N, " DE genes (log2 CPM, row-scaled)"),
  silent            = TRUE
)$gtable

# Make row annotation label ("Organism") horizontal
for (i in seq_along(ph_gt$grobs)) {
  g <- ph_gt$grobs[[i]]
  if (inherits(g, "text") && length(g$label) == 1 && g$label == "Organism") {
    ph_gt$grobs[[i]]$rot   <- 0
    ph_gt$grobs[[i]]$hjust <- 0
    ph_gt$grobs[[i]]$x     <- unit(0, "npc")
  }
}

p_heatmap <- wrap_elements(ph_gt)

p_heatmap

# ── 4. Top genes table ────────────────────────────────────────────────────────
top_table <- function(df, title) {
  df %>%
    dplyr::select(gene, product, organism, log2FoldChange, padj,
                  n_samples_with_counts, pct_total_reads) %>%
    mutate(across(c(log2FoldChange, pct_total_reads), ~ round(.x, 2)),
           padj = signif(padj, 3)) %>%
    slice_min(padj, n = 20, with_ties = FALSE) %>%
    knitr::kable(caption = title, digits = 3)
}

top_table(sig_healthy, "Top DE genes: CRC vs Healthy")
top_table(sig_paired,  "Top DE genes: CRC vs Paired")
