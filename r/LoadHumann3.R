
# meta and type_levels loaded in LoadData.R
source("/Users/alicehong/kvich_meta/r/LoadData.R")

# Script to analyze humann3 results

# Load .tsvs from 3.3!

KD_MERGED <- "kneaddata/functional_profiling/humann3_merged"

# Unstratified pathway abundance (CPM) - 390 pathways x 84 samples
path_abund <- read.delim(file.path(KD_MERGED, "all_pathabundance_unstratified.tsv"),check.names = FALSE) %>%
  rename_with(~"pathway", 1) %>%
  filter(!pathway %in% c("UNMAPPED", "UNINTEGRATED"))

# Pathway coverage (0-1) - same dimensions, used to filter low-confidence pathways
path_cov <- read.delim(file.path(KD_MERGED, "all_pathcoverage.tsv"),check.names = FALSE) %>%
  rename_with(~"pathway", 1) %>%
  filter(!pathway %in% c("UNMAPPED", "UNINTEGRATED"))

# Stratified pathway abundance - pathway|species rows
path_strat <- read.delim(file.path(KD_MERGED, "all_pathabundance_stratified.tsv"),check.names = FALSE) %>%
  rename_with(~"pathway_species", 1) %>%
  filter(!grepl("^UNMAPPED|^UNINTEGRATED", pathway_species))

# Clean up sample column names: strip "_cat_Abundance" suffix, remap Lasse -> La
clean_humann_cols <- function(df) {
  colnames(df) <- colnames(df) %>%
    str_remove("_cat_Abundance$") %>%
    str_remove("_cat_Coverage$") %>%
    str_replace("^Lasse(\\d+)$", "La\\1")
  df}
path_abund <- clean_humann_cols(path_abund)
path_cov   <- clean_humann_cols(path_cov)
path_strat <- clean_humann_cols(path_strat)

# Sample columns present in both HUMANn3 output and metadata
h3_samples <- intersect(colnames(path_abund)[-1], meta$sample)
meta_h3    <- meta %>% filter(sample %in% h3_samples) %>% mutate(type = factor(type, levels = type_levels))

# Coverage-filtered pathway matrix: keep pathways with mean coverage >= 0.1
cov_mat  <- path_cov  %>% column_to_rownames("pathway") %>% dplyr::select(all_of(h3_samples))
path_mat <- path_abund %>% column_to_rownames("pathway") %>% dplyr::select(all_of(h3_samples))
keep_paths <- rowMeans(cov_mat >= 0.1, na.rm = TRUE) >= 0.1
path_mat_filt <- path_mat[keep_paths, ]
path_mat_filt[is.na(path_mat_filt)] <- 0

# Long format for stats/plots
path_long <- path_mat_filt %>%
  rownames_to_column("pathway") %>%
  pivot_longer(-pathway, names_to = "sample", values_to = "cpm") %>%
  left_join(meta_h3, by = "sample")

############
## PCoA
path_mat_num <- path_mat_filt %>% as.matrix()

mode(path_mat_num) <- "numeric"
path_mat_num[is.na(path_mat_num)] <- 0

bc_dist <- vegdist(t(path_mat_num), method = "bray")

pcoa <- cmdscale(bc_dist, k = 2, eig = TRUE)
pct_var <- round(pcoa$eig[1:2] / sum(pcoa$eig[pcoa$eig > 0]) * 100, 1)

# plotting pcoa2
data.frame(sample = rownames(pcoa$points),PC1 = pcoa$points[,1],PC2 = pcoa$points[,2]) %>%
  left_join(meta_h3, by = "sample") %>%
  ggplot(aes(PC1, PC2, colour = type)) +
  geom_point(size = 2.5, alpha = 0.8) +
  stat_ellipse(level = 0.75, linetype = "dashed") +
  theme_classic()
