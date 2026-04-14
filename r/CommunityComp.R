# Community composition analysis
# All data loaded in LoadData.R

source("/Users/alicehong/kvich_meta/r/LoadData.R")

# Helper functions
fmt_p     <- function(p) ifelse(p < 0.001, sprintf("p = %.2e", p), sprintf("p = %.4f", p))
to_factor <- function(df) mutate(df, type = factor(type, levels = type_levels))

# Compute per-sample summary dfs (no_kneaddata)
total_df <- enframe(base::colSums(counts_scaled), name = "sample", value = "total_scaled") %>%
  left_join(meta %>% dplyr::select(sample, type, patient_ID), by = "sample")

shannon_df <- diversity(t(counts_unscaled), index = "shannon") %>%
  enframe(name = "sample", value = "shannon") %>%
  left_join(meta %>% dplyr::select(sample, type, patient_ID), by = "sample")

fuso_phylum_df <- base::colSums(counts_scaled[grepl("^Fusobacterium ", rownames(counts_scaled)), ]) %>%
  enframe(name = "sample", value = "fusobacteria") %>%
  left_join(meta %>% dplyr::select(sample, type, patient_ID), by = "sample")

species_df <- data.frame(
  sample    = colnames(counts_scaled),
  bfrag_rel = as.numeric(counts_scaled["Bacteroides fragilis", ]) / base::colSums(counts_scaled),
  fnuc_rel  = as.numeric(counts_scaled["Fusobacterium nucleatum", ]) / base::colSums(counts_scaled),
  check.names = FALSE
) %>% left_join(meta %>% dplyr::select(sample, type, patient_ID), by = "sample")

rel_long <- species_df %>%
  pivot_longer(c(bfrag_rel, fnuc_rel), names_to = "species", values_to = "rel") %>%
  mutate(
    group = paste0(ifelse(species == "bfrag_rel", "B. fragilis", "F. nucleatum"), "\n", type),
    group = factor(group, levels = c(
      "B. fragilis\nCRC", "B. fragilis\nPaired", "B. fragilis\nHealthy",
      "F. nucleatum\nCRC", "F. nucleatum\nPaired", "F. nucleatum\nHealthy"
    ))
  )

# KneadData versions
total_df_kd <- enframe(base::colSums(counts_scaled_kd), name = "sample", value = "total_scaled") %>%
  left_join(meta %>% dplyr::select(sample, type, patient_ID), by = "sample")

shannon_df_kd <- diversity(t(counts_unscaled_kd), index = "shannon") %>%
  enframe(name = "sample", value = "shannon") %>%
  left_join(meta %>% dplyr::select(sample, type, patient_ID), by = "sample")

fuso_phylum_df_kd <- base::colSums(counts_scaled_kd[grepl("^Fusobacterium ", rownames(counts_scaled_kd)), ]) %>%
  enframe(name = "sample", value = "fusobacteria") %>%
  left_join(meta %>% dplyr::select(sample, type, patient_ID), by = "sample")

species_df_kd <- data.frame(
  sample    = colnames(counts_scaled_kd),
  bfrag_rel = as.numeric(counts_scaled_kd["Bacteroides fragilis", ]) / base::colSums(counts_scaled_kd),
  fnuc_rel  = as.numeric(counts_scaled_kd["Fusobacterium nucleatum", ]) / base::colSums(counts_scaled_kd),
  check.names = FALSE
) %>% left_join(meta %>% dplyr::select(sample, type, patient_ID), by = "sample")

rel_long_kd <- species_df_kd %>%
  pivot_longer(c(bfrag_rel, fnuc_rel), names_to = "species", values_to = "rel") %>%
  mutate(
    group = paste0(ifelse(species == "bfrag_rel", "B. fragilis", "F. nucleatum"), "\n", type),
    group = factor(group, levels = c(
      "B. fragilis\nCRC", "B. fragilis\nPaired", "B. fragilis\nHealthy",
      "F. nucleatum\nCRC", "F. nucleatum\nPaired", "F. nucleatum\nHealthy"
    ))
  )

### Figure 2J recreation - Virulence gene expression vs F. nucleatum abundance

FOMA_ID <- "C7Y58_RS03545"
RADD_ID <- "C7Y58_RS00155"

# Extract raw featureCounts counts for fomA and radD, scale by sequencing depth
vir_counts <- counts_raw %>%
  filter(Geneid %in% c(FOMA_ID, RADD_ID)) %>%
  dplyr::select(Geneid, all_of(sample_cols)) %>%
  pivot_longer(-Geneid, names_to = "sample", values_to = "raw_count") %>%
  left_join(sf_df %>% dplyr::select(sample, scaling_factor), by = "sample") %>%
  mutate(scaled_count = raw_count * scaling_factor,
         gene = case_when(Geneid == FOMA_ID ~ "fomA", Geneid == RADD_ID ~ "radD"))

# Extract F. nucleatum abundance from the scaled Bracken table
fnuc_row   <- which(rownames(counts_scaled) == "Fusobacterium nucleatum")
fnuc_abund <- tibble(sample = colnames(counts_scaled),
                     fnuc_scaled = as.numeric(counts_scaled[fnuc_row, ]))

# Combine and label top 5% highest-F. nucleatum samples
plot_df <- vir_counts %>%
  left_join(fnuc_abund, by = "sample") %>%
  left_join(meta %>% dplyr::select(sample, type), by = "sample") %>%
  filter(!is.na(type), !is.na(fnuc_scaled)) %>%
  group_by(gene) %>%
  mutate(label = ifelse(fnuc_scaled > quantile(fnuc_scaled, 0.95), sample, NA_character_)) %>%
  ungroup()

p_fig2j <- ggplot(plot_df, aes(x = fnuc_scaled, y = scaled_count, colour = type)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_text_repel(aes(label = label), size = 2, max.overlaps = 20, show.legend = FALSE) +
  facet_wrap(~ gene, scales = "free_y") +
  scale_colour_manual(values = type_colors, name = "Sample type") +
  labs(x = "F. nucleatum abundance (Bracken scaled counts)",
       y = "Virulence gene expression\n(featureCounts scaled counts)") +
  theme_classic() + text_size

### Figure 4 recreation - Community composition by sample group

# Statistical test helpers
paired_w2 <- function(df, col) {
  df %>% filter(type %in% c("CRC", "Paired")) %>%
    dplyr::select(patient_ID, type, all_of(col)) %>%
    pivot_wider(names_from = type, values_from = all_of(col)) %>%
    drop_na() %>%
    summarise(p = wilcox.test(CRC, Paired, paired = TRUE)$p.value) %>%
    pull(p)
}

paired_t2 <- function(df, col) {
  df %>% filter(type %in% c("CRC", "Paired")) %>%
    dplyr::select(patient_ID, type, all_of(col)) %>%
    pivot_wider(names_from = type, values_from = all_of(col)) %>%
    drop_na() %>%
    summarise(p = t.test(CRC, Paired, paired = TRUE)$p.value) %>%
    pull(p)
}

ind_w2 <- function(df, col) wilcox.test(df[[col]][df$type == "CRC"],
                                         df[[col]][df$type == "Healthy"])$p.value
ind_t2 <- function(df, col) t.test(df[[col]][df$type == "CRC"],
                                    df[[col]][df$type == "Healthy"])$p.value

calc_fig4_stats <- function(total, shannon, fuso, species) {
  list(
    qA   = p.adjust(c(ind_w2(total,   "total_scaled"), paired_w2(total,   "total_scaled"))),
    qB   = p.adjust(c(ind_t2(shannon, "shannon"),      paired_t2(shannon, "shannon"))),
    qC   = p.adjust(c(ind_w2(fuso,    "fusobacteria"), paired_w2(fuso,    "fusobacteria"))),
    qD_b = p.adjust(c(ind_w2(species, "bfrag_rel"),    paired_w2(species, "bfrag_rel"))),
    qD_f = p.adjust(c(ind_w2(species, "fnuc_rel"),     paired_w2(species, "fnuc_rel")))
  )
}

# Apply factor ordering then compute stats
total_df       <- to_factor(total_df)
shannon_df     <- to_factor(shannon_df)
fuso_phylum_df <- to_factor(fuso_phylum_df)
species_df     <- to_factor(species_df)
rel_long       <- mutate(rel_long, type = factor(type, levels = type_levels))

total_df_kd       <- to_factor(total_df_kd)
shannon_df_kd     <- to_factor(shannon_df_kd)
fuso_phylum_df_kd <- to_factor(fuso_phylum_df_kd)
species_df_kd     <- to_factor(species_df_kd)
rel_long_kd       <- mutate(rel_long_kd, type = factor(type, levels = type_levels))

stats_raw <- calc_fig4_stats(total_df, shannon_df, fuso_phylum_df, species_df)
stats_kd  <- calc_fig4_stats(total_df_kd, shannon_df_kd, fuso_phylum_df_kd, species_df_kd)

make_fig4_panels <- function(total, shannon, fuso, rel, stats, main_title) {
  base_theme <- theme_classic(base_size = 7)

  pA <- ggplot(total, aes(x = type, y = total_scaled / 1e6)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 0.5) +
    geom_signif(comparisons = list(c("CRC", "Healthy"), c("CRC", "Paired")),
                annotations = fmt_p(stats$qA),
                step_increase = 0.12, tip_length = 0.02, textsize = 2) +
    ylim(0, max(total$total_scaled / 1e6) * 1.3) +
    labs(x = NULL, y = "Scaled Counts (x10^6)", title = "A  Bacterial scaled counts") + base_theme

  pB <- ggplot(shannon, aes(x = type, y = shannon)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 0.5) +
    geom_signif(comparisons = list(c("CRC", "Healthy"), c("CRC", "Paired")),
                annotations = fmt_p(stats$qB),
                step_increase = 0.1, tip_length = 0.02, textsize = 2) +
    ylim(0, max(shannon$shannon) * 1.3) +
    labs(x = NULL, y = "Shannon Index", title = "B  Alpha diversity") + base_theme

  pC <- ggplot(fuso, aes(x = type, y = fusobacteria)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 0.5) +
    geom_signif(comparisons = list(c("CRC", "Healthy"), c("CRC", "Paired")),
                annotations = fmt_p(stats$qC),
                step_increase = 0.12, tip_length = 0.02, textsize = 2) +
    ylim(0, max(fuso$fusobacteria) * 1.3) +
    labs(x = NULL, y = "Scaled Counts", title = "C  Fusobacterium counts") + base_theme

  pD <- ggplot(rel, aes(x = group, y = rel)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 0.5) +
    geom_signif(comparisons = list(
      c("B. fragilis\nCRC", "B. fragilis\nHealthy"),
      c("B. fragilis\nCRC", "B. fragilis\nPaired"),
      c("F. nucleatum\nCRC", "F. nucleatum\nHealthy"),
      c("F. nucleatum\nCRC", "F. nucleatum\nPaired")),
      annotations = fmt_p(c(stats$qD_b[1], stats$qD_b[2], stats$qD_f[1], stats$qD_f[2])),
      step_increase = 0.1, tip_length = 0.02, textsize = 2) +
    labs(x = NULL, y = "Relative abundance",
         title = "D  B. fragilis / F. nucleatum relative abundance") +
    base_theme +
    theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 6))

  (pA | pB) / (pC | pD) + plot_annotation(title = main_title)
}

p_fig4_raw <- make_fig4_panels(total_df, shannon_df, fuso_phylum_df, rel_long,
                                stats_raw, "Raw unmapped reads (n = 119)")
p_fig4_kd  <- make_fig4_panels(total_df_kd, shannon_df_kd, fuso_phylum_df_kd, rel_long_kd,
                                stats_kd, "KneadData-cleaned reads (n = 86)")
