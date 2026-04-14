# Checking bracken output
# kraken, bracken, meta (with mean_quality) loaded in LoadData.R

source("/Users/alicehong/kvich_meta/r/LoadData.R")

# Reads assigned at each major rank (unclassified + root-level domain breakdown)
kraken_levels <- kraken %>%
  filter(rank %in% c("U", "D", "P", "C", "O", "F", "G", "S")) %>%
  mutate(rank = factor(rank, levels = c("U", "D", "P", "C", "O", "F", "G", "S"),
                       labels = c("Unclassified", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")))

# Total classified vs unclassified reads per sample
kraken_summary <- kraken %>%
  filter(rank %in% c("U", "R")) %>%
  mutate(label = ifelse(rank == "U", "Unclassified", "Classified")) %>%
  dplyr::select(sample, type, label, clade_reads)

p_kraken_classified <- ggplot(kraken_summary, aes(x = sample, y = clade_reads / 1e6, fill = label)) +
  geom_col() +
  facet_grid(~ type, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = c("Classified" = "#4C9BE8", "Unclassified" = "#D3D3D3")) +
  labs(title = "Kraken2: Classified vs Unclassified Reads per Sample",
       x = NULL, y = "Reads (millions)", fill = NULL) +
  theme_classic() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

# % classified reads per sample
kraken_pct <- kraken_summary %>%
  pivot_wider(names_from = label, values_from = clade_reads) %>%
  mutate(
    total_reads    = (Classified + Unclassified) / 1e6,
    pct_classified = Classified / (Classified + Unclassified) * 100
  ) %>%
  dplyr::select(sample, type, total_reads, pct_classified)

# Compare % classified to mean sequence quality per sample
kraken_pct_qual <- kraken_pct %>%
  left_join(meta %>% dplyr::select(sample, mean_quality), by = "sample")

#####################
# Plots
p_pct_vs_depth <- ggplot(kraken_pct_qual, aes(x = total_reads, y = pct_classified, color = type)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.7) +
  scale_color_manual(values = type_colors) +
  labs(title = "% Classified vs Sequencing Depth",
       x = "Total Reads (millions)", y = "% Classified", color = NULL) +
  theme_classic()

p_pct_vs_qual <- ggplot(kraken_pct_qual, aes(x = mean_quality, y = pct_classified, color = type)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.7) +
  scale_color_manual(values = type_colors) +
  labs(title = "% Classified vs Mean Sequence Quality",
       x = "Mean Phred Quality Score", y = "% Classified", color = NULL) +
  theme_classic()

p_pct_vs_qual

p_kraken_pct <- ggplot(kraken_pct, aes(x = type, y = pct_classified, fill = type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(aes(color = type), width = 0.15, size = 1.5, alpha = 0.8) +
  scale_fill_manual(values = type_colors) +
  scale_color_manual(values = type_colors) +
  labs(title = "Kraken2: % Classified Reads by Condition",
       x = NULL, y = "% Classified") +
  theme_classic() +
  theme(legend.position = "none")

(p_kraken_classified / p_kraken_pct) +
  plot_annotation(title = "Kraken2 Classification Summary")

# --- Bracken output: species-level re-estimated reads per sample ---
bracken_dir <- "/Users/alicehong/kvich_meta/kneaddata/community_composition/bracken"

bracken_files <- list.files(bracken_dir, pattern = "\\.bracken$", full.names = TRUE)
bracken <- setNames(bracken_files, str_extract(clean_sample_ids(basename(bracken_files)), "La\\d+")) %>%
  map_dfr(~ read.delim(.x), .id = "sample") %>%
  left_join(meta, by = "sample")

# Total bracken-assigned reads per sample (sum of new_est_reads)
bracken_totals <- bracken %>%
  group_by(sample, type) %>%
  summarise(total_reads = sum(new_est_reads), .groups = "drop")

bracken_pct <- bracken_totals %>%
  left_join(kraken_pct %>% dplyr::select(sample, input_reads = total_reads), by = "sample") %>%
  mutate(pct_assigned = (total_reads / 1e6) / input_reads * 100) %>%
  dplyr::select(sample, type, input_reads, total_reads, pct_assigned)

bracken_totals_qual <- bracken_totals %>%
  left_join(meta %>% dplyr::select(sample, mean_quality), by = "sample")

p_bracken_vs_qual <- ggplot(bracken_totals_qual, aes(x = mean_quality, y = total_reads / 1e6, color = type)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.7) +
  scale_color_manual(values = type_colors) +
  labs(title = "Bracken Assigned Reads vs Mean Sequence Quality",
       x = "Mean Phred Quality Score", y = "Bracken Assigned Reads (millions)", color = NULL) +
  theme_classic()

p_bracken_vs_depth <- bracken_totals %>%
  left_join(kraken_pct %>% dplyr::select(sample, input_reads = total_reads), by = "sample") %>%
  ggplot(aes(x = input_reads, y = total_reads / 1e6, color = type)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.7) +
  scale_color_manual(values = type_colors) +
  labs(title = "Bracken Assigned Reads vs Sequencing Depth",
       x = "Total Input Reads (millions)", y = "Bracken Assigned Reads (millions)", color = NULL) +
  theme_classic()

(p_bracken_vs_depth | p_bracken_vs_qual)

p_bracken_depth <- ggplot(bracken_totals, aes(x = type, y = total_reads / 1e6, fill = type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(aes(color = type), width = 0.15, size = 1.5, alpha = 0.8) +
  scale_fill_manual(values = type_colors) +
  scale_color_manual(values = type_colors) +
  labs(title = "Bracken: Total Species-Assigned Reads by Condition",
       x = NULL, y = "Reads (millions)") +
  theme_classic() +
  theme(legend.position = "none")
