source("/Users/alicehong/kvich_meta/r/LoadData.R")

# flagstat loaded in LoadData.R

# Summary table: mean metrics per pipeline x condition
flagstat_summary <- flagstat %>%
  group_by(pipeline, type) %>%
  summarise(n                        = n(),
            mean_pct_mapped          = mean(pct_mapped,           na.rm = TRUE),
            mean_pct_properly_paired = mean(pct_properly_paired,  na.rm = TRUE),
            mean_pct_singletons      = mean(pct_singletons,       na.rm = TRUE),
            .groups = "drop")

p_mapped <- ggplot(flagstat, aes(x = type, y = pct_mapped, fill = type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(aes(color = type, shape = total_reads > 1e6), width = 0.15, size = 1.5, alpha = 0.8) +
  geom_jitter(data = filter(flagstat, total_reads > 2e6),
              aes(x = type, y = pct_mapped), shape = 20, size = 2.5,
              fill = NA, color = "black", stroke = .5) +
  scale_fill_manual(values = type_colors) +
  scale_color_manual(values = type_colors) +
  scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 16), guide = "none") +
  facet_wrap(~ pipeline) +
  labs(x = NULL, y = "% Mapped", caption = "Black = samples with >1M total reads") +
  theme_classic() + theme(legend.position = "none") + text_size
