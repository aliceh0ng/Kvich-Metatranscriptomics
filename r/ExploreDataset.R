source("/Users/alicehong/kvich_meta/r/LoadData.R")

# Read depth per sample by condition
p_depth <- ggplot(meta, aes(x = type, y = total_sequences, fill = type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(aes(color = type), width = 0.15, size = 1.5, alpha = 0.8) +
  scale_fill_manual(values = type_colors) +
  scale_color_manual(values = type_colors) +
  labs(title = "Read Depth by Condition", x = NULL, y = "Total Sequences (millions)") +
  theme_classic() + theme(legend.position = "none") + text_size

p_depth2 <- p_depth +
  coord_cartesian(ylim = c(0, 1.5)) +
  labs(title = "Read Depth (< 2M, zoomed)") + text_size

# Sequence quality by condition
p_qual <- ggplot(meta, aes(x = type, y = mean_quality, fill = type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(aes(color = type), width = 0.15, size = 1.5, alpha = 0.8) +
  scale_fill_manual(values = type_colors) +
  scale_color_manual(values = type_colors) +
  labs(title = "Mean Sequence Quality by Condition", x = NULL, y = "Mean Phred Quality Score") +
  theme_classic() + theme(legend.position = "none") + text_size

# % Duplicates by condition
p_dups <- ggplot(meta, aes(x = type, y = percent_dups, fill = type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(aes(color = type), width = 0.15, size = 1.5, alpha = 0.8) +
  scale_fill_manual(values = type_colors) +
  scale_color_manual(values = type_colors) +
  labs(title = "% Duplicate Reads by Condition", x = NULL, y = "% Duplicates") +
  theme_classic() + theme(legend.position = "none") + text_size
