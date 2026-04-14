source("/Users/alicehong/kvich_meta/r/LoadData.R")

# Stacked bar of read fate across samples
kd_long <- kd %>%
  dplyr::select(Sample, `final pair1`, reads_lost_trim, reads_lost_decon) %>%
  pivot_longer(-Sample, names_to = "stage", values_to = "reads") %>%
  mutate(stage = factor(stage,
                        levels = c("reads_lost_decon", "reads_lost_trim", "final pair1"),
                        labels = c("Removed (human)", "Removed (QC trim)", "Retained")))

p_bar <- ggplot(kd_long, aes(x = Sample, y = reads / 1e6, fill = stage)) +
  geom_col(width = 1) +
  scale_fill_manual(values = c("Retained" = "#ccebc5", "Removed (QC trim)" = "#b3cde3",
                                "Removed (human)" = "#e7298a"), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Sample (sorted by final read count)", y = "Reads (M)") +
  theme_classic() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position = "right") + text_size

p_bar_zoom <- p_bar +
  coord_cartesian(ylim = c(0, 0.6)) +
  labs(x = NULL, y = "Reads (M)", title = "Zoomed into low depth reads") +
  theme(legend.position = "none") + text_size

# Per-sample read loss summary
kd_loss <- kd %>%
  mutate(total_reads    = `raw pair1`,
         pct_lost_trim  = reads_lost_trim  / total_reads * 100,
         pct_lost_decon = reads_lost_decon / total_reads * 100,
         pct_lost_total = (reads_lost_trim + reads_lost_decon) / total_reads * 100,
         pct_retained   = `final pair1`    / total_reads * 100) %>%
  dplyr::select(Sample, total_reads, reads_lost_trim, pct_lost_trim,
                reads_lost_decon, pct_lost_decon, pct_lost_total, pct_retained)

kd_loss_summary <- kd_loss %>%
  summarise(mean_pct_lost_trim  = mean(pct_lost_trim),
            mean_pct_lost_decon = mean(pct_lost_decon),
            mean_pct_lost_total = mean(pct_lost_total),
            mean_pct_retained   = mean(pct_retained))
