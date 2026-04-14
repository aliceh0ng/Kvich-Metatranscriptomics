# Script to look at humann3 mapping results

library(readr)
library(dplyr)
library(tidyr)

# Read table
humann_sum <- read_table(
  "/Users/alicehong/kvich_meta/kneaddata/functional_profiling/humann3_alignment_summary.txt",
  comment = "#",
  col_types = cols()
)

# Add % MPA found in ChocoPhlAn
humann_sum <- humann_sum %>%
  mutate(pct_mpa_in_choco = (Choco_found / MPA_species) * 100)

# Summarize each metric across samples
metric_summary <- humann_sum %>%
  summarise(
    `% MPA in Choco` = list(c(
      mean = mean(pct_mpa_in_choco, na.rm = TRUE),
      sd = sd(pct_mpa_in_choco, na.rm = TRUE),
      median = median(pct_mpa_in_choco, na.rm = TRUE),
      min = min(pct_mpa_in_choco, na.rm = TRUE),
      max = max(pct_mpa_in_choco, na.rm = TRUE),
      IQR = IQR(pct_mpa_in_choco, na.rm = TRUE)
    )),
    
    `% Unaligned nucleotide` = list(c(
      mean = mean(`Unalign%_nucl`, na.rm = TRUE),
      sd = sd(`Unalign%_nucl`, na.rm = TRUE),
      median = median(`Unalign%_nucl`, na.rm = TRUE),
      min = min(`Unalign%_nucl`, na.rm = TRUE),
      max = max(`Unalign%_nucl`, na.rm = TRUE),
      IQR = IQR(`Unalign%_nucl`, na.rm = TRUE)
    )),
    
    `% Unaligned translated` = list(c(
      mean = mean(`Unalign%_transl`, na.rm = TRUE),
      sd = sd(`Unalign%_transl`, na.rm = TRUE),
      median = median(`Unalign%_transl`, na.rm = TRUE),
      min = min(`Unalign%_transl`, na.rm = TRUE),
      max = max(`Unalign%_transl`, na.rm = TRUE),
      IQR = IQR(`Unalign%_transl`, na.rm = TRUE)
    ))
  ) %>%
  pivot_longer(everything(), names_to = "metric", values_to = "summary") %>%
  unnest_wider(summary)

metric_summary