# Script to run Masalin

source("r/LoadLibs.R")
source("r/HelperFunctions.R")
source("r/LoadHumann3.R")

################
# Run Masalin

# Convert for MaAsLin, samples as rows, pathways as columns
feat <- path_mat_filt %>% as.data.frame() %>% t() %>% as.data.frame()
meta_maaslin <- meta_h3 %>% dplyr::select(sample, type, patient_ID) %>% column_to_rownames("sample")
feat <- feat[rownames(meta_maaslin), ]

rownames(feat) <- clean_sample_ids(rownames(feat))

# Metadata for MaAsLin
meta_maaslin <- meta_h3 %>%
  mutate(sample = clean_sample_ids(sample)) %>%
  dplyr::select(sample, type, patient_ID) %>%
  distinct() %>% column_to_rownames("sample")

# Keep only overlapping samples in the same order
common_samples <- intersect(rownames(feat), rownames(meta_maaslin))
feat <- feat[common_samples, , drop = FALSE]
meta_maaslin <- meta_maaslin[common_samples, , drop = FALSE]

# 1) CRC vs Paired (primary)
meta_cp <- meta_maaslin %>%
  rownames_to_column("sample") %>%
  filter(type %in% c("CRC", "Paired")) %>%
  mutate(type = factor(type, levels = c("Paired", "CRC"))) %>%
  column_to_rownames("sample")

feat_cp <- feat[rownames(meta_cp), , drop = FALSE]

# # UNCOMMENT TO RUN!
# fit_crc_vs_paired <- Maaslin2(
#   input_data      = feat_cp,
#   input_metadata  = meta_cp,
#   output          = "maaslin_crc_vs_paired",
#   fixed_effects   = "type",
#   random_effects  = "patient_ID",
#   normalization   = "NONE",
#   transform       = "LOG",
#   analysis_method = "LM")

# 2) CRC vs Healthy
meta_ch <- meta_maaslin %>%
  rownames_to_column("sample") %>%
  filter(type %in% c("CRC", "Healthy")) %>%
  mutate(type = factor(type, levels = c("Healthy", "CRC"))) %>%
  column_to_rownames("sample")

feat_ch <- feat[rownames(meta_ch), , drop = FALSE]

# # UNCOMMENT TO RUN!
# fit_crc_vs_healthy <- Maaslin2(
#   input_data      = feat_ch,
#   input_metadata  = meta_ch,
#   output          = "maaslin_crc_vs_healthy",
#   fixed_effects   = "type",
#   normalization   = "NONE",
#   transform       = "LOG",
#   analysis_method = "LM")

# 3) CRC vs all non-CRC
meta_comb <- meta_maaslin %>%
  rownames_to_column("sample") %>%
  mutate(group = ifelse(type == "CRC", "CRC", "NonCRC"),
         group = factor(group, levels = c("NonCRC", "CRC"))) %>%
  column_to_rownames("sample")

feat_comb <- feat[rownames(meta_comb), , drop = FALSE]

# # UNCOMMENT TO RUN!
# fit_crc_vs_noncrc <- Maaslin2(
#   input_data      = feat_comb,
#   input_metadata  = meta_comb,
#   output          = "maaslin_crc_vs_noncrc",
#   fixed_effects   = "group",
#   normalization   = "NONE",
#   transform       = "LOG",
#   analysis_method = "LM")