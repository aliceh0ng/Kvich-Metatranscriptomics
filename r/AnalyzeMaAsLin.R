
source("r/LoadLibs.R")
source("r/HelperFunctions.R")
source("r/LoadHumann3.R")
source("r/RunMaAsLin.R")

########################
### Plot MaAsLin2 results

res_cp   <- read.delim("kneaddata/maaslin_output/maaslin_crc_vs_paired/all_results.tsv")
res_ch   <- read.delim("kneaddata/maaslin_output/maaslin_crc_vs_healthy/all_results.tsv")
res_comb <- read.delim("kneaddata/maaslin_output/maaslin_crc_vs_noncrc/all_results.tsv")

top_cp <- res_cp %>%
  arrange(qval) %>%
  mutate(model = "CRC vs Paired") %>%
  dplyr::select(model, feature, value, pval, qval)

top_ch <- res_ch %>%
  arrange(qval) %>%
  mutate(model = "CRC vs Healthy") %>%
  dplyr::select(model, feature, value, pval, qval)

top_comb <- res_comb %>%
  arrange(qval) %>%
  mutate(model = "CRC vs NonCRC") %>%
  dplyr::select(model, feature, value, pval, qval)

# Top 10 CRC-associated pathways from paired comparison (controls for patient effects)
top_de <- res_cp %>%
  arrange(qval) %>%
  slice_head(n = 10)

maaslin_table <- top_de %>%
  transmute(
    Pathway    = feature %>% str_replace_all("\\.", " ") %>% str_trunc(100),
    Coefficient = round(coef, 3),
    `p-value`  = signif(pval, 3),
    `q-value`  = signif(qval, 3)
  )

### Slope plot: top 6 CRC-associated pathways (paired normal -> CRC)

# MaAsLin2 sanitizes pathway names with make.names(); map back to original HUMANn3 names
orig_features <- colnames(feat)
feature_map   <- tibble(
  pathway_original = orig_features,
  pathway_maaslin  = make.names(orig_features, unique = TRUE)
)

top6_lookup <- feature_map %>%
  filter(pathway_maaslin %in% (res_cp %>% arrange(qval) %>% slice_head(n = 6) %>% pull(feature)))

slope_df <- path_long %>%
  filter(pathway %in% top6_lookup$pathway_original,
         type %in% c("CRC", "Paired")) %>%
  mutate(
    label = str_remove(pathway, "^[^:]+: ") %>% str_trunc(75),
    type  = factor(type, levels = c("Paired", "CRC"))
  )

p_slope <- slope_df %>%
  ggplot(aes(x = type, y = cpm, group = patient_ID)) +
  geom_line(alpha = 0.4) +
  geom_point(size = 1.5) +
  facet_wrap(~ label, scales = "free_y", ncol = 2) +
  labs(x = NULL, y = "Pathway abundance (CPM)") +
  theme_classic() +
  theme(strip.text = element_text(size = 6), axis.text = element_text(size = 6))
