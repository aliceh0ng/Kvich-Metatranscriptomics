# Functions used throughout the project

# Sample names not consistent throughotu, La vs Lasse, has no meaning
clean_sample_ids <- function(x) {
  stringr::str_replace_all(x, "Lasse", "La")
}

# Load a featureCounts file: skip comment, clean sample names from BAM paths, keep Geneid + Chr + counts
load_featurecounts <- function(path) {
  df <- read.delim(path, comment.char = "#", check.names = FALSE)
  bam_cols    <- colnames(df)[-(1:6)]
  clean_names <- clean_sample_ids(stringr::str_extract(basename(bam_cols), "^[^.]+"))
  colnames(df) <- c(colnames(df)[1:6], clean_names)
  df %>% dplyr::select(Geneid, Chr, dplyr::all_of(clean_names))
}

# Parse a single samtools flagstat file; returns sample name + read1 count (= number of pairs)
parse_flagstat <- function(file) {
  lines <- readLines(file)
  data.frame(
    sample          = clean_sample_ids(tools::file_path_sans_ext(basename(file))),
    total_fragments = as.numeric(stringr::str_extract(grep("read1", lines, value = TRUE)[1], "^\\d+"))
  )
}

# Sort samples in a long dataframe by ascending total reads_M (for stacked bar plots)
sort_by_total <- function(df) {
  order <- df %>%
    dplyr::group_by(sample) %>%
    dplyr::summarise(t = sum(reads_M), .groups = "drop") %>%
    dplyr::arrange(t) %>%
    dplyr::pull(sample)
  dplyr::mutate(df, sample = factor(sample, levels = order))
}

# Color scheme for CRC vs Paired vs Healthy
type_colors <- c("CRC" = "#c51b7d", "Paired" = "#f1b6da", "Healthy" = "#7fbc41")

# Text size
text_size <- theme(
  axis.title = element_text(size = 7),
  axis.text = element_text(size = 7),
  plot.title = element_text(size = 7),
  strip.text = element_text(size = 6),
  legend.text=element_text(size=7)
  )