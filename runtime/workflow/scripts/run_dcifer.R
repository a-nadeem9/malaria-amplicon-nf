#!/usr/bin/env Rscript
# run_dcifer.R
# ---------------------------------------------------------------------------
# Runs dcifer on prepared SIMPLseq long-format CIGAR haplotypes.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- gsub("~+~", " ", sub("^--file=", "", script_arg[[1]]), fixed = TRUE)
script_dir <- dirname(normalizePath(script_path))
source(file.path(script_dir, "lib", "analysis_common.R"))

option_list <- list(
  make_option("--input", type = "character", help = "Path to dcifer input TSV"),
  make_option("--outdir", type = "character", default = ".",
              help = "Output directory for dcifer results [default %default]"),
  make_option("--coi_lrank", type = "integer", default = 2,
              help = "Ranked locus allele count used by getCOI [default %default]"),
  make_option("--ibd_grid_nr", type = "integer", default = 1000,
              help = "Grid resolution for ibdDat confidence intervals [default %default]"),
  make_option("--alpha", type = "double", default = 0.05,
              help = "Alpha for raw p-value flag and confidence intervals [default %default]"),
  make_option("--afreq_mode", type = "character", default = "current_run",
              help = "Allele-frequency source; currently only current_run is supported [default %default]")
)
args <- parse_args(OptionParser(option_list = option_list))

if (is.null(args$input)) {
  stop("Required argument: --input")
}
if (is.na(args$coi_lrank) || args$coi_lrank < 1) {
  stop("--coi_lrank must be at least 1.")
}
if (is.na(args$ibd_grid_nr) || args$ibd_grid_nr < 1) {
  stop("--ibd_grid_nr must be at least 1.")
}
if (is.na(args$alpha) || args$alpha <= 0 || args$alpha >= 1) {
  stop("--alpha must be greater than 0 and less than 1.")
}
args$afreq_mode <- tolower(trimws(as.character(args$afreq_mode)))
if (args$afreq_mode != "current_run") {
  stop("--afreq_mode currently supports only current_run.")
}

if (!requireNamespace("dcifer", quietly = TRUE)) {
  stop("[dcifer/run] ERROR: dcifer package not installed. ",
       "Install with: install.packages('dcifer')")
}

library(dcifer)

symmetrize_lower_matrix <- function(mat, diagonal = NA_real_) {
  out <- mat
  out[upper.tri(out)] <- t(out)[upper.tri(out)]
  diag(out) <- diagonal
  out
}

matrix_to_long <- function(mat, value_name) {
  out <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
  colnames(out) <- c("sample_a", "sample_b", value_name)
  out
}

write_matrix <- function(mat, path) {
  table <- data.frame(sample_id = rownames(mat), mat, check.names = FALSE)
  write.table(table, file = path, sep = "\t", row.names = FALSE, quote = FALSE)
}

cat("[dcifer/run] Input:", args$input, "\n")
cat("[dcifer/run] Output dir:", args$outdir, "\n")
cat("[dcifer/run] COI lrank:", args$coi_lrank, "\n")
cat("[dcifer/run] IBD grid:", args$ibd_grid_nr, "\n")

dir.create(args$outdir, showWarnings = FALSE, recursive = TRUE)
plots_dir <- file.path(args$outdir, "dcifer_plots")
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

dcifer_input <- read.delim(args$input, check.names = FALSE, stringsAsFactors = FALSE)
required_cols <- c("sample_id", "locus", "allele")
require_columns(dcifer_input, required_cols, "dcifer/run")
if (nrow(dcifer_input) == 0) {
  stop("[dcifer/run] Input has no allele rows.")
}

if (!"participant_id" %in% colnames(dcifer_input)) {
  dcifer_input$participant_id <- dcifer_input$sample_id
}
if (!"collection_date" %in% colnames(dcifer_input)) {
  dcifer_input$collection_date <- ""
}

metadata <- dcifer_input %>%
  distinct(.data$sample_id, .data$participant_id, .data$collection_date)

cat("[dcifer/run] Formatting data for dcifer...\n")
dsmp <- dcifer::formatDat(dcifer_input, svar = "sample_id", lvar = "locus", avar = "allele")

if (length(dsmp) < 2) {
  stop("[dcifer/run] At least two dcifer samples are required for pairwise relatedness.")
}

available_loci <- length(unique(dcifer_input$locus))
coi_lrank <- min(args$coi_lrank, available_loci)
if (coi_lrank < args$coi_lrank) {
  cat("[dcifer/run] Requested COI lrank", args$coi_lrank,
      "but only", available_loci, "locus/loci are available; using lrank",
      coi_lrank, "\n")
}

cat("[dcifer/run] Estimating COI...\n")
coi <- dcifer::getCOI(dsmp, lrank = coi_lrank)
if (any(is.na(coi))) {
  bad_samples <- names(coi)[is.na(coi)]
  stop("[dcifer/run] COI could not be estimated for: ", paste(bad_samples, collapse = ", "),
       ". Lower --coi_lrank or add more loci.")
}

coi_table <- data.frame(
  sample_id = names(coi),
  coi = as.integer(coi),
  stringsAsFactors = FALSE
) %>%
  left_join(metadata, by = "sample_id") %>%
  select(sample_id, participant_id, collection_date, coi) %>%
  arrange(participant_id, collection_date, sample_id)

coi_path <- file.path(args$outdir, "dcifer_coi.tsv")
write.table(coi_table, file = coi_path, sep = "\t", row.names = FALSE, quote = FALSE)
cat("[dcifer/run] Wrote COI:", coi_path, "\n")

cat("[dcifer/run] Estimating allele frequencies from current run...\n")
afreq <- dcifer::calcAfreq(dsmp, coi, tol = 1e-5)
afreq_table <- lapply(names(afreq), function(locus_name) {
  values <- afreq[[locus_name]]
  data.frame(
    locus = locus_name,
    allele = names(values),
    allele_frequency = as.numeric(values),
    stringsAsFactors = FALSE
  )
}) %>%
  bind_rows() %>%
  arrange(.data$locus, desc(.data$allele_frequency), .data$allele)

afreq_path <- file.path(args$outdir, "dcifer_allele_frequencies.tsv")
write.table(afreq_table, file = afreq_path, sep = "\t", row.names = FALSE, quote = FALSE)
cat("[dcifer/run] Wrote allele frequencies:", afreq_path, "\n")

cat("[dcifer/run] Estimating pairwise relatedness...\n")
dres <- dcifer::ibdDat(
  dsmp,
  coi,
  afreq,
  pval = TRUE,
  confint = TRUE,
  rnull = 0,
  alpha = args$alpha,
  nr = args$ibd_grid_nr
)

samples <- names(dsmp)
pairwise_relatedness <- lapply(seq_along(samples), function(i) {
  lapply(seq_along(samples), function(j) {
    if (i <= j) {
      return(NULL)
    }
    data.frame(
      sample_a = samples[i],
      sample_b = samples[j],
      estimate = dres[i, j, "estimate"],
      p_value = dres[i, j, "p_value"],
      ci_lower = dres[i, j, "CI_lower"],
      ci_upper = dres[i, j, "CI_upper"],
      stringsAsFactors = FALSE
    )
  }) %>%
    bind_rows()
}) %>%
  bind_rows() %>%
  left_join(
    metadata %>%
      select(sample_a = sample_id,
             participant_a = participant_id,
             collection_date_a = collection_date),
    by = "sample_a"
  ) %>%
  left_join(
    metadata %>%
      select(sample_b = sample_id,
             participant_b = participant_id,
             collection_date_b = collection_date),
    by = "sample_b"
  ) %>%
  mutate(
    comparison_type = ifelse(.data$participant_a == .data$participant_b,
                             "within_patient", "between_patient"),
    q_value = p.adjust(.data$p_value, method = "BH"),
    raw_p_le_alpha = .data$p_value <= args$alpha,
    q_le_alpha = .data$q_value <= args$alpha
  ) %>%
  arrange(desc(.data$estimate), .data$p_value, .data$sample_a, .data$sample_b)

pairs_path <- file.path(args$outdir, "dcifer_pairwise_relatedness.tsv")
write.table(pairwise_relatedness, file = pairs_path, sep = "\t", row.names = FALSE, quote = FALSE)
cat("[dcifer/run] Wrote pairwise relatedness:", pairs_path, "\n")

estimate_mat <- symmetrize_lower_matrix(dres[, , "estimate"], diagonal = 1)
pvalue_mat <- symmetrize_lower_matrix(dres[, , "p_value"], diagonal = NA_real_)

estimate_matrix_path <- file.path(args$outdir, "dcifer_relatedness_matrix.tsv")
pvalue_matrix_path <- file.path(args$outdir, "dcifer_pvalue_matrix.tsv")
write_matrix(estimate_mat, estimate_matrix_path)
write_matrix(pvalue_mat, pvalue_matrix_path)

estimate_long <- matrix_to_long(estimate_mat, "estimate") %>%
  mutate(
    sample_a = factor(.data$sample_a, levels = samples),
    sample_b = factor(.data$sample_b, levels = rev(samples))
  )

heatmap_inches <- min(24, max(12, length(samples) * 0.065))
axis_text_size <- if (length(samples) > 150) 5 else if (length(samples) > 75) 6 else 8

relatedness_heatmap <- ggplot(estimate_long, aes(x = .data$sample_a, y = .data$sample_b, fill = .data$estimate)) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_gradient(low = analysis_palette$teal_pale, high = analysis_palette$teal_dark, limits = c(0, 1), na.value = analysis_palette$neutral) +
  coord_fixed() +
  labs(
    title = "dcifer pairwise relatedness estimates",
    x = "Sample",
    y = "Sample",
    fill = "Relatedness"
  ) +
  analysis_theme(base_size = 13) +
  theme(axis.text = element_text(size = axis_text_size),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 16),
        panel.grid = element_blank())

relatedness_plot_path <- file.path(plots_dir, "dcifer_relatedness_heatmap.png")
ggsave(relatedness_plot_path, plot = relatedness_heatmap, width = heatmap_inches, height = heatmap_inches, dpi = 220, limitsize = FALSE)

pvalue_long <- matrix_to_long(pvalue_mat, "p_value") %>%
  mutate(
    sample_a = factor(.data$sample_a, levels = samples),
    sample_b = factor(.data$sample_b, levels = rev(samples)),
    neg_log10_p = ifelse(is.na(.data$p_value), NA_real_, -log10(pmax(.data$p_value, 1e-16)))
  )

pvalue_heatmap <- ggplot(pvalue_long, aes(x = .data$sample_a, y = .data$sample_b, fill = .data$neg_log10_p)) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_gradient(low = analysis_palette$pink_pale, high = analysis_palette$pink, na.value = analysis_palette$neutral) +
  coord_fixed() +
  labs(
    title = "dcifer raw p-values for H0: relatedness = 0",
    x = "Sample",
    y = "Sample",
    fill = "-log10(p)"
  ) +
  analysis_theme(base_size = 13) +
  theme(axis.text = element_text(size = axis_text_size),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 16),
        panel.grid = element_blank())

pvalue_plot_path <- file.path(plots_dir, "dcifer_pvalue_heatmap.png")
ggsave(pvalue_plot_path, plot = pvalue_heatmap, width = heatmap_inches, height = heatmap_inches, dpi = 220, limitsize = FALSE)
cat("[dcifer/run] Wrote plots:", plots_dir, "\n")

summary_path <- file.path(args$outdir, "dcifer_summary.json")
max_relatedness <- if (nrow(pairwise_relatedness) > 0) {
  max(pairwise_relatedness$estimate, na.rm = TRUE)
} else {
  NA_real_
}
raw_p_count <- sum(pairwise_relatedness$raw_p_le_alpha, na.rm = TRUE)
write_json_file(list(
  status = "complete",
  afreq_mode = args$afreq_mode,
  coi_lrank = coi_lrank,
  ibd_grid_nr = args$ibd_grid_nr,
  alpha = args$alpha,
  samples = length(samples),
  loci = available_loci,
  pairs = nrow(pairwise_relatedness),
  max_relatedness = if (is.na(max_relatedness)) NA_real_ else round(max_relatedness, 6),
  raw_p_le_alpha = raw_p_count,
  q_le_alpha = sum(pairwise_relatedness$q_le_alpha, na.rm = TRUE),
  adequacy = if (available_loci < 2) "not_recommended_single_locus" else "multi_locus",
  caveat = "Current-run allele frequencies and raw p-values are exploratory; interpret relatedness only with adequate representative background data."
), summary_path)
cat("[dcifer/run] Wrote summary:", summary_path, "\n")
cat("[dcifer/run] Done.\n")
