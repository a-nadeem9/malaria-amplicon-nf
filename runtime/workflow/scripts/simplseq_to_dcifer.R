#!/usr/bin/env Rscript
# simplseq_to_dcifer.R
# ---------------------------------------------------------------------------
# Bridge between SIMPLseq CIGAR output and dcifer long-format input.
#
# Reads the wide seqtab_cigar.tsv and samples.csv, applies an allele abundance
# filter per sequencing sample, merges technical replicates by intersection,
# and writes one row per retained patient/date/locus/allele.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(tidyr)
})

option_list <- list(
  make_option("--cigar", type = "character", help = "Path to seqtab_cigar.tsv"),
  make_option("--samples", type = "character", help = "Path to samples.csv"),
  make_option("--out", type = "character", help = "Output path for dcifer input TSV"),
  make_option("--filter_summary", type = "character", default = "dcifer_filter_summary.tsv",
              help = "Output path for per-replicate filter summary [default %default]"),
  make_option("--replicate_summary", type = "character", default = "dcifer_replicate_summary.tsv",
              help = "Output path for replicate-intersection summary [default %default]"),
  make_option("--min_abundance_pct", type = "double", default = 0.3,
              help = "Minimum allele abundance within each sequencing sample [default %default]"),
  make_option("--abundance_denominator", type = "character", default = "locus",
              help = "Allele abundance denominator: locus or sample [default %default]")
)
args <- parse_args(OptionParser(option_list = option_list))

if (is.null(args$cigar) || is.null(args$samples) || is.null(args$out)) {
  stop("Required arguments: --cigar, --samples, --out")
}
if (is.na(args$min_abundance_pct) || args$min_abundance_pct < 0 || args$min_abundance_pct > 100) {
  stop("--min_abundance_pct must be between 0 and 100. Got: ", args$min_abundance_pct)
}
args$abundance_denominator <- tolower(trimws(as.character(args$abundance_denominator)))
if (!args$abundance_denominator %in% c("locus", "sample")) {
  stop("--abundance_denominator must be either locus or sample. Got: ",
       args$abundance_denominator)
}

cat("[dcifer/bridge] Reading CIGAR table:", args$cigar, "\n")
cat("[dcifer/bridge] Reading samples:", args$samples, "\n")
cat("[dcifer/bridge] Allele abundance filter:",
    args$min_abundance_pct, "% of",
    ifelse(args$abundance_denominator == "locus", "sample+locus reads", "total sample reads"),
    "\n")

cigar_wide <- read.delim(args$cigar, header = TRUE, sep = "\t",
                         check.names = FALSE, stringsAsFactors = FALSE)
samples <- read.csv(args$samples, header = TRUE, stringsAsFactors = FALSE)
colnames(samples) <- tolower(trimws(colnames(samples)))

required_sample_cols <- c("sample_id", "participant_id", "collection_date", "replicate")
missing_cols <- setdiff(required_sample_cols, colnames(samples))
if (length(missing_cols) > 0) {
  stop("[dcifer/bridge] Missing required samples.csv columns: ",
       paste(missing_cols, collapse = ", "))
}

haplotype_cols <- setdiff(colnames(cigar_wide), "sample")
if (length(haplotype_cols) == 0) {
  stop("[dcifer/bridge] No CIGAR haplotype columns found.")
}

parsed_cols <- data.frame(
  col_name = haplotype_cols,
  locus = sub(",.*", "", haplotype_cols),
  allele = sub("^[^,]*,", "", haplotype_cols),
  stringsAsFactors = FALSE
)

cigar_long <- cigar_wide %>%
  pivot_longer(
    cols = all_of(haplotype_cols),
    names_to = "haplotype_col",
    values_to = "reads"
  ) %>%
  left_join(parsed_cols, by = c("haplotype_col" = "col_name")) %>%
  mutate(reads = suppressWarnings(as.numeric(.data$reads))) %>%
  mutate(reads = ifelse(is.na(.data$reads), 0, .data$reads)) %>%
  select(sample_id = sample, locus, allele, reads)

sample_type <- if ("sample_type" %in% colnames(samples)) {
  tolower(trimws(samples$sample_type))
} else {
  rep("sample", nrow(samples))
}

meta <- samples %>%
  mutate(
    .sample_type_for_filter = sample_type,
    sample_id = ifelse(is.na(.data$sample_id), "", trimws(as.character(.data$sample_id))),
    participant_id = ifelse(is.na(.data$participant_id), "", trimws(as.character(.data$participant_id))),
    collection_date = ifelse(is.na(.data$collection_date), "", trimws(as.character(.data$collection_date))),
    replicate = ifelse(is.na(.data$replicate), "", trimws(as.character(.data$replicate)))
  ) %>%
  filter(
    .data$.sample_type_for_filter == "sample" &
      nchar(.data$participant_id) > 0 &
      nchar(.data$collection_date) > 0
  ) %>%
  mutate(
    participant_id = ifelse(nchar(.data$participant_id) > 0, .data$participant_id, .data$sample_id),
    replicate = ifelse(nchar(.data$replicate) > 0, .data$replicate, "sample"),
    analysis_sample_id = ifelse(
      nchar(.data$collection_date) > 0,
      paste(.data$participant_id, .data$collection_date, sep = "_"),
      .data$sample_id
    )
  ) %>%
  select(sample_id, analysis_sample_id, participant_id, collection_date, replicate)

matched_meta <- meta %>%
  semi_join(cigar_long %>% distinct(sample_id), by = "sample_id")

if (nrow(matched_meta) < nrow(meta)) {
  cat("[dcifer/bridge] WARNING:", nrow(meta) - nrow(matched_meta),
      "sample sheet rows are not present in the CIGAR table and will be ignored.\n")
}
if (any(nchar(matched_meta$collection_date) == 0)) {
  cat("[dcifer/bridge] WARNING:",
      sum(nchar(matched_meta$collection_date) == 0),
      "matched sample rows have no collection_date; those rows will be analyzed as separate sequencing samples.\n")
}

cigar_meta <- cigar_long %>%
  inner_join(matched_meta, by = "sample_id")

if (nrow(cigar_meta) == 0) {
  stop("[dcifer/bridge] No sample_id values matched between the CIGAR table and samples.csv.")
}

if (args$abundance_denominator == "locus") {
  cigar_meta <- cigar_meta %>%
    group_by(.data$sample_id, .data$locus) %>%
    mutate(total_reads_for_filter = sum(.data$reads, na.rm = TRUE)) %>%
    ungroup()
} else {
  cigar_meta <- cigar_meta %>%
    group_by(.data$sample_id) %>%
    mutate(total_reads_for_filter = sum(.data$reads, na.rm = TRUE)) %>%
    ungroup()
}

cigar_meta <- cigar_meta %>%
  mutate(
    min_reads_required = ifelse(
      .data$total_reads_for_filter > 0 & args$min_abundance_pct > 0,
      pmax(1, ceiling(.data$total_reads_for_filter * (args$min_abundance_pct / 100))),
      1
    ),
    allele_abundance_pct = ifelse(
      .data$total_reads_for_filter > 0,
      100 * .data$reads / .data$total_reads_for_filter,
      0
    ),
    present = as.integer(.data$reads > 0 & .data$reads >= .data$min_reads_required)
  )

filter_summary <- cigar_meta %>%
  group_by(.data$analysis_sample_id, .data$participant_id, .data$collection_date, .data$replicate) %>%
  summarise(
    total_nonzero_alleles = sum(.data$reads > 0),
    passing_filter = sum(.data$present > 0),
    .groups = "drop"
  ) %>%
  mutate(removed_by_threshold = .data$total_nonzero_alleles - .data$passing_filter) %>%
  arrange(.data$participant_id, .data$collection_date, .data$replicate)

replicate_merged <- cigar_meta %>%
  group_by(.data$analysis_sample_id, .data$participant_id, .data$collection_date, .data$locus, .data$allele) %>%
  summarise(
    n_replicates = n(),
    n_present = sum(.data$present),
    read_count = max(.data$reads, na.rm = TRUE),
    min_reads_required = max(.data$min_reads_required, na.rm = TRUE),
    max_abundance_pct = max(.data$allele_abundance_pct, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(.data$n_present == .data$n_replicates & .data$n_present > 0)

replicate_summary <- replicate_merged %>%
  group_by(.data$analysis_sample_id, .data$participant_id, .data$collection_date, .data$locus) %>%
  summarise(
    alleles_after_intersection = n(),
    .groups = "drop"
  ) %>%
  arrange(.data$participant_id, .data$collection_date, .data$locus)

dcifer_input <- replicate_merged %>%
  mutate(sample_id = .data$analysis_sample_id) %>%
  select(
    sample_id,
    participant_id,
    collection_date,
    locus,
    allele,
    read_count,
    max_abundance_pct,
    n_replicates,
    n_present,
    min_reads_required
  ) %>%
  arrange(.data$sample_id, .data$locus, .data$allele)

dir.create(dirname(args$out), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(args$filter_summary), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(args$replicate_summary), showWarnings = FALSE, recursive = TRUE)

write.table(dcifer_input, file = args$out, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(filter_summary, file = args$filter_summary, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(replicate_summary, file = args$replicate_summary, sep = "\t", row.names = FALSE, quote = FALSE)

cat("[dcifer/bridge] Output rows:", nrow(dcifer_input), "\n")
cat("[dcifer/bridge] dcifer samples:", length(unique(dcifer_input$sample_id)), "\n")
cat("[dcifer/bridge] Wrote", args$out, "\n")
cat("[dcifer/bridge] Done.\n")
