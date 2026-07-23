#!/usr/bin/env Rscript
# simplseq_to_dinemites.R
# ---------------------------------------------------------------------------
# Bridge between SIMPLseq CIGAR output and DINEMITES longitudinal input format.
#
# Reads the wide seqtab_cigar.tsv and samples.csv, performs replicate
# intersection merging (keep only haplotypes present in BOTH replicates),
# converts exact YYYY-MM-DD dates as-is, or YYYY-MM dates using the default day, and outputs
# a tab-delimited file with columns: allele, time, subject, locus.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(tidyr)
})

collection_date_to_full_date <- function(collection_date, day_of_month) {
  raw_values <- trimws(as.character(collection_date))
  parsed <- rep(as.Date(NA), length(raw_values))

  has_full_day <- grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$", raw_values)
  has_month_only <- grepl("^[0-9]{4}-[0-9]{2}$", raw_values)

  parsed[has_full_day] <- as.Date(raw_values[has_full_day], format = "%Y-%m-%d")
  parsed[has_month_only] <- as.Date(
    paste0(raw_values[has_month_only], "-", sprintf("%02d", day_of_month)),
    format = "%Y-%m-%d"
  )

  parsed
}

first_non_empty <- function(values) {
  values <- trimws(as.character(values))
  values <- values[!is.na(values) & nchar(values) > 0]
  if (length(values) == 0) {
    return("")
  }
  values[[1]]
}

clean_numeric <- function(values) {
  suppressWarnings(as.numeric(gsub(",", "", as.character(values))))
}

fill_numeric_covariate <- function(data, source_col, target_col) {
  if (!source_col %in% colnames(data)) {
    return(data)
  }
  raw_values <- clean_numeric(data[[source_col]])
  missing <- is.na(raw_values)
  fallback <- suppressWarnings(stats::median(raw_values[!missing], na.rm = TRUE))
  if (is.na(fallback) || !is.finite(fallback)) {
    fallback <- 0
  }
  data[[target_col]] <- ifelse(missing, fallback, raw_values)
  data[[paste0(target_col, "_missing")]] <- as.integer(missing)
  data
}

add_binary_covariate <- function(data, source_col, target_col, positive_pattern) {
  if (!source_col %in% colnames(data)) {
    return(data)
  }
  raw_values <- tolower(trimws(as.character(data[[source_col]])))
  missing <- is.na(raw_values) | nchar(raw_values) == 0
  data[[target_col]] <- as.integer(!missing & grepl(positive_pattern, raw_values))
  data[[paste0(target_col, "_missing")]] <- as.integer(missing)
  data
}

add_gender_covariate <- function(data) {
  if (!"metadata_gender" %in% colnames(data)) {
    return(data)
  }
  raw_values <- tolower(trimws(as.character(data$metadata_gender)))
  female <- raw_values %in% c("female", "f")
  male <- raw_values %in% c("male", "m")
  missing <- is.na(raw_values) | nchar(raw_values) == 0 | (!female & !male)
  data$covariate_gender <- as.integer(female)
  data$covariate_gender_missing <- as.integer(missing)
  data
}

add_metadata_covariates <- function(data) {
  data <- add_binary_covariate(data, "metadata_season", "covariate_season", "wet|rain")
  data <- fill_numeric_covariate(data, "metadata_age", "covariate_age")
  data <- add_gender_covariate(data)
  data
}

# --- CLI arguments ---
option_list <- list(
  make_option("--cigar", type = "character", help = "Path to seqtab_cigar.tsv"),
  make_option("--samples", type = "character", help = "Path to samples.csv"),
  make_option("--out", type = "character", help = "Output path for DINEMITES input TSV"),
  make_option("--preparation_summary", type = "character", default = NULL,
              help = "Optional JSON audit of filtering and replicate merging"),
  make_option("--min_abundance_pct", type = "double", default = 0.3,
              help = "Minimum allele abundance within each sequencing sample [default %default]"),
  make_option("--abundance_denominator", type = "character", default = "locus",
              help = "Allele abundance denominator: locus or sample [default %default]"),
  make_option("--day_of_month", type = "integer", default = 27,
              help = "Day to assume for YYYY-MM collection dates [default %default]")
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

cat("[DINEMITES/bridge] Reading CIGAR table:", args$cigar, "\n")
cat("[DINEMITES/bridge] Reading samples:", args$samples, "\n")
cat("[DINEMITES/bridge] Allele abundance filter:",
    args$min_abundance_pct, "% of",
    ifelse(args$abundance_denominator == "locus", "sample+locus reads", "total sample reads"),
    "\n")

# --- 1. Read inputs ---
cigar_wide <- read.delim(args$cigar, header = TRUE, sep = "\t",
                         check.names = FALSE, stringsAsFactors = FALSE)

samples <- read.csv(args$samples, header = TRUE, stringsAsFactors = FALSE)

# Normalise column names for samples.csv (handle minor variations)
colnames(samples) <- tolower(trimws(colnames(samples)))
metadata_cols <- grep("^(metadata_|covariate_)", colnames(samples), value = TRUE)

required_sample_cols <- c("sample_id", "participant_id", "collection_date", "replicate")
missing_cols <- setdiff(required_sample_cols, colnames(samples))
if (length(missing_cols) > 0) {
  stop("[DINEMITES/bridge] ERROR: Missing required samples.csv columns: ",
       paste(missing_cols, collapse = ", "))
}

# --- 2. Parse column headers: "LOCUS,CIGAR" to locus + cigar ---
# The first column is "sample"; remaining are "LOCUS,CIGAR" encoded haplotypes
haplotype_cols <- setdiff(colnames(cigar_wide), "sample")

parsed_cols <- data.frame(
  col_name = haplotype_cols,
  locus    = sub(",.*", "", haplotype_cols),
  cigar    = sub("^[^,]*,", "", haplotype_cols),
  stringsAsFactors = FALSE
)

cat("[DINEMITES/bridge]", nrow(parsed_cols), "haplotype columns across",
    length(unique(parsed_cols$locus)), "loci\n")

# --- 3. Pivot wide to long ---
cigar_long <- cigar_wide %>%
  pivot_longer(
    cols      = all_of(haplotype_cols),
    names_to  = "haplotype_col",
    values_to = "reads"
  ) %>%
  left_join(parsed_cols, by = c("haplotype_col" = "col_name")) %>%
  mutate(reads = suppressWarnings(as.numeric(.data$reads))) %>%
  mutate(reads = ifelse(is.na(.data$reads), 0, .data$reads)) %>%
  select(sample_id = sample, locus, cigar, reads)

# --- 4. Join metadata ---
# Match samples.csv on sample_id
sample_type <- if ("sample_type" %in% colnames(samples)) {
  tolower(trimws(samples$sample_type))
} else {
  rep("sample", nrow(samples))
}
control_source <- if ("collection_date_source" %in% colnames(samples)) {
  tolower(trimws(samples$collection_date_source))
} else {
  rep("", nrow(samples))
}

meta_candidates <- samples %>%
  mutate(
    .sample_type_for_filter = sample_type,
    .control_source_for_filter = control_source,
    sample_id = ifelse(is.na(.data$sample_id), "", trimws(as.character(.data$sample_id))),
    participant_id = ifelse(is.na(.data$participant_id), "", trimws(as.character(.data$participant_id))),
    collection_date = ifelse(is.na(.data$collection_date), "", trimws(as.character(.data$collection_date))),
    replicate = ifelse(is.na(.data$replicate), "", trimws(as.character(.data$replicate)))
  ) %>%
  filter(
    .data$.sample_type_for_filter == "sample" &
      nchar(.data$participant_id) > 0 &
      nchar(.data$collection_date) > 0 &
      .data$.control_source_for_filter != "control_excluded_from_metadata"
  )

excluded_controls <- nrow(samples) - nrow(meta_candidates)
if (excluded_controls > 0) {
  cat("[DINEMITES/bridge] Excluded", excluded_controls,
      "sequencing control rows from longitudinal analysis.\n")
}

matched_candidates <- meta_candidates %>%
  semi_join(cigar_long %>% distinct(sample_id), by = "sample_id")

if (nrow(matched_candidates) < nrow(meta_candidates)) {
  cat("[DINEMITES/bridge] WARNING:", nrow(meta_candidates) - nrow(matched_candidates),
      "sample sheet rows are not present in the CIGAR table and will be ignored.\n")
}

if (nrow(matched_candidates) == 0) {
  stop("[DINEMITES/bridge] ERROR: No sample_id values matched between the CIGAR table and samples.csv.")
}

missing_participant <- sum(nchar(matched_candidates$participant_id) == 0)
missing_date <- sum(nchar(matched_candidates$collection_date) == 0)
if (missing_participant > 0) {
  stop("[DINEMITES/bridge] ERROR: DINEMITES requires participant_id for every matched sample. ",
       missing_participant, " matched rows are missing participant_id.")
}
if (missing_date > 0) {
  stop("[DINEMITES/bridge] ERROR: DINEMITES requires collection_date values in YYYY-MM or YYYY-MM-DD format. ",
       missing_date, " matched sample rows are missing collection_date. Add metadata or rescan with a fallback year.")
}

meta <- matched_candidates %>%
  select(sample_id, participant_id, collection_date, replicate,
         dplyr::all_of(metadata_cols))

matched_meta <- meta

visit_replicates <- matched_meta %>%
  distinct(sample_id, participant_id, collection_date) %>%
  count(participant_id, collection_date, name = "n_replicates")

sample_timepoints <- matched_meta %>%
  group_by(participant_id, collection_date) %>%
  summarise(across(dplyr::all_of(metadata_cols), first_non_empty),
            .groups = "drop") %>%
  mutate(
    date_full = collection_date_to_full_date(.data$collection_date, args$day_of_month)
  )

if (any(is.na(sample_timepoints$date_full))) {
  stop("[DINEMITES/bridge] ERROR: Could not parse one or more collection_date values. ",
       "Expected YYYY-MM or YYYY-MM-DD.")
}

sample_timepoints <- sample_timepoints %>%
  mutate(
    time = as.integer(date_full - min(date_full)),
    date_label = format(date_full, "%d %b %Y")
  )
sample_timepoints <- add_metadata_covariates(sample_timepoints)

cigar_meta <- cigar_long %>%
  inner_join(matched_meta, by = "sample_id")

if (nrow(cigar_meta) == 0) {
  stop("[DINEMITES/bridge] ERROR: No samples matched between CIGAR table and samples.csv. ",
       "Check that sample IDs match.")
}

cat("[DINEMITES/bridge]", nrow(cigar_meta), "rows after metadata join (",
    length(unique(cigar_meta$participant_id)), " participants)\n")

# --- 5. Determine presence (reads > 0) ---
cigar_meta <- cigar_meta %>%
  group_by(.data$sample_id, dplyr::across(dplyr::all_of(
    if (args$abundance_denominator == "locus") "locus" else character()
  ))) %>%
  mutate(total_sample_reads = sum(.data$reads, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(
    min_reads_required = ifelse(
      .data$total_sample_reads > 0 & args$min_abundance_pct > 0,
      pmax(1, ceiling(.data$total_sample_reads * (args$min_abundance_pct / 100))),
      1
    ),
    allele_abundance_pct = ifelse(
      .data$total_sample_reads > 0,
      100 * .data$reads / .data$total_sample_reads,
      0
    ),
    present = as.integer(.data$reads > 0 & .data$reads >= min_reads_required)
  )

positive_calls_before_filter <- sum(cigar_meta$reads > 0, na.rm = TRUE)
positive_calls_after_filter <- sum(cigar_meta$present > 0, na.rm = TRUE)

# --- 6. Replicate intersection merging ---
# For each participant + date + locus + cigar combination:
# Keep haplotype ONLY if present in ALL replicates for that time point.
# If only one replicate exists, keep as-is.

merged <- cigar_meta %>%
  group_by(participant_id, collection_date, locus, cigar) %>%
  summarise(
    n_replicates     = n_distinct(sample_id),
    n_present        = n_distinct(sample_id[present > 0]),
    min_reads_required = max(min_reads_required, na.rm = TRUE),
    max_abundance_pct = max(allele_abundance_pct, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # Intersection: keep only if present in ALL replicates
  filter(n_present == n_replicates & n_present > 0)

cat("[DINEMITES/bridge]", nrow(merged),
    "haplotype-timepoints after replicate intersection merging\n")

if (nrow(merged) == 0) {
  cat("[DINEMITES/bridge] WARNING: No haplotypes survived replicate intersection. ",
      "Writing sampled visits as allele=NA.\n", sep = "")
}

# --- 7. Date conversion: exact YYYY-MM-DD, or YYYY-MM plus default day ---
merged <- merged %>%
  left_join(sample_timepoints,
            by = c("participant_id", "collection_date"))

# Report date range
date_range <- range(sample_timepoints$date_full)
cat("[DINEMITES/bridge] Date range:",
    format(date_range[1]), "to", format(date_range[2]),
    "(", max(sample_timepoints$time), "days )\n")

# --- 8. Build DINEMITES input ---
# Allele naming: LOCUS:CIGAR (e.g., "KELT:139C", "SERA8:.")
timepoint_extra_cols <- setdiff(
  colnames(sample_timepoints),
  c("participant_id", "collection_date", "date_full", "time", "date_label")
)

present_input <- merged %>%
  mutate(
    allele  = paste0(locus, ":", cigar),
    subject = participant_id
  ) %>%
  select(allele, time, subject, locus, collection_date, date_full, date_label,
         dplyr::all_of(timepoint_extra_cols))

empty_timepoints <- sample_timepoints %>%
  anti_join(merged %>% distinct(participant_id, collection_date),
            by = c("participant_id", "collection_date")) %>%
  transmute(
    allele = NA_character_,
    time,
    subject = participant_id,
    locus = NA_character_,
    collection_date,
    date_full,
    date_label,
    across(dplyr::all_of(timepoint_extra_cols))
  )

dinemites_input <- bind_rows(present_input, empty_timepoints) %>%
  arrange(subject, time, locus, allele)

# Report summary
n_subjects  <- length(unique(dinemites_input$subject))
n_timepoints <- length(unique(dinemites_input$time))
n_alleles   <- length(unique(dinemites_input$allele[!is.na(dinemites_input$allele)]))
n_loci      <- length(unique(dinemites_input$locus[!is.na(dinemites_input$locus)]))

cat("[DINEMITES/bridge] Output summary:\n")
cat("  Subjects:    ", n_subjects, "\n")
cat("  Time points: ", n_timepoints, "\n")
cat("  Unique alleles:", n_alleles, "\n")
cat("  Loci:        ", n_loci, "\n")
cat("  Total rows:  ", nrow(dinemites_input), "\n")

# --- 9. Write output ---
dir.create(dirname(args$out), showWarnings = FALSE, recursive = TRUE)
write.table(dinemites_input, file = args$out, sep = "\t",
            row.names = FALSE, quote = FALSE)

if (!is.null(args$preparation_summary) && nchar(trimws(args$preparation_summary)) > 0) {
  preparation_summary <- list(
    schema_version = 1,
    source = list(
      table = basename(args$cigar),
      matched_sequencing_samples = n_distinct(matched_meta$sample_id),
      controls_excluded = excluded_controls,
      sample_sheet_rows_not_in_cigar = nrow(meta_candidates) - nrow(matched_candidates)
    ),
    abundance_filter = list(
      threshold_percent = args$min_abundance_pct,
      denominator = args$abundance_denominator,
      positive_calls_before = positive_calls_before_filter,
      positive_calls_after = positive_calls_after_filter
    ),
    replicate_merge = list(
      rule = "intersection",
      visits_with_multiple_replicates = sum(visit_replicates$n_replicates > 1),
      maximum_replicates_per_visit = max(visit_replicates$n_replicates),
      retained_participant_visit_allele_calls = nrow(merged)
    ),
    submitted = list(
      subjects = n_subjects,
      participant_visits = nrow(sample_timepoints),
      unique_time_values = n_timepoints,
      loci = n_loci,
      unique_alleles = n_alleles,
      genotype_rows = nrow(present_input),
      empty_visit_rows = nrow(empty_timepoints),
      total_rows = nrow(dinemites_input)
    )
  )
  dir.create(dirname(args$preparation_summary), showWarnings = FALSE, recursive = TRUE)
  jsonlite::write_json(
    preparation_summary,
    path = args$preparation_summary,
    auto_unbox = TRUE,
    pretty = TRUE,
    na = "null"
  )
  cat("[DINEMITES/bridge] Wrote", args$preparation_summary, "\n")
}

cat("[DINEMITES/bridge] Wrote", args$out, "\n")
cat("[DINEMITES/bridge] Done.\n")
