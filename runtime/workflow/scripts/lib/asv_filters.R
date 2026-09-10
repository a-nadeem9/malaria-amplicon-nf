# Reference filtering, independent of alignment and sequencing dependencies.
apply_reference_filters <- function(asvdf, variant_counts, reference_lengths,
                                    strain = "3D7", length_ratio_min = 0.90) {
  cutoff <- suppressWarnings(as.numeric(length_ratio_min))
  if (length(cutoff) != 1L || !is.finite(cutoff) || cutoff <= 0 || cutoff > 1) {
    stop("--indel_filter must be a length-ratio cutoff greater than 0 and at most 1 (default: 0.90)")
  }

  required <- c("id", "SNP")
  missing <- setdiff(required, colnames(variant_counts))
  if (length(missing) > 0) {
    stop(paste("SNP filter table is missing columns:", paste(missing, collapse = ", ")))
  }
  variant_counts$id <- trimws(as.character(variant_counts$id))
  variant_counts$SNP <- suppressWarnings(as.numeric(variant_counts$SNP))
  if (anyNA(variant_counts$id) || any(variant_counts$id == "") ||
      anyDuplicated(variant_counts$id) || any(!is.finite(variant_counts$SNP)) ||
      any(variant_counts$SNP < 0)) {
    stop("SNP filter table IDs must be unique and thresholds must be finite non-negative values")
  }
  if (is.null(names(reference_lengths)) || anyNA(names(reference_lengths)) ||
      any(names(reference_lengths) == "") || anyDuplicated(names(reference_lengths)) ||
      any(!is.finite(reference_lengths)) || any(reference_lengths <= 0)) {
    stop("Reference lengths must be positive and have unique, non-empty reference IDs")
  }
  if (any(!is.finite(asvdf$haplength)) || any(asvdf$haplength <= 0)) {
    stop("ASV lengths must be finite and positive")
  }

  reference_column <- paste0("refid_", strain)
  snv_distance_column <- paste0("snv_dist_from_", strain)
  tie_column <- paste0("best_reference_ties_", strain)
  asvdf[[paste0("best_refid_", strain)]] <- as.character(asvdf[[reference_column]])
  asvdf$snv_threshold <- NA_real_
  asvdf$reference_length <- NA_real_
  asvdf$indel_length_ratio <- NA_real_
  asvdf$indel_ratio_min <- rep(cutoff, nrow(asvdf))
  asvdf$indel_ratio_max <- rep(1 / cutoff, nrow(asvdf))
  asvdf$indel_filter_method <- rep("broad_length_ratio", nrow(asvdf))
  asvdf$snv_filter <- "FAIL"
  asvdf$indel_filter <- "FAIL"
  asvdf$mapping_status <- "UNMAPPED"

  for (i in seq_len(nrow(asvdf))) {
    refid <- as.character(asvdf[[reference_column]][i])
    threshold_index <- match(refid, variant_counts$id)
    if (is.na(threshold_index)) {
      asvdf$mapping_status[i] <- "MISSING_THRESHOLD"
      next
    }
    reference_index <- match(refid, names(reference_lengths))
    if (is.na(reference_index)) {
      stop(paste("No reference length available for aligned target:", refid))
    }
    asvdf$snv_threshold[i] <- variant_counts$SNP[threshold_index]
    asvdf$reference_length[i] <- reference_lengths[reference_index]
    ratio <- asvdf$haplength[i] / asvdf$reference_length[i]
    asvdf$indel_length_ratio[i] <- ratio
    asvdf$snv_filter[i] <- if (asvdf[[snv_distance_column]][i] <= asvdf$snv_threshold[i]) "PASS" else "FAIL"

    # Broad postProc_dada2.R, commit 6183028737d2c69b619c57b76ef62aa4160ac675:
    # reject lprop < ind OR lprop > 1/ind. This uses net sequence lengths,
    # not aligned gap counts or the reference table's unused INDEL column.
    asvdf$indel_filter[i] <- if (ratio < cutoff || ratio > 1 / cutoff) "FAIL" else "PASS"
    if (asvdf[[tie_column]][i] > 1) {
      asvdf$mapping_status[i] <- "AMBIGUOUS"
    } else if (asvdf$snv_filter[i] == "PASS" && asvdf$indel_filter[i] == "PASS") {
      asvdf$mapping_status[i] <- "MAPPED"
    }
  }
  asvdf[[reference_column]][asvdf$mapping_status != "MAPPED"] <- NA_character_
  asvdf
}
