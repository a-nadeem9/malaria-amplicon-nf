#!/usr/bin/env Rscript

library(argparse)

parser <- ArgumentParser()
parser$add_argument("-s", "--seqtab", 
                    help="Path to input")
parser$add_argument("--samples",
                    help="Sample sheet used to distinguish biological samples from controls")
parser$add_argument("-ref", "--reference",
                    help="Path to reference fasta sequences")
parser$add_argument("-ref2", "--reference2")
parser$add_argument("-noref","--no_reference", action='store_true', help="specify no reference provided")
parser$add_argument("-b","--bimera", help="ASV File with identifed bimeras")
parser$add_argument("-o", "--output", 
                    help="Path to output file")
parser$add_argument("--fasta", action='store_true', help="Write ASV sequences separately into fasta file")
parser$add_argument("-snv", "--snv_filter", 
                    help="Path to file for filtering ASVs based on edit distance")
parser$add_argument("--indel_filter", 
                    help="Optional absolute INDEL-distance override. By default, use the INDEL column in --snv_filter")
parser$add_argument("--strain", default="3D7", help="Name of Specific strain to map to. Defaults to 3D7")
parser$add_argument("--strain2", help="Name of second strain if mapping to 2 different strains")
parser$add_argument("--parallel", action='store_true', help="Enable parallel processing")

args <- parser$parse_args()

# Required packages
library(limma)
library(data.table)
library(stringr)
library(seqinr)
library(parallel)
library(doMC)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- gsub("~+~", " ", sub("^--file=", "", script_arg[[1]]), fixed = TRUE)
script_dir <- dirname(normalizePath(script_path))
source(file.path(script_dir, "lib", "workflow_common.R"))

get_alignment_fun <- function(fun_name) {
  for (namespace in c("pwalign", "Biostrings", "BiocGenerics")) {
    if (requireNamespace(namespace, quietly = TRUE) &&
        exists(fun_name, envir = asNamespace(namespace), inherits = FALSE)) {
      return(get(fun_name, envir = asNamespace(namespace)))
    }
  }
  stop(paste("Alignment function not found:", fun_name))
}

pw_nucleotideSubstitutionMatrix <- get_alignment_fun("nucleotideSubstitutionMatrix")
pw_pairwiseAlignment <- get_alignment_fun("pairwiseAlignment")
pw_score <- get_alignment_fun("score")
pw_alignedPattern <- get_alignment_fun("alignedPattern")
pw_alignedSubject <- get_alignment_fun("alignedSubject")

# Pairwise Alignment
alignment_row <- function(aln, hapid)
{
  scores <- as.numeric(pw_score(aln))
  if (length(scores) == 0 || all(is.na(scores))) {
    stop(paste("No usable reference alignment score for", hapid))
  }
  best_score <- max(scores, na.rm = TRUE)
  best_indices <- which(scores == best_score)
  num <- best_indices[[1]]
  remaining <- scores[-num]
  second_score <- if (length(remaining) == 0 || all(is.na(remaining))) NA_real_ else max(remaining, na.rm = TRUE)
  patt <- c(pw_alignedPattern(aln[num]), pw_alignedSubject(aln[num]))
  dist <- as.numeric(adist(as.character(patt)[1], as.character(patt)[2]))
  ind <- as.numeric(sum(str_count(as.character(patt), "-")))
  data.frame(
    hapid = hapid,
    hapseq = as.character(patt)[2],
    refseq = as.character(patt)[1],
    refid = names(patt)[1],
    aln_score = best_score,
    second_aln_score = second_score,
    alignment_score_margin = if (is.na(second_score)) NA_real_ else best_score - second_score,
    best_reference_ties = length(best_indices),
    snv_dist = (dist - ind),
    indel_dist = ind
  )
}

seq_align <- function(seqs_df, path_to_ref, overlap = TRUE, parallel = TRUE)
{
  library(Biostrings)
  align_df <- data.frame()
  if (file.exists(path_to_ref)) {
    ref <- read.fasta(path_to_ref)
    ref_str <- toupper(sapply(ref, c2s))
  } else {
    stop(paste("File",path_to_ref,"not found!"))
  }
  
  sigma <- pw_nucleotideSubstitutionMatrix(match = 2, mismatch = -1, baseOnly = FALSE)
  if (!overlap) {
    split <- paste0(rep("N",10), collapse = "")
    seq_split <- strsplit2(x = seqs_df[,1], split = split)
    seq_all <- data.frame(sequence=paste0(seq_split[,1],seq_split[,2]), hapid = seqs_df[,2])
  } else {
    seq_all <- data.frame(sequence=seqs_df[,1], hapid = seqs_df[,2])
  }
  if (parallel) {
    workers <- register_bounded_backend()
    message("Using ", workers, " bounded alignment worker(s)")
    align_df <- foreach(seq_1=seq_along(seq_all$sequence), .combine = "rbind") %dopar% {
      aln <- pw_pairwiseAlignment(ref_str, seq_all$sequence[seq_1], substitutionMatrix = sigma, gapOpening = -8, gapExtension = -5, scoreOnly = FALSE)
      alignment_row(aln, seq_all$hapid[seq_1])
    }
  } else {
    for (seq_1 in seq_along(seq_all$sequence)) {
      aln <- pw_pairwiseAlignment(ref_str, seq_all$sequence[seq_1], substitutionMatrix = sigma, gapOpening = -8, gapExtension = -5, scoreOnly = FALSE)
      df <- alignment_row(aln, seq_all$hapid[seq_1])
      align_df <- rbind(align_df,df)
    }
  }
  return(align_df)
}


# postProc Begin
if (!is.null(args$output)) {
  output <- args$output
  } else {
    stop("Output filename not provided!")
  }

if (!is.null(args$seqtab)) {
  seqfile <- args$seqtab
  if (file.exists(seqfile)) {
	seqtab_input <- fread(seqfile, data.table = FALSE, check.names = FALSE)
    if (ncol(seqtab_input) < 2) {
      stop("ASV sequence table must contain a sample column and at least one ASV column")
    }
    sample_ids <- trimws(as.character(seqtab_input[[1]]))
    if (any(sample_ids == "") || anyDuplicated(sample_ids)) {
      stop("ASV sequence table sample IDs must be non-empty and unique")
    }
    seqtab_input[[1]] <- NULL
    seqtab <- as.matrix(seqtab_input)
    suppressWarnings(storage.mode(seqtab) <- "numeric")
    if (any(!is.finite(seqtab)) || any(seqtab < 0) || any(seqtab != floor(seqtab))) {
      stop("ASV sequence table counts must be finite, non-negative integers")
    }
    rownames(seqtab) <- sample_ids
    seqs <- colnames(seqtab)
	nsample=nrow(seqtab)
	hapid <- paste0("ASV",1:length(seqs))
  	# DataFrame for aligning to truth set
  	seqs_df <- data.frame(sequence = seqs, hapid = hapid)
  	# Change colnames of ASV from sequences to ASV ids
	seqtab_haps <- seqtab
	colnames(seqtab_haps) <- hapid

    if (is.null(args$samples) || !file.exists(args$samples)) {
      stop("Sample sheet (--samples) is required to calculate biological-sample ASV totals")
    }
    sample_sheet <- fread(args$samples, data.table = FALSE, colClasses = "character")
    required_sample_columns <- c(
      "sample_id", "biological_sample_id", "sample_type", "participant_id", "collection_date"
    )
    missing_sample_columns <- setdiff(required_sample_columns, colnames(sample_sheet))
    if (length(missing_sample_columns) > 0) {
      stop(paste("Sample sheet is missing required columns:", paste(missing_sample_columns, collapse = ", ")))
    }
    sample_sheet$sample_id <- trimws(as.character(sample_sheet$sample_id))
    if (any(sample_sheet$sample_id == "") || anyDuplicated(sample_sheet$sample_id)) {
      stop("Sample sheet sample_id values must be non-empty and unique")
    }
    sample_match <- match(rownames(seqtab_haps), sample_sheet$sample_id)
    if (any(is.na(sample_match))) {
      stop(paste(
        "ASV sequence table contains sample IDs absent from samples.csv:",
        paste(rownames(seqtab_haps)[is.na(sample_match)], collapse = ", ")
      ))
    }
    extra_sample_ids <- setdiff(sample_sheet$sample_id, rownames(seqtab_haps))
    if (length(extra_sample_ids) > 0) {
      stop(paste(
        "samples.csv contains sample IDs absent from the ASV sequence table:",
        paste(extra_sample_ids, collapse = ", ")
      ))
    }
    matched_samples <- sample_sheet[sample_match, , drop = FALSE]
    biological_mask <- (
      tolower(trimws(matched_samples$sample_type)) == "sample" &
      trimws(matched_samples$participant_id) != "" &
      trimws(matched_samples$collection_date) != ""
    )
    if (!any(biological_mask)) {
      stop("No dated biological samples remain after excluding controls and unassigned rows")
    }

    biological_ids <- trimws(as.character(matched_samples$biological_sample_id[biological_mask]))
    if (any(biological_ids == "")) {
      stop("Dated biological sample rows must have a biological_sample_id")
    }
    biological_identity <- unique(matched_samples[
      biological_mask,
      c("biological_sample_id", "participant_id", "collection_date"),
      drop = FALSE
    ])
    if (anyDuplicated(biological_identity$biological_sample_id)) {
      stop("A biological_sample_id maps to more than one participant or collection date")
    }

	## ASV summary table
	all_total_reads <- apply(seqtab_haps, 2, sum)
	all_total_samples <- apply(seqtab_haps, 2, function(x) sum(x != 0))
    biological_seqtab <- seqtab_haps[biological_mask, , drop = FALSE]
    biological_sample_seqtab <- rowsum(
      biological_seqtab,
      group = biological_ids,
      reorder = FALSE
    )
	total_reads <- apply(biological_sample_seqtab, 2, sum)
	total_samples <- apply(biological_sample_seqtab, 2, function(x) sum(x != 0))
	asvdf <- data.frame(hapid = hapid,
					  haplength = nchar(seqs),
                      total_reads = total_reads,
                      total_samples = total_samples,
                      qc_total_reads_all_rows = all_total_reads,
                      qc_total_samples_all_rows = all_total_samples,
                      biological_rows_used = sum(biological_mask),
                      biological_samples_used = nrow(biological_sample_seqtab),
                      all_rows_observed = nrow(seqtab_haps),
                      strain = "N")
  	asvdf$hapid <- as.character(asvdf$hapid)
  	asvdf$strain <- as.character(asvdf$strain)

  } else {
    stop(paste("ASV sequence table file",seqfile,"not found!"))
  }
} else {
  stop("Sequence table file (--seqtab) is required")
}

if (!args$no_reference) {
	if (!is.null(args$reference)) {
  		path_to_refseq <- args$reference
  		if (!is.null(args$strain)) {
    		strains <- args$strain
  		} else {
    		stop("Name of target strain (--strain) is required")
  		}
  		if (!is.null(args$reference2)) {
    		path_to_refseq <- c(args$reference,args$reference2)
    		if (!is.null(args$strain2)) {
      			strains <- c(args$strain,args$strain2)
    		} else {
      			stop("Name of second target strain (--strain2) required if --reference2 is given")
    		}
  		}
	} else {
  		stop("Reference fasta file with target sequences (--reference) is required")
	}
	for (p in seq_along(path_to_refseq)) {
  		refasta <- read.fasta(path_to_refseq[p])
  		refseq <- toupper(sapply(refasta,c2s))
  		amplicons <- as.character(names(refseq))
  		# Alignment with RefSet
		align_df <- seq_align(seqs_df = seqs_df, path_to_ref = path_to_refseq[p], overlap = TRUE, parallel = args$parallel)
		align_df$refid <- as.character(align_df$refid)
		align_df$hapid <- as.character(align_df$hapid)

		## Map truthset onto ASV summary table based on exact and inexact matches to truth set
		alignment_columns <- c(
        "hapid", "refid", "snv_dist", "indel_dist", "aln_score",
        "second_aln_score", "alignment_score_margin", "best_reference_ties"
      )
      missing_alignment_columns <- setdiff(alignment_columns, colnames(align_df))
      if (length(missing_alignment_columns) > 0) {
        stop(paste0(
          "Alignment result is missing required columns: ",
          paste(missing_alignment_columns, collapse = ", "),
          ". Available columns: ",
          paste(colnames(align_df), collapse = ", ")
        ))
      }
		df <- align_df[, alignment_columns, drop = FALSE]
		colnames(df) <- c(
        "hapid",
        paste0("refid_", strains[p]),
        paste0("snv_dist_from_", strains[p]),
        paste0("indel_dist_from_", strains[p]),
        paste0("aln_score_", strains[p]),
        paste0("second_aln_score_", strains[p]),
        paste0("alignment_score_margin_", strains[p]),
        paste0("best_reference_ties_", strains[p])
      )
		asvdf <- merge(asvdf, df, by = "hapid", sort = FALSE)
		aligned_index <- match(asvdf$hapid, align_df$hapid)
		asvdf$strain[(align_df$snv_dist[aligned_index] == 0 & align_df$indel_dist[aligned_index] == 0)] <- as.character(strains[p])
	}

  if (is.null(args$snv_filter) || !file.exists(args$snv_filter)) {
    stop("A readable --snv_filter table with id, SNP, and INDEL columns is required")
  }
  VariantCounts <- fread(args$snv_filter, data.table = FALSE)
  required_filter_columns <- c("id", "SNP", "INDEL")
  missing_filter_columns <- setdiff(required_filter_columns, colnames(VariantCounts))
  if (length(missing_filter_columns) > 0) {
    stop(paste("SNV/INDEL filter table is missing columns:", paste(missing_filter_columns, collapse = ", ")))
  }
  VariantCounts$id <- trimws(as.character(VariantCounts$id))
  VariantCounts$SNP <- suppressWarnings(as.numeric(VariantCounts$SNP))
  VariantCounts$INDEL <- suppressWarnings(as.numeric(VariantCounts$INDEL))
  if (any(VariantCounts$id == "") || anyDuplicated(VariantCounts$id) ||
      any(!is.finite(VariantCounts$SNP)) || any(!is.finite(VariantCounts$INDEL)) ||
      any(VariantCounts$SNP < 0) || any(VariantCounts$INDEL < 0)) {
    stop("SNV/INDEL filter table IDs must be unique and thresholds must be finite non-negative values")
  }

  indel_override <- NULL
  if (!is.null(args$indel_filter) && nzchar(trimws(args$indel_filter))) {
    indel_override <- suppressWarnings(as.numeric(args$indel_filter))
    if (length(indel_override) != 1 || !is.finite(indel_override) || indel_override < 0 || indel_override != floor(indel_override)) {
      stop("--indel_filter must be a non-negative absolute INDEL count; fractional length-ratio thresholds are not supported")
    }
  }

  primary_strain <- strains[1]
  reference_column <- paste0("refid_", primary_strain)
  snv_distance_column <- paste0("snv_dist_from_", primary_strain)
  indel_distance_column <- paste0("indel_dist_from_", primary_strain)
  tie_column <- paste0("best_reference_ties_", primary_strain)
  asvdf[[paste0("best_refid_", primary_strain)]] <- as.character(asvdf[[reference_column]])
  asvdf$snv_threshold <- NA_real_
  asvdf$indel_threshold <- NA_real_
  asvdf$snv_filter <- "FAIL"
  asvdf$indel_filter <- "FAIL"
  asvdf$mapping_status <- "UNMAPPED"

  for (i in seq_len(nrow(asvdf))) {
    refid <- as.character(asvdf[[reference_column]][i])
    threshold_index <- match(refid, VariantCounts$id)
    if (is.na(threshold_index)) {
      asvdf$mapping_status[i] <- "MISSING_THRESHOLD"
      next
    }
    snv_threshold <- VariantCounts$SNP[threshold_index]
    indel_threshold <- if (is.null(indel_override)) VariantCounts$INDEL[threshold_index] else indel_override
    asvdf$snv_threshold[i] <- snv_threshold
    asvdf$indel_threshold[i] <- indel_threshold
    asvdf$snv_filter[i] <- if (asvdf[[snv_distance_column]][i] <= snv_threshold) "PASS" else "FAIL"
    asvdf$indel_filter[i] <- if (asvdf[[indel_distance_column]][i] <= indel_threshold) "PASS" else "FAIL"
    if (asvdf[[tie_column]][i] > 1) {
      asvdf$mapping_status[i] <- "AMBIGUOUS"
    } else if (asvdf$snv_filter[i] == "PASS" && asvdf$indel_filter[i] == "PASS") {
      asvdf$mapping_status[i] <- "MAPPED"
    }
  }
  asvdf[[reference_column]][asvdf$mapping_status != "MAPPED"] <- NA_character_
}

if (!is.null(args$bimera)) {
  if (file.exists(args$bimera)) {
    bimera <- fread(args$bimera)
    seqs_df <- merge(seqs_df,bimera, by = "sequence", all = TRUE, sort = FALSE)
    asvdf <- merge(asvdf,seqs_df[,2:3], by = "hapid", all = TRUE, sort = FALSE)
  } else {
    warning(paste("File",args$bimera,"not found. Skipping bimera flag.."))
  }
} else {
  print("Bimeric ASV file not given. Skipping bimera flag..")
}

write.table(asvdf, file = output, sep = "\t", quote = FALSE, row.names = FALSE)

if (args$fasta) {
  write.fasta(lapply(seqs, s2c), names = hapid, file.out = paste0(dirname(output),"/ASVSeqs.fasta"), nbchar = 600)
}
