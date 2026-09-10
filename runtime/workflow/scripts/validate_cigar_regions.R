# Shared input contract for the standalone downstream R bridges.
validate_cigar_regions <- function(columns) {
  if (any(!grepl("^[^,]+,.+$", columns))) {
    stop("Allele columns must be encoded as LOCUS,ALLELE.")
  }
  loci <- sub(",.*", "", columns)
  alleles <- sub("^[^,]*,", "", columns)
  is_partial <- startsWith(alleles, "partial@")
  if (any(is_partial & !grepl("^partial@[0-9a-f]{32}:.+$", alleles))) {
    stop("Malformed partial-region allele.")
  }
  definitions <- rep("full", length(alleles))
  definitions[is_partial] <- sub("^partial@([0-9a-f]{32}):.*$", "\\1", alleles[is_partial])
  by_locus <- split(definitions, loci)
  mixed <- names(by_locus)[vapply(by_locus, function(x) length(unique(x)) > 1, logical(1))]
  if (length(mixed)) {
    stop("Cannot compare different observation regions at the same locus: ",
         paste(sort(mixed), collapse = ", "),
         ". Reprocess the samples together to define one shared comparison region.")
  }
  invisible(NULL)
}
