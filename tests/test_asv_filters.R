args <- commandArgs(trailingOnly = TRUE)
source(file.path(args[[1]], "workflow", "scripts", "lib", "asv_filters.R"))

expect_error <- function(expr, message) {
  error <- tryCatch({ force(expr); NULL }, error = identity)
  stopifnot(inherits(error, "error"), grepl(message, conditionMessage(error), fixed = TRUE))
}

asvs <- data.frame(
  hapid = paste0("ASV", 1:10),
  haplength = c(179, 180, 200, 222, 223, 90, 100, 111, 112, 200),
  refid_3D7 = c(rep("long", 5), rep("short", 4), "exact_upper"),
  snv_dist_from_3D7 = 0,
  indel_dist_from_3D7 = 0,
  best_reference_ties_3D7 = 1
)
# Deliberately different orders test matching by ID rather than row position.
limits <- data.frame(id = c("short", "exact_upper", "long"), SNP = c(8, 12, 28))
lengths <- c(exact_upper = 180, long = 200, short = 100)
filtered <- apply_reference_filters(asvs, limits, lengths)
stopifnot(identical(filtered$indel_filter,
                    c("FAIL", "PASS", "PASS", "PASS", "FAIL", "PASS", "PASS", "PASS", "FAIL", "PASS")))
stopifnot(identical(filtered$reference_length, c(rep(200, 5), rep(100, 4), 180)))
stopifnot(all(filtered$indel_ratio_min == 0.90), all(filtered$indel_ratio_max == 1 / 0.90))
stopifnot(all(filtered$indel_filter_method == "broad_length_ratio"))
stopifnot(all(is.na(filtered$refid_3D7[c(1, 5, 9)])))
stopifnot(identical(filtered$best_refid_3D7, asvs$refid_3D7))

# An unused INDEL column must not alter either filtering decision.
limits$INDEL <- c("not a cutoff", "0", "999")
stopifnot(identical(filtered, apply_reference_filters(asvs, limits, lengths)))

# An equal-length ASV passes the length check even with compensating gaps.
diagnostic <- asvs[3, ]
diagnostic$indel_dist_from_3D7 <- 40
result <- apply_reference_filters(diagnostic, limits, lengths)
stopifnot(result$indel_filter == "PASS", result$indel_dist_from_3D7 == 40)
diagnostic$snv_dist_from_3D7 <- 29
result <- apply_reference_filters(diagnostic, limits, lengths)
stopifnot(result$snv_filter == "FAIL", result$indel_filter == "PASS", result$mapping_status == "UNMAPPED")

# Preserve ambiguous-target and missing-threshold handling.
diagnostic$snv_dist_from_3D7 <- 0
diagnostic$best_reference_ties_3D7 <- 2
result <- apply_reference_filters(diagnostic, limits, lengths)
stopifnot(result$mapping_status == "AMBIGUOUS", is.na(result$refid_3D7))
result <- apply_reference_filters(asvs, limits[limits$id != "long", ], lengths)
stopifnot(all(result$mapping_status[1:5] == "MISSING_THRESHOLD"))

# Explicit overrides are length ratios, including exact-length-only mode.
result <- apply_reference_filters(asvs, limits, lengths, length_ratio_min = 1)
stopifnot(identical(which(result$indel_filter == "PASS"), c(3L, 7L)))
result <- apply_reference_filters(asvs, limits, lengths, length_ratio_min = "0.95")
stopifnot(result$indel_filter[2] == "FAIL", result$indel_filter[3] == "PASS")

for (value in list(0, -0.1, 1.1, 3, NA_real_, Inf, "per_locus_snv_filters", "", numeric(0), c(0.9, 0.95))) {
  expect_error(apply_reference_filters(asvs, limits, lengths, length_ratio_min = value), "length-ratio cutoff")
}
expect_error(apply_reference_filters(asvs, limits["id"], lengths), "missing columns: SNP")
expect_error(apply_reference_filters(asvs, rbind(limits, limits[1, ]), lengths), "IDs must be unique")
bad_limits <- limits
bad_limits$SNP[1] <- -1
expect_error(apply_reference_filters(asvs, bad_limits, lengths), "thresholds must be finite non-negative")
expect_error(apply_reference_filters(asvs, limits, lengths[-2]), "No reference length available")
expect_error(apply_reference_filters(asvs, limits, c(long = 0)), "Reference lengths must be positive")
bad_asvs <- asvs
bad_asvs$haplength[1] <- 0
expect_error(apply_reference_filters(bad_asvs, limits, lengths), "ASV lengths must be finite and positive")

cat("PASS: ratio boundaries, ID matching, unused INDEL column, diagnostic gaps, SNP filtering, mapping status, overrides, invalid inputs\n")
