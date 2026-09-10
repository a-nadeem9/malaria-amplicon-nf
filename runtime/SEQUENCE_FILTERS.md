# Reference sequence filters

Reference mapping applies two separate filters before homopolymer masking and
CIGAR conversion:

- The `SNP` column of `reference/snv_filters.txt` limits mismatching bases per locus.
- The INDEL filter compares the ASV's sequence length with the selected primary
  reference's length. Its default cutoff is `0.90`:

  `0.90 <= ASV length / reference length <= 1 / 0.90`

Both boundaries are inclusive. For a 200-base reference, lengths 180 through 222
pass; 179 and 223 fail. This is Broad's original reciprocal length-ratio check,
not a sum of inserted and deleted alignment bases. Equal-length sequences can
pass this check even when their alignment contains insertions and deletions.
The ratio is calculated from the sequence supplied to postprocessing, after the
workflow's merging or non-overlap preparation, not from the raw read length.

The original INDEL column remains in the reference table for preservation. It
is not used as an absolute cutoff. Tables containing only `id` and `SNP` are also
accepted. The aligned `indel_dist_from_*` counts remain available as diagnostics.

The mapped ASV output records `reference_length`, `indel_length_ratio`,
`indel_ratio_min`, `indel_ratio_max`, and `indel_filter_method`, as well as the
existing SNP/INDEL PASS/FAIL flags and mapping status. The run's recorded
`analysis_parameters.indel_filter` is the numerical ratio cutoff. The old
`indel_threshold` output field is replaced by the explicitly named ratio bounds.

This implements the length filter in Broad's
[original postprocessor](https://github.com/broadinstitute/malaria-amplicon-pipeline/blob/6183028737d2c69b619c57b76ef62aa4160ac675/postProc_dada2.R#L177-L195).
[SIMPLseq's methods](https://link.springer.com/article/10.1186/s12936-026-05796-1)
report a 10% INDEL exclusion and cite that pipeline; the study's exact executed
command has not been independently recovered. This implementation is therefore
documented as Broad-compatible rather than a verified reproduction of that run.

Existing saved results are not rewritten. Reprocessing changes sequence-QC
decisions and requires new downstream CIGAR and allele outputs. An AMA1-only
comparison does not validate the six SIMPLseq panel targets.
