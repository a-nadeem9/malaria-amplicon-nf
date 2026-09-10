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

## Unknown bases in non-overlapping reads

Non-overlapping read pairs can contain an unread middle, represented by `N`
after reference-guided adjustment. Those positions are unknown, not SNP
observations. For a partial ASV, the postprocessor counts substitutions only in
aligned columns where both the reference and query are A/C/G/T. Fully observed
ASVs retain the existing distance calculation. The output adds
`unknown_bases_*` and `comparable_bases_*` for each reference strain. A partial
ASV matching the reference in its observed positions is not assigned an exact
reference-strain identity.

This follows the intent of the newer Broad
[postprocessor](https://github.com/broadinstitute/malaria_ampseq/blob/ee3fc96368172d02112e2d0b1ff28f46cbe6bb58/Code/postProc_dada2.R#L43),
which subtracts Ns before SNP filtering, and its default downstream
[N masking](https://github.com/broadinstitute/malaria_ampseq/blob/ee3fc96368172d02112e2d0b1ff28f46cbe6bb58/Code/MHap/amplseq_functions.R#L33).
The direct observed-column count avoids negative SNP distances reproduced with
the upstream subtraction formula in partial deletion examples.

### Comparable partial alleles

Deleting each ASV's N tokens independently is insufficient: a variant visible
in a longer observation can be unknown in a shorter observation of the same
haplotype. Treating the two remaining strings as complete alleles can create
false differences or reference matches.

For any locus with unknown positions, CIGAR conversion therefore defines one
shared set of callable reference positions from the same MUSCLE alignment used
for variant calls. It projects **every ASV at that locus** onto that set before
summing allele counts and applying downstream abundance or recurrence rules.
Existing coordinate and homopolymer masks also contribute to this set. Gap
runs touching unknown sequence are uncertain; deletions crossing a masked
boundary are excluded as whole events, and an insertion requires both
reference flanks to be callable. Observed differences within the remaining
positions are retained.

Partial labels have the form `partial@<comparison_id>:<CIGAR>`. The identifier
contains a hash of the method version, reference ID, reference sequence hash,
and callable coordinates. Thus `partial@...:.` means reference agreement only
over that declared comparison region. It does not establish complete-locus
identity. A locus with no unknown positions keeps its existing CIGAR labels.

The primary analysis preparation and both standalone DINEMITES/dcifer R
bridges reject different comparison IDs, or full and partial definitions,
mixed at the same locus. Reprocess the underlying observations together to
obtain a common region; do not strip the prefix to combine them. A full AMA1
locus and a partial CSP locus may coexist in one table.

### Coverage audit and limits

Nextflow publishes `run_dada2/cigar_coverage.json`. It records the reference
identity, retained coordinate intervals and length, uncertainty per ASV,
boundary-event mask drivers, and original CIGAR strings. The original ASV
sequences remain in `run_dada2/intermediate/ASVSeqs.fasta`. A direct invocation
of `ASV_to_CIGAR.py` writes `<output>.coverage.json` by default, or the path
given by `--coverage-report`.

The shared mask is derived from ASVs passing the CIGAR input filters, before
biological abundance and recurrence filtering. Under production defaults,
control-only ASVs do not meet the required biological read/sample totals.
A low-support but eligible ASV can still shorten the common region for all
samples. Inspect retained coverage and its ASV drivers when assessing
discrimination; nonempty coverage alone is not evidence of adequate resolution.
If no shared callable reference positions remain, conversion stops with a
clear error rather than assigning a reference allele.

This is a reduced-locus adaptation for partial reads. Differences outside the
common region are deliberately excluded even if some samples observed them;
the method cannot recover the unread sequence. Its allele frequencies and
relatedness resolution are not automatically equivalent to full-locus results.
SIMPLseq reports longer paired sequencing and merged microhaplotypes; its exact
executed N-handling revision remains unverified.

These changes do not add fixed read-end clipping. They act after denoising and
preserve the pooling, `MAX_CONSIST`, and length-ratio INDEL settings. Previously
saved allele tables require reprocessing to use the new definition.
