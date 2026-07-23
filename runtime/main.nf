#!/usr/bin/env nextflow

process PREFLIGHT {
  tag 'preflight'
  publishDir "${params.outdir}/qc", mode: 'copy', overwrite: true

  input:
  path samples
  path amplicons

  output:
  path 'preflight.tsv', emit: report
  path 'amplicon_geometry.tsv', emit: geometry
  path 'barcode_readiness.tsv', emit: barcode

  script:
  def samples_root = file(params.samples).isAbsolute() ? file(params.samples).parent : "${projectDir}/${file(params.samples).parent}"
  """
  python -m simplseq.progress emit --file "${params.outdir}/progress.jsonl" --stage prepare_inputs --status started --message "Checking sample sheet"
  python "\${SIMPLSEQ_NF_HELPERS}" preflight \
    --samples "${samples}" \
    --amplicons "${amplicons}" \
    --report preflight.tsv \
    --geometry amplicon_geometry.tsv \
    --barcode barcode_readiness.tsv \
    --inline-barcodes-enabled ${params.inline_barcodes_enabled} \
    --sentinel-locus "${params.sentinel_locus}" \
    --samples-root "${samples_root}"
  python -m simplseq.progress emit --file "${params.outdir}/progress.jsonl" --stage prepare_inputs --status complete
  """
}

process WRITE_META {
  tag 'meta'
  publishDir "${params.outdir}", mode: 'copy', overwrite: true

  input:
  path samples

  output:
  path 'meta.tsv', emit: meta

  script:
  def samples_root = file(params.samples).isAbsolute() ? file(params.samples).parent : "${projectDir}/${file(params.samples).parent}"
  """
  python "\${SIMPLSEQ_NF_HELPERS}" write-meta --samples "${samples}" --out meta.tsv --samples-root "${samples_root}"
  """
}

process WRITE_PIPELINE_JSON {
  tag 'pipeline_json'
  publishDir "${params.outdir}/config", mode: 'copy', overwrite: true

  input:
  path meta

  output:
  path 'pipeline_inputs.json', emit: json

  script:
  def dada2_randomize_arg = params.dada2_randomize ? "--dada2-randomize ${params.dada2_randomize}" : ''
  def dada2_multithread_arg = params.dada2_multithread ? "--dada2-multithread ${params.dada2_multithread}" : ''
  def dada2_seed_arg = params.dada2_seed ? "--dada2-seed ${params.dada2_seed}" : ''
  """
  python "\${SIMPLSEQ_NF_HELPERS}" write-pipeline-json \
    --meta "${meta}" \
    --out pipeline_inputs.json \
    --pipeline-class "${params.pipeline_class}" \
    --max-ee ${params.max_ee} \
    --trim-right ${params.trim_right} \
    --min-len ${params.min_len} \
    --trunc-q ${params.trunc_q} \
    --max-consist ${params.max_consist} \
    --omega-a ${params.omega_a} \
    --just-concatenate ${params.just_concatenate} \
    --save-rdata ${params.save_rdata} \
    ${dada2_randomize_arg} \
    ${dada2_multithread_arg} \
    ${dada2_seed_arg} \
    --primers-fwd "${projectDir}/${params.primers_fwd}" \
    --primers-rev "${projectDir}/${params.primers_rev}" \
    --overlap-primers-fwd "${projectDir}/${params.overlap_primers_fwd}" \
    --overlap-primers-rev "${projectDir}/${params.overlap_primers_rev}"
  """
}

process DADA2_PIPELINE {
  tag 'dada2'
  cpus params.dada2_cpus
  memory params.dada2_memory
  time '6h'
  publishDir "${params.outdir}/run_dada2/intermediate/dada2_op", mode: 'copy', overwrite: true, pattern: 'op_ASVBimeras.txt', saveAs: { 'ASVBimeras.txt' }
  publishDir "${params.outdir}/run_dada2/intermediate/dada2_nop", mode: 'copy', overwrite: true, pattern: 'nop_ASVBimeras.txt', saveAs: { 'ASVBimeras.txt' }
  publishDir "${params.outdir}/run_dada2/intermediate/dada2_nop", mode: 'copy', overwrite: true, pattern: 'correctedASV.txt'
  publishDir "${params.outdir}/logs", mode: 'copy', overwrite: true, pattern: 'dada2_pipeline.log'

  input:
  path json
  path meta
  path reference
  path overlap_fwd
  path overlap_rev
  path preflight_report
  path input_fastq_md5s

  output:
  path 'seqtab_iseq.tsv', emit: seqtab
  path 'op_ASVBimeras.txt', emit: op_bimera
  path 'nop_ASVBimeras.txt', emit: nop_bimera
  path 'correctedASV.txt', emit: corrected_asv
  path 'preprocess_fq', emit: preprocessed_fastqs
  path 'dada2_pipeline.log', emit: log

  script:
  def bimera_flag = params.cigar_exclude_bimeras ? '--exclude-bimeras' : ''
  """
  python -m simplseq.progress emit --file "${params.outdir}/progress.jsonl" --stage dada2 --status started --message "Running DADA2"
  {
    echo "[SIMPLseq/Nextflow] \$(date) starting DADA2/AmpliconPipeline"
    test -s "${input_fastq_md5s}"
    SIMPLSEQ_THREADS=${task.cpus} python "\${SIMPLSEQ_AMP_PIPELINE}" \
      --json "${json}" \
      --workers ${task.cpus} \
      --iseq \
      --reference "${reference}" \
      --overlap_pr1 "${overlap_fwd}" \
      --overlap_pr2 "${overlap_rev}"
    echo "[SIMPLseq/Nextflow] \$(date) finished DADA2/AmpliconPipeline"
  } 2>&1 | tee dada2_pipeline.log

  cp run_dada2/seqtab_iseq.tsv seqtab_iseq.tsv
  cp run_dada2/dada2_op/ASVBimeras.txt op_ASVBimeras.txt
  cp run_dada2/dada2_nop/ASVBimeras.txt nop_ASVBimeras.txt
  cp run_dada2/dada2_nop/correctedASV.txt correctedASV.txt
  python -m simplseq.progress emit --file "${params.outdir}/progress.jsonl" --stage dada2 --status complete
  """
}

process KELT_CONTAMINATION {
  tag 'kelt_qc'
  cpus params.kelt_cpus
  memory '4 GB'
  time '2h'
  publishDir "${params.outdir}/qc/kelt", mode: 'copy', overwrite: true

  input:
  path preprocessed_fastqs
  path samples

  output:
  path 'kelt_contamination_calls.tsv', emit: calls
  path 'kelt_barcode_counts.tsv', emit: counts
  path 'kelt_contamination_summary.json', emit: summary
  path 'kelt_contamination_report.html', emit: report
  path 'kelt_contamination.log', emit: log

  script:
  """
  {
    echo "[SIMPLseq/Nextflow] \$(date) starting KELT contamination QC"
    python "\${SIMPLSEQ_KELT_DETECTION}" \
      --fastq-dir "${preprocessed_fastqs}" \
      --samples "${samples}" \
      --outdir . \
      --mismatches ${params.kelt_barcode_mismatches} \
      --trace-max-reads ${params.kelt_trace_max_reads} \
      --workers ${task.cpus}
    echo "[SIMPLseq/Nextflow] \$(date) finished KELT contamination QC"
  } 2>&1 | tee kelt_contamination.log
  """
}

process PREPARE_STAGE2 {
  tag 'prepare_stage2'
  publishDir "${params.outdir}/run_dada2/intermediate", mode: 'copy', overwrite: true, pattern: 'ASVBimeras_strict.txt'

  input:
  path seqtab
  path op_bimera
  path nop_bimera
  path corrected_asv

  output:
  path 'seqtab_iseq_strict_fixed.tsv', emit: strict_seqtab
  path 'ASVBimeras_strict.txt', emit: strict_bimera

  script:
  """
  python -m simplseq.progress emit --file "${params.outdir}/progress.jsonl" --stage prepare_stage2 --status started --message "Cleaning ASV table"
  python "\${SIMPLSEQ_NF_HELPERS}" prepare-stage2 \
    --seqtab "${seqtab}" \
    --op-bimera "${op_bimera}" \
    --nop-bimera "${nop_bimera}" \
    --corrected-asv "${corrected_asv}" \
    --strict-seqtab seqtab_iseq_strict_fixed.tsv \
    --strict-bimera ASVBimeras_strict.txt \
    --strict-min-asv-length ${params.strict_min_asv_length}
  python -m simplseq.progress emit --file "${params.outdir}/progress.jsonl" --stage prepare_stage2 --status complete
  """
}

process POSTPROCESS_DADA2 {
  tag 'postprocess'
  cpus params.postprocess_cpus
  memory params.postprocess_memory
  time '4h'
  publishDir "${params.outdir}/run_dada2/intermediate", mode: 'copy', overwrite: true, pattern: 'ASVSeqs.fasta'
  publishDir "${params.outdir}/logs", mode: 'copy', overwrite: true, pattern: 'postprocess_dada2.log'

  input:
  path seqtab
  path bimera
  path reference
  path snv_filters
  path samples

  output:
  path 'ASV_mapped_table.tsv', emit: mapped
  path 'ASVSeqs.fasta', emit: fasta
  path 'postprocess_dada2.log', emit: log

  script:
  """
  python -m simplseq.progress emit --file "${params.outdir}/progress.jsonl" --stage asv_mapping --status started --message "Mapping ASVs"
  {
    echo "[SIMPLseq/Nextflow] \$(date) starting postprocess_dada2"
    SIMPLSEQ_THREADS=${task.cpus} Rscript "\${SIMPLSEQ_POSTPROC_DADA2}" \
      -s "${seqtab}" \
      --samples "${samples}" \
      -ref "${reference}" \
      -b "${bimera}" \
      -o ASV_mapped_table.tsv \
      --fasta \
      --parallel \
      --snv_filter "${snv_filters}"
    echo "[SIMPLseq/Nextflow] \$(date) finished postprocess_dada2"
  } 2>&1 | tee postprocess_dada2.log
  python -m simplseq.progress emit --file "${params.outdir}/progress.jsonl" --stage asv_mapping --status complete
  """
}

process PREPARE_STAGE3 {
  tag 'prepare_stage3'
  publishDir "${params.outdir}/run_dada2", mode: 'copy', overwrite: true

  input:
  path mapped
  path seqtab

  output:
  path 'ASV_mapped_table.tsv', emit: filtered_mapped
  path 'seqtab_iseq.tsv', emit: fixed_seqtab

  script:
  """
  python -m simplseq.progress emit --file "${params.outdir}/progress.jsonl" --stage prepare_stage3 --status started --message "Preparing CIGAR inputs"
  python "\${SIMPLSEQ_NF_HELPERS}" prepare-stage3 \
    --mapped "${mapped}" \
    --seqtab "${seqtab}" \
    --filtered-mapped ASV_mapped_table.tsv \
    --fixed-seqtab seqtab_iseq.tsv
  python -m simplseq.progress emit --file "${params.outdir}/progress.jsonl" --stage prepare_stage3 --status complete
  """
}

process CHECK_CIGAR_INPUTS {
  tag 'check_cigar_inputs'
  publishDir "${params.outdir}/qc", mode: 'copy', overwrite: true

  input:
  path table

  output:
  path 'cigar_input_summary.tsv', emit: summary

  script:
  def bimera_flag = params.cigar_exclude_bimeras ? '--exclude-bimeras' : ''
  """
  python -m simplseq.progress emit --file "${params.outdir}/progress.jsonl" --stage cigar_check --status started --message "Checking CIGAR inputs"
  python "\${SIMPLSEQ_NF_HELPERS}" check-cigar-inputs \
    --table "${table}" \
    --summary cigar_input_summary.tsv \
    --min-reads ${params.cigar_min_total_reads} \
    --min-samples ${params.cigar_min_samples} \
    ${bimera_flag}
  python -m simplseq.progress emit --file "${params.outdir}/progress.jsonl" --stage cigar_check --status complete
  """
}

process CIGAR_CONVERSION {
  tag 'cigar'
  cpus params.cigar_cpus
  memory params.cigar_memory
  time '1h'
  publishDir "${params.outdir}/run_dada2", mode: 'copy', overwrite: true, pattern: 'seqtab_cigar.tsv'
  publishDir "${params.outdir}/run_dada2", mode: 'copy', overwrite: true, pattern: 'asv_to_cigar.tsv'
  publishDir "${params.outdir}/logs", mode: 'copy', overwrite: true, pattern: 'cigar_conversion.log'

  input:
  path fasta
  path table
  path seqtab
  path reference
  path summary

  output:
  path 'seqtab_cigar.tsv', emit: seqtab_cigar
  path 'asv_to_cigar.tsv', emit: asv_to_cigar
  path 'cigar_conversion.log', emit: log

  script:
  def bimera_flag = params.cigar_exclude_bimeras ? '-b' : ''
  """
  python -m simplseq.progress emit --file "${params.outdir}/progress.jsonl" --stage asv_to_cigar --status started --message "Converting ASVs to CIGAR"
  {
    echo "[SIMPLseq/Nextflow] \$(date) starting CIGAR conversion"
    python "\${SIMPLSEQ_ASV_TO_CIGAR}" \
      "${fasta}" \
      "${table}" \
      "${seqtab}" \
      seqtab_cigar.tsv \
      --asv_to_cigar asv_to_cigar.tsv \
      --amp_db "${reference}" \
      --alignments alignments \
      -p ${params.cigar_homopolymer_mask_length} \
      -r ${params.cigar_min_total_reads} \
      -n ${params.cigar_min_samples} \
      ${bimera_flag}
    echo "[SIMPLseq/Nextflow] \$(date) finished CIGAR conversion"
  } 2>&1 | tee cigar_conversion.log
  python -m simplseq.progress emit --file "${params.outdir}/progress.jsonl" --stage asv_to_cigar --status complete
  """
}

process CDHIT_CLUSTERING {
  tag 'cdhit'
  cpus params.cdhit_cpus
  memory '2 GB'
  time '1h'
  publishDir "${params.outdir}/cdhit", mode: 'copy', overwrite: true

  input:
  path fasta
  path mapped
  path seqtab
  path asv_to_cigar
  path samples

  output:
  path 'filtered_ASVSeqs.fasta', emit: filtered_fasta
  path 'cdhit_representatives.fasta', emit: representatives
  path 'cdhit_representatives.fasta.clstr', emit: raw_clusters
  path 'cdhit_cluster_membership.tsv', emit: membership
  path 'cdhit_cluster_counts.tsv', emit: counts
  path 'cdhit_summary.json', emit: summary
  path 'cdhit.log', emit: log

  script:
  def bimera_flag = params.cigar_exclude_bimeras ? '--exclude-bimeras' : ''
  """
  python -m simplseq.progress emit --file "${params.outdir}/progress.jsonl" --stage cdhit --status started --message "Clustering filter-pass ASVs at ${params.cdhit_identity} identity"
  python "\${SIMPLSEQ_CDHIT_CLUSTERING}" \
    --fasta "${fasta}" \
    --mapped "${mapped}" \
    --seqtab "${seqtab}" \
    --samples "${samples}" \
    --asv-to-cigar "${asv_to_cigar}" \
    --outdir . \
    --identity ${params.cdhit_identity} \
    --word-size ${params.cdhit_word_size} \
    --threads ${task.cpus} \
    --min-total-reads ${params.cigar_min_total_reads} \
    --min-samples ${params.cigar_min_samples} \
    ${bimera_flag}
  python -m simplseq.progress emit --file "${params.outdir}/progress.jsonl" --stage cdhit --status complete
  """
}

process RUN_REPORT {
  tag 'run_report'
  publishDir "${params.outdir}/reports", mode: 'copy', overwrite: true

  input:
  path preflight
  path geometry
  path cigar_summary
  path mapped
  path asv_to_cigar
  path cigar
  path cdhit_summary

  output:
  path 'run_summary.html', emit: html
  path 'asv_filtering_summary.txt', emit: asv_filtering_summary

  script:
  def bimera_flag = params.cigar_exclude_bimeras ? '--exclude-bimeras' : ''
  """
  python -m simplseq.progress emit --file "${params.outdir}/progress.jsonl" --stage report --status started --message "Writing report"
  python "\${SIMPLSEQ_NF_HELPERS}" make-report \
    --project-name "${params.project_name}" \
    --preflight "${preflight}" \
    --geometry "${geometry}" \
    --cigar-summary "${cigar_summary}" \
    --mapped "${mapped}" \
    --asv-to-cigar "${asv_to_cigar}" \
    --cigar "${cigar}" \
    --cdhit-summary "${cdhit_summary}" \
    --min-total-reads ${params.cigar_min_total_reads} \
    --min-samples ${params.cigar_min_samples} \
    ${bimera_flag} \
    --samples "${params.samples}" \
    --out run_summary.html
  python -m simplseq.progress emit --file "${params.outdir}/progress.jsonl" --stage report --status complete
  """
}

process DINEMITES_ANALYSIS {
  tag 'dinemites'
  memory '4 GB'
  time '2h'
  publishDir "${params.outdir}/dinemites", mode: 'copy', overwrite: true

  input:
  path seqtab_cigar
  path samples

  output:
  path 'dinemites_allele_probabilities.tsv', emit: allele_probs, optional: true
  path 'dinemites_allele_key.tsv', emit: allele_key, optional: true
  path 'dinemites_molfoi.tsv', emit: molfoi, optional: true
  path 'dinemites_new_infections.tsv', emit: new_infections, optional: true
  path 'dinemites_input.tsv', emit: input, optional: true
  path 'dinemites_preparation.json', emit: preparation, optional: true
  path 'dinemites_summary.json', emit: summary, optional: true
  path 'dinemites_plots/*.png', emit: plots, optional: true
  path 'dinemites_subjects', emit: subject_outputs, optional: true
  path 'dinemites.log', emit: log

  script:
  def samples_path = file(params.samples).isAbsolute() ? params.samples : "${projectDir}/${params.samples}"
  """
  python -m simplseq.progress emit --file "${params.outdir}/progress.jsonl" --stage dinemites --status started --message "Running DINEMITES analysis"
  {
    echo "[SIMPLseq/Nextflow] \$(date) starting DINEMITES data bridge"
    Rscript "\${SIMPLSEQ_DINEMITES_BRIDGE}" \
      --cigar "${seqtab_cigar}" \
      --samples "${samples_path}" \
      --out dinemites_input.tsv \
      --preparation_summary dinemites_preparation.json \
      --min_abundance_pct ${params.dinemites_min_abundance_pct} \
      --abundance_denominator ${params.dinemites_abundance_denominator}

    echo "[SIMPLseq/Nextflow] \$(date) starting DINEMITES model (${params.dinemites_model})"
    Rscript "\${SIMPLSEQ_DINEMITES_RUN}" \
      --input dinemites_input.tsv \
      --model ${params.dinemites_model} \
      --outdir . \
      --n_lags ${params.dinemites_n_lags} \
      --t_lag ${params.dinemites_t_lag} \
      --seed ${params.dinemites_seed} \
      --refresh ${params.dinemites_refresh} \
      --bayesian_lag_days ${params.dinemites_bayesian_lag_days} \
      --bayesian_chains ${params.dinemites_bayesian_chains} \
      --bayesian_parallel_chains ${params.dinemites_bayesian_parallel_chains} \
      --bayesian_iter_warmup ${params.dinemites_bayesian_iter_warmup} \
      --bayesian_iter_sampling ${params.dinemites_bayesian_iter_sampling} \
      --bayesian_adapt_delta ${params.dinemites_bayesian_adapt_delta} \
      --bayesian_drop_out ${params.dinemites_bayesian_drop_out} \
      --plot_width ${params.dinemites_plot_width} \
      --plot_height ${params.dinemites_plot_height}

    echo "[SIMPLseq/Nextflow] \$(date) finished DINEMITES"
  } 2>&1 | tee dinemites.log
  python -m simplseq.progress emit --file "${params.outdir}/progress.jsonl" --stage dinemites --status complete
  """
}

process DCIFER_ANALYSIS {
  tag 'dcifer'
  memory '4 GB'
  time '2h'
  publishDir "${params.outdir}/dcifer", mode: 'copy', overwrite: true

  input:
  path seqtab_cigar
  path samples

  output:
  path 'dcifer_input_long.tsv', emit: input, optional: true
  path 'dcifer_filter_summary.tsv', emit: filter_summary, optional: true
  path 'dcifer_replicate_summary.tsv', emit: replicate_summary, optional: true
  path 'dcifer_coi.tsv', emit: coi, optional: true
  path 'dcifer_allele_frequencies.tsv', emit: allele_frequencies, optional: true
  path 'dcifer_pairwise_relatedness.tsv', emit: pairwise_relatedness, optional: true
  path 'dcifer_relatedness_matrix.tsv', emit: relatedness_matrix, optional: true
  path 'dcifer_pvalue_matrix.tsv', emit: pvalue_matrix, optional: true
  path 'dcifer_summary.json', emit: summary, optional: true
  path 'dcifer_plots/*.png', emit: plots, optional: true
  path 'dcifer.log', emit: log

  script:
  def samples_path = file(params.samples).isAbsolute() ? params.samples : "${projectDir}/${params.samples}"
  """
  python -m simplseq.progress emit --file "${params.outdir}/progress.jsonl" --stage dcifer --status started --message "Running dcifer relatedness analysis"
  {
    echo "[SIMPLseq/Nextflow] \$(date) starting dcifer data bridge"
    Rscript "\${SIMPLSEQ_DCIFER_BRIDGE}" \
      --cigar "${seqtab_cigar}" \
      --samples "${samples_path}" \
      --out dcifer_input_long.tsv \
      --filter_summary dcifer_filter_summary.tsv \
      --replicate_summary dcifer_replicate_summary.tsv \
      --min_abundance_pct ${params.dcifer_min_abundance_pct} \
      --abundance_denominator ${params.dcifer_abundance_denominator}

    echo "[SIMPLseq/Nextflow] \$(date) starting dcifer model"
    Rscript "\${SIMPLSEQ_DCIFER_RUN}" \
      --input dcifer_input_long.tsv \
      --outdir . \
      --coi_lrank ${params.dcifer_coi_lrank} \
      --ibd_grid_nr ${params.dcifer_ibd_grid_nr} \
      --alpha ${params.dcifer_alpha} \
      --afreq_mode ${params.dcifer_afreq_mode}

    echo "[SIMPLseq/Nextflow] \$(date) finished dcifer"
  } 2>&1 | tee dcifer.log
  python -m simplseq.progress emit --file "${params.outdir}/progress.jsonl" --stage dcifer --status complete
  """
}

workflow {
  if (!params.input_fastq_md5s) {
    error "--input_fastq_md5s is required. Start the workflow with simplseq run so FASTQ content fingerprints are recorded."
  }
  samples = file(params.samples)
  input_fastq_md5s = file(params.input_fastq_md5s)
  reference = file(params.amplicons_noprimers)
  snv_filters = file(params.snv_filters)
  overlap_fwd = file(params.overlap_primers_fwd)
  overlap_rev = file(params.overlap_primers_rev)

  preflight = PREFLIGHT(samples, reference)
  meta = WRITE_META(samples)
  pipeline_json = WRITE_PIPELINE_JSON(meta.meta)

  dada2 = DADA2_PIPELINE(
    pipeline_json.json,
    meta.meta,
    reference,
    overlap_fwd,
    overlap_rev,
    preflight.report,
    input_fastq_md5s
  )
  if (params.inline_barcodes_enabled) {
    KELT_CONTAMINATION(dada2.preprocessed_fastqs, samples)
  }
  stage2 = PREPARE_STAGE2(dada2.seqtab, dada2.op_bimera, dada2.nop_bimera, dada2.corrected_asv)
  post = POSTPROCESS_DADA2(stage2.strict_seqtab, stage2.strict_bimera, reference, snv_filters, samples)
  stage3 = PREPARE_STAGE3(post.mapped, stage2.strict_seqtab)
  cigar_check = CHECK_CIGAR_INPUTS(stage3.filtered_mapped)
  cigar = CIGAR_CONVERSION(post.fasta, stage3.filtered_mapped, stage3.fixed_seqtab, reference, cigar_check.summary)
  cdhit = CDHIT_CLUSTERING(post.fasta, stage3.filtered_mapped, stage3.fixed_seqtab, cigar.asv_to_cigar, samples)
  report = RUN_REPORT(preflight.report, preflight.geometry, cigar_check.summary, post.mapped, cigar.asv_to_cigar, cigar.seqtab_cigar, cdhit.summary)
  if (params.dinemites_enabled) {
    DINEMITES_ANALYSIS(cigar.seqtab_cigar, samples)
  }
  if (params.dcifer_enabled) {
    DCIFER_ANALYSIS(cigar.seqtab_cigar, samples)
  }
}
