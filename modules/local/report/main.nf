nextflow.enable.dsl = 2

/*
 * QC_REPORT — Aggregate all tool outputs into a per-assembly PDF report.
 *
 * Pages:
 *   1. Summary stats table (FASTA N50, QUAST, Merqury, compleasm, LAI,
 *      kplexity asymptote, Flagger percentages)
 *   2. kplexity curve with double-sigmoid fit (if kplex CSV present)
 *   3+. Merqury spectrum plots (spectra-asm + spectra-cn, if dir present)
 *   N+. NucFreq plots (if dir present)
 *   N+. Nchart plots (if dir present)
 *
 * All tool inputs are optional — the report degrades gracefully when a tool
 * was disabled (NO_FILE sentinel is passed and silently skipped).
 *
 * Requires: numpy, scipy, matplotlib in the execution environment.
 */

process QC_REPORT {
    tag "${meta.id} [${meta.type}]"
    label 'process_low'

    container 'quay.io/biocontainers/mulled-v2-ad9dd5f398966bf899ae05f8e7c54d0fb10cdfa3:077a1b86b8f7a99cc04dc7f9c1fd4ee4e12e8ed3-0'

    publishDir "${params.outdir}/${meta.id}/${meta.type}/report", mode: 'copy'

    input:
    tuple val(meta),
          path(assembly_fasta),
          path(quast_tsv),
          path(merqury_qv),
          path(merqury_completeness),
          path(merqury_dir),
          path(compleasm_summary),
          path(lai_output),
          path(kplex_csv),
          path(flagger_bed),
          path(nucfreq_dir),
          path(nchart_dir)

    output:
    tuple val(meta), path("${meta.id}_${meta.type}_qc_report.pdf"), emit: pdf
    path "versions.yml",                                            emit: versions

    script:
    def out_pdf  = "${meta.id}_${meta.type}_qc_report.pdf"
    def opt_quast   = quast_tsv.name            != 'NO_FILE' ? "--quast_tsv ${quast_tsv}"                         : ''
    def opt_mqv     = merqury_qv.name           != 'NO_FILE' ? "--merqury_qv ${merqury_qv}"                       : ''
    def opt_mqc     = merqury_completeness.name != 'NO_FILE' ? "--merqury_completeness ${merqury_completeness}"   : ''
    def opt_mqdir   = merqury_dir.name          != 'NO_FILE' ? "--merqury_dir ${merqury_dir}"                     : ''
    def opt_compl   = compleasm_summary.name    != 'NO_FILE' ? "--compleasm_summary ${compleasm_summary}"          : ''
    def opt_lai     = lai_output.name           != 'NO_FILE' ? "--lai_output ${lai_output}"                       : ''
    def opt_kplex   = kplex_csv.name            != 'NO_FILE' ? "--kplex_csv ${kplex_csv}"                         : ''
    def opt_flag    = flagger_bed.name          != 'NO_FILE' ? "--flagger_bed ${flagger_bed}"                     : ''
    def opt_nfdir   = nucfreq_dir.name          != 'NO_FILE' ? "--nucfreq_dir ${nucfreq_dir}"                     : ''
    def opt_ncdir   = nchart_dir.name           != 'NO_FILE' ? "--nchart_dir ${nchart_dir}"                       : ''
    """
    render_qc_report.py \\
        --sample         ${meta.id} \\
        --type           ${meta.type} \\
        --assembly_fasta ${assembly_fasta} \\
        --outfile        ${out_pdf} \\
        ${opt_quast} \\
        ${opt_mqv} \\
        ${opt_mqc} \\
        ${opt_mqdir} \\
        ${opt_compl} \\
        ${opt_lai} \\
        ${opt_kplex} \\
        ${opt_flag} \\
        ${opt_nfdir} \\
        ${opt_ncdir}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        render_qc_report: 1.1.0
        python: \$(python3 --version 2>&1 | awk '{print \$2}')
        scipy: \$(python3 -c "import scipy; print(scipy.__version__)")
    END_VERSIONS
    """
}
