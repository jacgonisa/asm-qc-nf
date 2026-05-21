nextflow.enable.dsl = 2

process QUAST {
    tag "${meta.id}"
    label 'process_medium'

    container 'quay.io/biocontainers/quast:5.2.0--py39pl5321h2add14b_1'

    publishDir "${params.outdir}/${meta.id}/quast", mode: 'symlink'

    input:
    tuple val(meta), path(assembly_fasta), path(reference_fasta)

    output:
    tuple val(meta), path("quast_out/"),                   emit: results
    tuple val(meta), path("quast_out/report.tsv"),         emit: report_tsv
    tuple val(meta), path("quast_out/report.html"),        emit: report_html, optional: true
    path "versions.yml",                                   emit: versions

    script:
    """
    quast.py \\
        ${assembly_fasta} \\
        -r ${reference_fasta} \\
        -o quast_out \\
        -t ${task.cpus} \\
        --large

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quast: \$(quast.py --version 2>&1 | grep -oP 'QUAST v[\\d.]+')
    END_VERSIONS
    """
}
