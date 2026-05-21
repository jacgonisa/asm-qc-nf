nextflow.enable.dsl = 2

process MERQURY {
    tag "${meta.id}"
    label 'process_high'

    container 'quay.io/biocontainers/merqury:1.3--hdfd78af_2'

    publishDir "${params.outdir}/${meta.id}/${meta.type}/merqury", mode: 'symlink'

    input:
    tuple val(meta), path(meryl_db), path(assembly_fasta)

    output:
    tuple val(meta), path("merqury_out/"),                emit: results
    tuple val(meta), path("merqury_out/*.qv"),            emit: qv
    tuple val(meta), path("merqury_out/*.completeness.stats"), emit: completeness
    path "versions.yml",                                  emit: versions

    script:
    """
    mkdir merqury_out
    cd merqury_out
    merqury.sh \\
        ../${meryl_db} \\
        ../${assembly_fasta} \\
        ${meta.id}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        merqury: \$(merqury.sh 2>&1 | grep -i version | head -1 || echo 'merqury 1.3')
    END_VERSIONS
    """
}
