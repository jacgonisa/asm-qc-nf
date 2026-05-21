nextflow.enable.dsl = 2

process COMPLEASM {
    tag "${meta.id}"
    label 'process_medium'

    container 'quay.io/biocontainers/compleasm:0.2.6--pyh7cba7a3_0'

    publishDir "${params.outdir}/${meta.id}/${meta.type}/compleasm", mode: 'symlink'

    input:
    tuple val(meta), path(assembly_fasta)

    output:
    tuple val(meta), path("compleasm_out/"),           emit: results
    tuple val(meta), path("compleasm_out/summary.txt"), emit: summary
    path "versions.yml",                               emit: versions

    script:
    def lineage = meta.busco_lineage ?: params.busco_lineage
    """
    compleasm run \\
        -a ${assembly_fasta} \\
        -o compleasm_out \\
        -l ${lineage} \\
        -t ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        compleasm: \$(compleasm --version 2>&1 | head -1)
    END_VERSIONS
    """
}
