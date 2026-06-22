nextflow.enable.dsl = 2

/*
 * Flagger — read-alignment based assembly error detection.
 * Workflow: minimap2 align reads → samtools sort/index → flagger
 *
 * Uses minimap2 + samtools containers; flagger via conda or local install.
 */

process FLAGGER_ALIGN {
    tag "${meta.id}"
    label 'process_high'

    container 'quay.io/biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:3161f532a5ea6f1dec9be5667c9efc2afdac6104-0'

    input:
    tuple val(meta), path(assembly_fasta), path(reference_fasta)

    output:
    tuple val(meta), path(assembly_fasta), path("${meta.id}.sorted.bam"), path("${meta.id}.sorted.bam.bai"), emit: bam

    script:
    """
    minimap2 \\
        -ax map-hifi \\
        -t ${task.cpus} \\
        ${assembly_fasta} \\
        ${reference_fasta} \\
    | samtools sort -@ ${task.cpus} -o ${meta.id}.sorted.bam

    samtools index ${meta.id}.sorted.bam
    """
}

process FLAGGER_RUN {
    tag "${meta.id}"
    label 'process_medium'

    conda 'flagger'

    publishDir "${params.outdir}/${meta.id}/${meta.type}/flagger", mode: 'symlink'

    input:
    tuple val(meta), path(assembly_fasta), path(bam), path(bai)

    output:
    tuple val(meta), path("flagger_out/"),            emit: results
    tuple val(meta), path("flagger_out/*.bed"),       emit: bed, optional: true
    path "versions.yml",                              emit: versions

    script:
    """
    mkdir -p flagger_out

    flagger \\
        --bam ${bam} \\
        --fasta ${assembly_fasta} \\
        --output flagger_out/${meta.id} \\
        --threads ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        flagger: \$(flagger --version 2>&1 | head -1 || echo 'flagger')
    END_VERSIONS
    """
}

workflow FLAGGER {
    take:
    ch_input  // [ meta, assembly_fasta, reference_fasta ]

    main:
    FLAGGER_ALIGN(ch_input)
    FLAGGER_RUN(FLAGGER_ALIGN.out.bam)

    emit:
    results  = FLAGGER_RUN.out.results
    bed      = FLAGGER_RUN.out.bed
    versions = FLAGGER_RUN.out.versions
}
