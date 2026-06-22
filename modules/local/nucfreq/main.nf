nextflow.enable.dsl = 2

/*
 * NucFreq — nucleotide frequency per position, reveals collapses / heterozygosity.
 * Aligns reference reads to assembly, then runs nucfreq on each chromosome window.
 *
 * Input: assembly FASTA + reference FASTA (used as the "read" input for alignment
 *        in assembly-vs-reference mode, or pass HiFi reads directly if available).
 */

process NUCFREQ_ALIGN {
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
        -ax asm5 \\
        -t ${task.cpus} \\
        ${assembly_fasta} \\
        ${reference_fasta} \\
    | samtools sort -@ ${task.cpus} -o ${meta.id}.sorted.bam

    samtools index ${meta.id}.sorted.bam
    """
}

process NUCFREQ_RUN {
    tag "${meta.id}"
    label 'process_medium'

    container 'quay.io/biocontainers/nucfreq:0.1.3--py38h1eb1b71_0'

    publishDir "${params.outdir}/${meta.id}/${meta.type}/nucfreq", mode: 'symlink'

    input:
    tuple val(meta), path(assembly_fasta), path(bam), path(bai)

    output:
    tuple val(meta), path("nucfreq_out/"),           emit: results
    tuple val(meta), path("nucfreq_out/*.png"),      emit: plots, optional: true
    path "versions.yml",                             emit: versions

    script:
    """
    mkdir -p nucfreq_out

    # Build region list from FASTA index
    samtools faidx ${assembly_fasta}
    awk '{print \$1":1-"\$2}' ${assembly_fasta}.fai > regions.txt

    nucfreq \\
        -b ${bam} \\
        -f regions.txt \\
        -o nucfreq_out/${meta.id}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nucfreq: \$(nucfreq --version 2>&1 | head -1 || echo 'nucfreq')
    END_VERSIONS
    """
}

workflow NUCFREQ {
    take:
    ch_input  // [ meta, assembly_fasta, reference_fasta ]

    main:
    NUCFREQ_ALIGN(ch_input)
    NUCFREQ_RUN(NUCFREQ_ALIGN.out.bam)

    emit:
    results  = NUCFREQ_RUN.out.results
    versions = NUCFREQ_RUN.out.versions
}
