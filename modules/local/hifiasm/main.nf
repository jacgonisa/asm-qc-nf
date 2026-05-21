nextflow.enable.dsl = 2

process HIFIASM {
    tag "${meta.id}"
    label 'process_high'

    container 'quay.io/biocontainers/hifiasm:0.19.8--h43eeafb_0'

    publishDir "${params.outdir}/${meta.id}/hifiasm", mode: 'symlink'

    input:
    tuple val(meta), path(hifi_reads), path(pat_yak), path(mat_yak), path(reference)

    output:
    tuple val(meta), path("${meta.id}.primary.fa"),  emit: primary_fasta
    tuple val(meta), path("${meta.id}.hap1.fa"),     emit: hap1_fasta,  optional: true
    tuple val(meta), path("${meta.id}.hap2.fa"),     emit: hap2_fasta,  optional: true
    tuple val(meta), path("*.gfa"),                  emit: gfa
    path "versions.yml",                             emit: versions

    script:
    def trio = (pat_yak && mat_yak && pat_yak.name != 'NO_FILE' && mat_yak.name != 'NO_FILE')
    def extra = params.hifiasm_extra ?: ''
    if (trio)
        """
        hifiasm -o ${meta.id} \\
            -t ${task.cpus} \\
            --h1 ${pat_yak} --h2 ${mat_yak} \\
            ${extra} \\
            ${hifi_reads}

        # Convert trio hap GFAs to FASTA
        awk '/^S/{print ">"\$2"\\n"\$3}' ${meta.id}.dip.hap1.p_ctg.gfa > ${meta.id}.hap1.fa
        awk '/^S/{print ">"\$2"\\n"\$3}' ${meta.id}.dip.hap2.p_ctg.gfa > ${meta.id}.hap2.fa
        # Use hap1 as primary in trio mode
        cp ${meta.id}.hap1.fa ${meta.id}.primary.fa

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            hifiasm: \$(hifiasm --version 2>&1 | head -1)
        END_VERSIONS
        """
    else
        """
        hifiasm -o ${meta.id} \\
            -t ${task.cpus} \\
            ${extra} \\
            ${hifi_reads}

        awk '/^S/{print ">"\$2"\\n"\$3}' ${meta.id}.bp.p_ctg.gfa > ${meta.id}.primary.fa

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            hifiasm: \$(hifiasm --version 2>&1 | head -1)
        END_VERSIONS
        """
}
