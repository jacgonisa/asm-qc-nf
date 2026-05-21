nextflow.enable.dsl = 2

/*
 * LAI (LTR Assembly Index) — three sequential steps:
 *   1. gt suffixerator  — build index
 *   2. gt ltrharvest    — find LTR-RT candidates
 *   3. gt ltrdigest     — annotate domains (optional tRNA HMM)
 *   4. LTR_retriever    — build NR library, compute LAI score
 *
 * Requires local installations:
 *   - GenomeTools  (params.genometools_bin, default: 'gt')
 *   - LTR_retriever (params.ltr_retriever)
 *
 * ltrharvest is slow (~60-120 min for 130 Mb Arabidopsis) and single-threaded.
 * Use label 'process_long' with a single CPU for the harvest step.
 */

process LAI_HARVEST {
    tag "${meta.id}"
    label 'process_long'

    publishDir "${params.outdir}/${meta.id}/${meta.type}/lai", mode: 'symlink', saveAs: { fn ->
        fn.endsWith('.harvest.scn') ? fn : null
    }

    input:
    tuple val(meta), path(assembly_fasta)

    output:
    tuple val(meta), path(assembly_fasta), path("${meta.id}.harvest.scn"), emit: harvest
    path "versions.yml",                                                    emit: versions

    script:
    def gt = params.genometools_bin ?: 'gt'
    """
    # Build suffix array index
    ${gt} suffixerator \\
        -db ${assembly_fasta} \\
        -indexname ${meta.id} \\
        -tis -suf -lcp -des -ssp -sds -dna \\
        2> suffixerator.log

    # Find LTR retrotransposon candidates
    ${gt} ltrharvest \\
        -index ${meta.id} \\
        -minlenltr 100 \\
        -maxlenltr 7000 \\
        -mintsd 4 \\
        -maxtsd 6 \\
        -motif TGCA \\
        -motifmis 1 \\
        -similar 85 \\
        -vic 10 \\
        -seed 20 \\
        -seqids yes \\
        > ${meta.id}.harvest.scn \\
        2> ltrharvest.log

    echo "[ltrharvest] candidates: \$(grep -v '^#' ${meta.id}.harvest.scn | wc -l)"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genometools: \$(${gt} --version 2>&1 | head -1)
    END_VERSIONS
    """
}

process LAI_DIGEST {
    tag "${meta.id}"
    label 'process_low'

    publishDir "${params.outdir}/${meta.id}/${meta.type}/lai", mode: 'symlink', saveAs: { fn ->
        fn.endsWith('.harvest.digest.scn') ? fn : null
    }

    input:
    tuple val(meta), path(assembly_fasta), path(harvest_scn)

    output:
    tuple val(meta), path(assembly_fasta), path(harvest_scn), path("${meta.id}.harvest.digest.scn"), emit: digest
    path "versions.yml",                                                                              emit: versions

    script:
    def gt = params.genometools_bin ?: 'gt'
    """
    ${gt} ltrdigest \\
        -trnas \$(find \$CONDA_PREFIX -name 'codon.hmm' 2>/dev/null | head -1 || echo '') \\
        -seqfile ${assembly_fasta} \\
        -matchdescstart \\
        ${harvest_scn} \\
        ${meta.id} \\
        > ${meta.id}.harvest.digest.scn \\
        2> ltrdigest.log \\
    || {
        echo "[WARN] ltrdigest failed (no tRNA HMM?), using harvest.scn directly"
        cp ${harvest_scn} ${meta.id}.harvest.digest.scn
    }

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genometools: \$(${gt} --version 2>&1 | head -1)
    END_VERSIONS
    """
}

process LAI_RETRIEVER {
    tag "${meta.id}"
    label 'process_high'

    publishDir "${params.outdir}/${meta.id}/${meta.type}/lai", mode: 'symlink'

    input:
    tuple val(meta), path(assembly_fasta), path(harvest_scn), path(digest_scn)

    output:
    tuple val(meta), path("${assembly_fasta}.LAI.output"), emit: lai_output, optional: true
    tuple val(meta), path("*.LAI.output"),                 emit: lai_output_alt, optional: true
    tuple val(meta), path("*.pass.list"),                  emit: pass_list, optional: true
    path "versions.yml",                                   emit: versions

    script:
    def ltr_ret = params.ltr_retriever ?: 'LTR_retriever'
    """
    ${ltr_ret} \\
        -genome ${assembly_fasta} \\
        -inharvest ${harvest_scn} \\
        -threads ${task.cpus} \\
        -u 1.3e-8 \\
        2>&1 | tee ltr_retriever.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        LTR_retriever: \$(${ltr_ret} -h 2>&1 | grep -i version | head -1 || echo 'LTR_retriever')
    END_VERSIONS
    """
}

workflow LAI {
    take:
    ch_input  // [ meta, assembly_fasta ]

    main:
    LAI_HARVEST(ch_input)
    LAI_DIGEST(LAI_HARVEST.out.harvest)
    LAI_RETRIEVER(LAI_DIGEST.out.digest)

    emit:
    lai_output = LAI_RETRIEVER.out.lai_output.mix(LAI_RETRIEVER.out.lai_output_alt)
    versions   = LAI_HARVEST.out.versions.mix(LAI_DIGEST.out.versions).mix(LAI_RETRIEVER.out.versions)
}
