nextflow.enable.dsl = 2

/*
 * k-plexity curve (fraction_unique(k) for k = 5..151)
 *
 * Requires local Kplex binary: params.kplexity_bin
 * FastK must also be in PATH (same directory as Kplex, or set in nextflow.config).
 *
 * Output: CSV with columns k, total_kmers, unique_kmers, fraction_unique
 */

process KPLEXITY {
    tag "${meta.id}"
    label 'process_high'

    publishDir "${params.outdir}/${meta.id}/${meta.type}/kplexity", mode: 'symlink'

    input:
    tuple val(meta), path(assembly_fasta)

    output:
    tuple val(meta), path("${meta.id}.kplex.csv"), emit: csv
    path "versions.yml",                           emit: versions

    script:
    def kplex_bin = params.kplexity_bin
    if (!kplex_bin) error "params.kplexity_bin is required for KPLEXITY process"
    def kplex_dir = kplex_bin.replaceAll(/\/[^\/]+$/, '')   // parent dir for FastK
    def k_range   = params.kplexity_k ?: '5:151:1'
    def h_range   = params.kplexity_h ?: '1:10000000'
    """
    export PATH="${kplex_dir}:\${PATH}"

    ${kplex_bin} \\
        -i ${assembly_fasta} \\
        -k ${k_range} \\
        -h ${h_range} \\
        -T ${task.cpus} \\
        -o tmp_kplex_${meta.id} \\
        -c ${meta.id}.kplex.csv

    rm -f tmp_kplex_${meta.id}.*

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kplexity: \$(${kplex_bin} --version 2>&1 | head -1 || echo 'kplexity (custom)')
    END_VERSIONS
    """
}
