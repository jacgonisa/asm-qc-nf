nextflow.enable.dsl = 2

/*
 * Nchart — contig/scaffold size distribution chart (N50 visual).
 * Tool: fastalengths binary → Nchart_script.R
 *
 * Requires:
 *   params.fastalengths_bin  — path to fastalengths binary
 *   params.nchart_script     — path to Nchart_script.R
 */

process NCHART {
    tag "${meta.id}"
    label 'process_single'

    container 'quay.io/biocontainers/r-base:4.3.1'

    publishDir "${params.outdir}/${meta.id}/${meta.type}/nchart", mode: 'symlink'

    input:
    tuple val(meta), path(assembly_fasta), path(reference_fasta)

    output:
    tuple val(meta), path("*.lengths"),   emit: lengths
    tuple val(meta), path("*.png"),       emit: plot, optional: true
    path "versions.yml",                  emit: versions

    script:
    def fastalengths = params.fastalengths_bin
    def nchart_r     = params.nchart_script
    if (!fastalengths) error "params.fastalengths_bin is required for NCHART"
    if (!nchart_r)     error "params.nchart_script is required for NCHART"
    """
    chmod +x ${fastalengths}

    ${fastalengths} ${assembly_fasta}  > ${meta.id}.assembly.lengths
    ${fastalengths} ${reference_fasta} > ${meta.id}.reference.lengths

    Rscript - <<'RSCRIPT'
    # Minimal inline Nchart: read lengths files, compute cumulative size, plot
    asm_len <- sort(read.table("${meta.id}.assembly.lengths", sep="\t")\$V2, decreasing=TRUE)
    ref_len <- sort(read.table("${meta.id}.reference.lengths", sep="\t")\$V2, decreasing=TRUE)

    png("${meta.id}.nchart.png", width=900, height=600, res=120)
    plot(cumsum(asm_len)/1e6, type="l", col="#2980b9", lwd=2,
         xlab="Contig rank", ylab="Cumulative size (Mb)",
         main="${meta.id} — Nchart", ylim=c(0, max(cumsum(asm_len),cumsum(ref_len))/1e6))
    lines(cumsum(ref_len)/1e6, col="#228833", lwd=2, lty=2)
    legend("bottomright", c("Assembly","Reference"), col=c("#2980b9","#228833"),
           lwd=2, lty=c(1,2), bty="n")
    dev.off()
    RSCRIPT

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nchart: fastalengths + R inline
    END_VERSIONS
    """
}
