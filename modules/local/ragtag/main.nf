nextflow.enable.dsl = 2

process RAGTAG {
    tag "${meta.id}"
    label 'process_medium'

    container 'quay.io/biocontainers/ragtag:2.1.0--pyhb7b1952_0'

    publishDir "${params.outdir}/${meta.id}/ragtag", mode: 'symlink'

    input:
    tuple val(meta), path(query_fasta), path(reference_fasta)

    output:
    tuple val(meta), path("${meta.id}.scaffold.chr_only.fasta"), emit: chr_fasta
    tuple val(meta), path("ragtag_output/"),                     emit: ragtag_dir
    path "versions.yml",                                         emit: versions

    script:
    def min_len = params.ragtag_min_len ?: 1000000
    """
    ragtag.py scaffold \\
        ${reference_fasta} \\
        ${query_fasta} \\
        -o ragtag_output \\
        -t ${task.cpus} \\
        --mm2-params "-x asm5"

    # Keep only chr-scale scaffolds (>= min_len)
    python3 - <<'PY'
from Bio import SeqIO
min_len = ${min_len}
with open("ragtag_output/ragtag.scaffold.fasta") as fh, \
     open("${meta.id}.scaffold.chr_only.fasta", "w") as out:
    for rec in SeqIO.parse(fh, "fasta"):
        if len(rec) >= min_len:
            SeqIO.write(rec, out, "fasta")
PY

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ragtag: \$(ragtag.py --version 2>&1 | head -1)
    END_VERSIONS
    """
}
