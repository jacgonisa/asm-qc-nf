/*
 * ASSEMBLE workflow: HiFi reads → hifiasm contigs → RagTag scaffolds
 */

nextflow.enable.dsl = 2

include { HIFIASM  } from '../modules/local/hifiasm/main'
include { RAGTAG   } from '../modules/local/ragtag/main'

workflow ASSEMBLE {
    take:
    ch_input  // [ meta, hifi_reads, pat_yak, mat_yak, reference_fasta ]

    main:
    // ── hifiasm ───────────────────────────────────────────────────────────────
    HIFIASM(ch_input)

    // ── RagTag scaffold ───────────────────────────────────────────────────────
    ch_ragtag_in = HIFIASM.out.primary_fasta
        .join(
            ch_input.map { meta, reads, pat, mat, ref -> [ meta, ref ] },
            by: 0
        )
    RAGTAG(ch_ragtag_in)

    emit:
    assembly     = RAGTAG.out.chr_fasta   // [ meta, chr_only.fasta ]
    contigs      = HIFIASM.out.primary_fasta
    versions     = HIFIASM.out.versions.mix(RAGTAG.out.versions)
}
