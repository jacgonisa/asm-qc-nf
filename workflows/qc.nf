/*
 * QC workflow — all quality-assessment tools.
 *
 * Input channel: [ meta, assembly_fasta, reference_fasta, meryl_db (or []) ]
 *   meta.type = 'contigs' | 'scaffolds'  — propagated to publishDir and outputs
 *
 * Each tool is gated by a params.run_* toggle.
 * Tools run identically on contigs and scaffolds; output lands in:
 *   results/<sample_id>/<type>/<tool>/
 */

nextflow.enable.dsl = 2

include { COMPLEASM } from '../modules/local/compleasm/main'
include { QUAST     } from '../modules/local/quast/main'
include { MERQURY   } from '../modules/local/merqury/main'
include { LAI       } from '../modules/local/lai/main'
include { KPLEXITY  } from '../modules/local/kplexity/main'
include { FLAGGER   } from '../modules/local/flagger/main'
include { NUCFREQ   } from '../modules/local/nucfreq/main'
include { NCHART    } from '../modules/local/nchart/main'

workflow QC {
    take:
    ch_input  // [ meta, assembly_fasta, reference_fasta, meryl_db ]

    main:
    ch_versions = Channel.empty()

    if (params.run_compleasm) {
        COMPLEASM(ch_input.map { meta, fa, ref, meryl -> [ meta, fa ] })
        ch_versions = ch_versions.mix(COMPLEASM.out.versions)
    }

    if (params.run_quast) {
        QUAST(ch_input.map { meta, fa, ref, meryl -> [ meta, fa, ref ] })
        ch_versions = ch_versions.mix(QUAST.out.versions)
    }

    if (params.run_merqury) {
        MERQURY(
            ch_input
                .filter { meta, fa, ref, meryl -> meryl as Boolean }
                .map    { meta, fa, ref, meryl -> [ meta, meryl, fa ] }
        )
        ch_versions = ch_versions.mix(MERQURY.out.versions)
    }

    if (params.run_lai) {
        LAI(ch_input.map { meta, fa, ref, meryl -> [ meta, fa ] })
        ch_versions = ch_versions.mix(LAI.out.versions)
    }

    if (params.run_kplexity) {
        KPLEXITY(ch_input.map { meta, fa, ref, meryl -> [ meta, fa ] })
        ch_versions = ch_versions.mix(KPLEXITY.out.versions)
    }

    if (params.run_flagger) {
        FLAGGER(ch_input.map { meta, fa, ref, meryl -> [ meta, fa, ref ] })
        ch_versions = ch_versions.mix(FLAGGER.out.versions)
    }

    if (params.run_nucfreq) {
        NUCFREQ(ch_input.map { meta, fa, ref, meryl -> [ meta, fa, ref ] })
        ch_versions = ch_versions.mix(NUCFREQ.out.versions)
    }

    if (params.run_nchart) {
        NCHART(ch_input.map { meta, fa, ref, meryl -> [ meta, fa, ref ] })
        ch_versions = ch_versions.mix(NCHART.out.versions)
    }

    emit:
    versions = ch_versions
}
