/*
 * QC workflow: all quality-assessment tools run on a finished assembly FASTA.
 *
 * Input channel: [ meta, assembly_fasta, reference_fasta, meryl_db (or []) ]
 *
 * Tools (each gated by a params toggle):
 *   compleasm  — BUSCO-based completeness
 *   quast      — assembly statistics vs. reference
 *   merqury    — k-mer QV and completeness
 *   lai        — LTR Assembly Index (ltrharvest → ltrdigest → LTR_retriever)
 *   kplexity   — k-plexity complexity curve (k=5..151)
 *   flagger    — read-alignment based error flagging
 *   nucfreq    — nucleotide frequency (heterozygosity / collapse)
 *   nchart     — contig-size Nchart plot
 */

nextflow.enable.dsl = 2

include { COMPLEASM     } from '../modules/local/compleasm/main'
include { QUAST         } from '../modules/local/quast/main'
include { MERQURY       } from '../modules/local/merqury/main'
include { LAI           } from '../modules/local/lai/main'
include { KPLEXITY      } from '../modules/local/kplexity/main'
include { FLAGGER       } from '../modules/local/flagger/main'
include { NUCFREQ       } from '../modules/local/nucfreq/main'
include { NCHART        } from '../modules/local/nchart/main'

workflow QC {
    take:
    ch_input  // [ meta, assembly_fasta, reference_fasta, meryl_db ]

    main:
    ch_versions = Channel.empty()

    // ── compleasm ─────────────────────────────────────────────────────────────
    if (params.run_compleasm) {
        COMPLEASM(
            ch_input.map { meta, asm, ref, meryl -> [ meta, asm ] }
        )
        ch_versions = ch_versions.mix(COMPLEASM.out.versions)
    }

    // ── QUAST ─────────────────────────────────────────────────────────────────
    if (params.run_quast) {
        QUAST(
            ch_input.map { meta, asm, ref, meryl -> [ meta, asm, ref ] }
        )
        ch_versions = ch_versions.mix(QUAST.out.versions)
    }

    // ── Merqury ───────────────────────────────────────────────────────────────
    if (params.run_merqury) {
        ch_merqury = ch_input.filter { meta, asm, ref, meryl -> meryl as Boolean }
            .map { meta, asm, ref, meryl -> [ meta, meryl, asm ] }
        MERQURY(ch_merqury)
        ch_versions = ch_versions.mix(MERQURY.out.versions)
    }

    // ── LAI ───────────────────────────────────────────────────────────────────
    if (params.run_lai) {
        LAI(
            ch_input.map { meta, asm, ref, meryl -> [ meta, asm ] }
        )
        ch_versions = ch_versions.mix(LAI.out.versions)
    }

    // ── kplexity ──────────────────────────────────────────────────────────────
    if (params.run_kplexity) {
        KPLEXITY(
            ch_input.map { meta, asm, ref, meryl -> [ meta, asm ] }
        )
        ch_versions = ch_versions.mix(KPLEXITY.out.versions)
    }

    // ── Flagger ───────────────────────────────────────────────────────────────
    if (params.run_flagger) {
        FLAGGER(
            ch_input.map { meta, asm, ref, meryl -> [ meta, asm, ref ] }
        )
        ch_versions = ch_versions.mix(FLAGGER.out.versions)
    }

    // ── NucFreq ───────────────────────────────────────────────────────────────
    if (params.run_nucfreq) {
        NUCFREQ(
            ch_input.map { meta, asm, ref, meryl -> [ meta, asm, ref ] }
        )
        ch_versions = ch_versions.mix(NUCFREQ.out.versions)
    }

    // ── Nchart ────────────────────────────────────────────────────────────────
    if (params.run_nchart) {
        NCHART(
            ch_input.map { meta, asm, ref, meryl -> [ meta, asm, ref ] }
        )
        ch_versions = ch_versions.mix(NCHART.out.versions)
    }

    emit:
    versions = ch_versions
}
