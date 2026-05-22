/*
 * QC workflow — all quality-assessment tools + summary report.
 *
 * Input channel: [ meta, assembly_fasta, reference_fasta, meryl_db (or []) ]
 *   meta.type = 'contigs' | 'scaffolds'  — propagated to publishDir and outputs
 *
 * Each tool is gated by a params.run_* toggle.
 * Tools run identically on contigs and scaffolds; output lands in:
 *   results/<sample_id>/<type>/<tool>/
 *
 * A final QC_REPORT step aggregates all available outputs into a per-assembly
 * PDF (params.run_report = true by default).
 */

nextflow.enable.dsl = 2

include { COMPLEASM  } from '../modules/local/compleasm/main'
include { QUAST      } from '../modules/local/quast/main'
include { MERQURY    } from '../modules/local/merqury/main'
include { LAI        } from '../modules/local/lai/main'
include { KPLEXITY   } from '../modules/local/kplexity/main'
include { FLAGGER    } from '../modules/local/flagger/main'
include { NUCFREQ    } from '../modules/local/nucfreq/main'
include { NCHART     } from '../modules/local/nchart/main'
include { QC_REPORT  } from '../modules/local/report/main'

// Unique string key for joining channels across tools
def _key = { meta -> "${meta.id}__${meta.type}" }

// Attach a key to a [meta, file] channel for use with join()
def keyed = { ch -> ch.map { meta, f -> [_key(meta), f] } }

// After joining, replace null (no match) with the NO_FILE sentinel
def no_file = file('NO_FILE')


workflow QC {
    take:
    ch_input  // [ meta, assembly_fasta, reference_fasta, meryl_db ]

    main:
    ch_versions = Channel.empty()

    // ── Run QC tools ─────────────────────────────────────────────────────────

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

    // ── Build report input channel ────────────────────────────────────────────
    if (params.run_report) {
        // Base: [key, meta, assembly_fasta]
        ch_base = ch_input.map { meta, fa, ref, meryl -> [ _key(meta), meta, fa ] }

        // Optional tool outputs as [key, file/dir] — Channel.empty() when tool disabled
        ch_quast_tsv   = params.run_quast     ? keyed(QUAST.out.report_tsv)            : Channel.empty()
        ch_mqv         = params.run_merqury   ? keyed(MERQURY.out.qv)                  : Channel.empty()
        ch_mqc         = params.run_merqury   ? keyed(MERQURY.out.completeness)         : Channel.empty()
        ch_mqdir       = params.run_merqury   ? keyed(MERQURY.out.results)              : Channel.empty()
        ch_compl       = params.run_compleasm ? keyed(COMPLEASM.out.summary)            : Channel.empty()
        ch_lai         = params.run_lai       ? keyed(LAI.out.lai_output)               : Channel.empty()
        ch_kplex       = params.run_kplexity  ? keyed(KPLEXITY.out.csv)                : Channel.empty()
        ch_flagger_bed = params.run_flagger   ? keyed(FLAGGER.out.bed)                 : Channel.empty()
        ch_nfdir       = params.run_nucfreq   ? keyed(NUCFREQ.out.results)              : Channel.empty()
        ch_ncdir       = params.run_nchart    ? keyed(NCHART.out.plot.map { m, p -> [m, p instanceof List ? p[0] : p] })
                                             : Channel.empty()

        // Join all optional outputs to the base channel.
        // remainder: true keeps base entries that have no matching tool output.
        // Nulls (unmatched) are replaced with the NO_FILE sentinel.
        ch_report = ch_base
            .join(ch_quast_tsv,   remainder: true)
            .join(ch_mqv,         remainder: true)
            .join(ch_mqc,         remainder: true)
            .join(ch_mqdir,       remainder: true)
            .join(ch_compl,       remainder: true)
            .join(ch_lai,         remainder: true)
            .join(ch_kplex,       remainder: true)
            .join(ch_flagger_bed, remainder: true)
            .join(ch_nfdir,       remainder: true)
            .join(ch_ncdir,       remainder: true)
            .map { key, meta, fa, quast, mqv, mqc, mqdir, compl, lai, kplex, fla, nfdir, ncdir ->
                [
                    meta,
                    fa,
                    quast  ?: no_file,
                    mqv    ?: no_file,
                    mqc    ?: no_file,
                    mqdir  ?: no_file,
                    compl  ?: no_file,
                    lai    ?: no_file,
                    kplex  ?: no_file,
                    fla    ?: no_file,
                    nfdir  ?: no_file,
                    ncdir  ?: no_file,
                ]
            }

        QC_REPORT(ch_report)
        ch_versions = ch_versions.mix(QC_REPORT.out.versions)
    }

    emit:
    versions = ch_versions
}
