#!/usr/bin/env nextflow
/*
 * asm-qc-nf  —  HiFi genome assembly + quality-control pipeline
 *
 * Inputs (via samplesheet CSV):
 *   hifi_reads          → runs hifiasm + RagTag, then all QC
 *   hifi_reads + yak    → trio hifiasm (pat_yak / mat_yak)
 *   assembly_fasta only → QC-only mode (skip hifiasm/RagTag)
 *
 * Run:
 *   nextflow run main.nf -profile singularity --input samplesheet.csv --outdir results
 */

nextflow.enable.dsl = 2

include { ASSEMBLE } from './workflows/assemble'
include { QC       } from './workflows/qc'

// ── Samplesheet parsing ───────────────────────────────────────────────────────
def parse_samplesheet(csv_path) {
    Channel
        .fromPath(csv_path, checkIfExists: true)
        .splitCsv(header: true, strip: true)
        .map { row ->
            def meta = [
                id:             row.sample,
                busco_lineage:  row.busco_lineage  ?: params.busco_lineage,
            ]
            def hifi_reads = row.hifi_reads      ? file(row.hifi_reads,      checkIfExists: true) : []
            def pat_yak    = row.pat_yak          ? file(row.pat_yak,         checkIfExists: true) : []
            def mat_yak    = row.mat_yak          ? file(row.mat_yak,         checkIfExists: true) : []
            def asm_fasta  = row.assembly_fasta   ? file(row.assembly_fasta,  checkIfExists: true) : []
            def ref_fasta  = row.reference_fasta  ? file(row.reference_fasta, checkIfExists: true) : null
            def meryl_db   = row.meryl_db         ? file(row.meryl_db,        checkIfExists: true) : []

            if (!ref_fasta) error "Sample '${meta.id}': reference_fasta is required"
            if (!hifi_reads && !asm_fasta) error "Sample '${meta.id}': provide hifi_reads or assembly_fasta"

            [ meta, hifi_reads, pat_yak, mat_yak, asm_fasta, ref_fasta, meryl_db ]
        }
}

// ── Main ──────────────────────────────────────────────────────────────────────
workflow {

    ch_input = parse_samplesheet(params.input)

    // Split into assembly-mode and QC-only
    ch_needs_assembly = ch_input.filter { meta, reads, pat, mat, asm, ref, meryl ->
        reads && params.run_assembly
    }
    ch_qc_only = ch_input.filter { meta, reads, pat, mat, asm, ref, meryl ->
        asm as Boolean
    }

    // ── Assembly path ─────────────────────────────────────────────────────────
    ch_assembled = Channel.empty()
    if (params.run_assembly) {
        ASSEMBLE(
            ch_needs_assembly.map { meta, reads, pat, mat, asm, ref, meryl ->
                [ meta, reads, pat, mat, ref ]
            }
        )
        ch_assembled = ASSEMBLE.out.assembly
            .join(
                ch_needs_assembly.map { meta, reads, pat, mat, asm, ref, meryl ->
                    [ meta, ref, meryl ]
                },
                by: 0
            )
            .map { meta, asm, ref, meryl -> [ meta, asm, ref, meryl ] }
    }

    // ── QC-only path: combine with existing FASTAs ────────────────────────────
    ch_qc_direct = ch_qc_only.map { meta, reads, pat, mat, asm, ref, meryl ->
        [ meta, asm, ref, meryl ]
    }

    ch_all_for_qc = ch_assembled.mix(ch_qc_direct)

    // ── QC ────────────────────────────────────────────────────────────────────
    QC(ch_all_for_qc)
}
