#!/usr/bin/env nextflow
/*
 * asm-qc-nf  —  Genome assembly quality-control pipeline
 *
 * Input (samplesheet CSV):
 *   contigs_fasta   — raw assembler output (hifiasm .fa / .gfa-derived)
 *   scaffolds_fasta — reference-guided scaffolds (RagTag, etc.)  [optional]
 *   reference_fasta — reference genome for QUAST, Flagger, NucFreq
 *   meryl_db        — Meryl k-mer database for Merqury             [optional]
 *
 * At least one of contigs_fasta / scaffolds_fasta is required per row.
 * All QC tools run independently on contigs and scaffolds, labelled in output.
 *
 * Run:
 *   nextflow run main.nf -profile singularity \
 *     --input samplesheet.csv --outdir results \
 *     --kplexity_bin /path/to/Kplex \
 *     --ltr_retriever /path/to/LTR_retriever
 */

nextflow.enable.dsl = 2

include { QC } from './workflows/qc'

// ── Samplesheet parsing ───────────────────────────────────────────────────────
def parse_samplesheet(csv_path) {
    Channel
        .fromPath(csv_path, checkIfExists: true)
        .splitCsv(header: true, strip: true)
        .flatMap { row ->
            def base_meta = [
                id:            row.sample,
                busco_lineage: row.busco_lineage ?: params.busco_lineage,
            ]
            def ref   = row.reference_fasta ? file(row.reference_fasta, checkIfExists: true) : null
            def meryl = row.meryl_db        ? file(row.meryl_db,        checkIfExists: true) : []

            if (!ref) error "Sample '${base_meta.id}': reference_fasta is required"

            def entries = []

            if (row.contigs_fasta?.trim()) {
                def fa = file(row.contigs_fasta, checkIfExists: true)
                entries << [ base_meta + [type: 'contigs'],   fa, ref, meryl ]
            }
            if (row.scaffolds_fasta?.trim()) {
                def fa = file(row.scaffolds_fasta, checkIfExists: true)
                entries << [ base_meta + [type: 'scaffolds'], fa, ref, meryl ]
            }

            if (entries.isEmpty())
                error "Sample '${base_meta.id}': provide contigs_fasta and/or scaffolds_fasta"

            entries
        }
}

// ── Main ──────────────────────────────────────────────────────────────────────
workflow {
    ch_input = parse_samplesheet(params.input)
    QC(ch_input)
}
