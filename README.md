# asm-qc-nf

Nextflow DSL2 pipeline for genome assembly quality control.

The user provides **contigs** and/or **scaffolds** (already assembled); the pipeline runs all QC tools on each. Contigs and scaffolds are treated independently and outputs are labelled accordingly — no assembler is included.

---

## Quick start

```bash
nextflow run main.nf \
  -profile singularity \
  --input assets/samplesheet.csv \
  --outdir results \
  --kplexity_bin /path/to/FASTK/Kplex \
  --ltr_retriever /path/to/LTR_retriever/LTR_retriever
```

---

## Samplesheet format

`assets/samplesheet.csv`:

| column | required | description |
|---|---|---|
| `sample` | yes | unique sample ID |
| `contigs_fasta` | one required | raw assembler contigs (hifiasm .fa) |
| `scaffolds_fasta` | one required | reference-guided scaffolds (RagTag, etc.) |
| `reference_fasta` | yes | reference genome for QUAST, alignment tools |
| `meryl_db` | for Merqury | path to Meryl k-mer database (see below) |
| `busco_lineage` | no | overrides `params.busco_lineage` per sample |

At least one of `contigs_fasta` or `scaffolds_fasta` must be provided. When both are given, all tools run on each independently and outputs are labelled `contigs/` or `scaffolds/` in the results directory.

---

## Building a Meryl database (required for Merqury)

Merqury computes k-mer–based QV and completeness scores. It requires a Meryl
database built from the **same HiFi reads** used for assembly (or from the reference,
for reference-based QV).

### From HiFi reads (recommended)

```bash
conda activate merqury   # or: module load merqury

# Single reads file
meryl k=21 count reads.hifi.fastq.gz output sample.meryl

# Multiple read files
meryl k=21 count *.fastq.gz output sample.meryl
```

Use `k=21` for most genomes. For very large or highly repetitive genomes, `k=19` or
`k=25` may be more appropriate.

### From a reference genome (alternative)

If you don't have the original reads but want reference-based completeness:

```bash
meryl k=21 count reference.fa output reference.meryl
```

### For trio Merqury (hapmers)

When assessing haplotype-phased assemblies, build haplotype-specific k-mer sets:

```bash
# Build individual databases
meryl k=21 count paternal_reads.fastq.gz output pat.meryl
meryl k=21 count maternal_reads.fastq.gz output mat.meryl

# Build hapmers (k-mers specific to each haplotype)
$MERQURY/trio/hapmers.sh pat.meryl mat.meryl child_reads.fastq.gz
```

This produces `pat.hapmer.meryl` and `mat.hapmer.meryl` — pass these as `meryl_db`
in the samplesheet.

### Quick reference: Meryl in the CENH3ox project

The Col-HiFi reference Meryl DB is at:
```
/mnt/ssd-4tb/HIFI_NAMIL/01_genomes/Col-HiFi/Col-0.ragtag_scaffolds.merylDB/
```

Build for a new sample (example — 40× parent reads):
```bash
meryl k=21 count /path/to/cenh3ox_col_parent_40x.fastq.gz \
    output /mnt/ssd-4tb/HIFI_NAMIL/reassembly/assembly_qc/cenh3ox_col_parent_40x/meryl/reads.meryl
```

---

## Parameters

### Required for tools without containers

| param | description |
|---|---|
| `--kplexity_bin` | Full path to `Kplex` binary (FastK must be in the same directory) |
| `--ltr_retriever` | Full path to `LTR_retriever` script |
| `--fastalengths_bin` | Full path to `fastalengths` binary (Nchart) |
| `--nchart_script` | Full path to `Nchart_script.R` |

### Key optional parameters

| param | default | description |
|---|---|---|
| `--busco_lineage` | `embryophyta_odb10` | BUSCO lineage for compleasm |
| `--kplexity_k` | `5:151:1` | k-mer range (start:end:step) |
| `--hifiasm_extra` | `''` | Extra args passed to hifiasm |
| `--ragtag_min_len` | `1000000` | Min scaffold length to keep (chr filter) |
| `--run_assembly` | `true` | Run hifiasm + RagTag if reads provided |
| `--run_compleasm` | `true` | Toggle each QC tool individually |
| `--run_merqury` | `true` | Requires `meryl_db` in samplesheet |
| `--run_lai` | `true` | Slow (~2 h per assembly); uses local GT + LTR_retriever |
| `--run_kplexity` | `true` | Requires `--kplexity_bin` |
| `--run_flagger` | `true` | |
| `--run_nucfreq` | `true` | |
| `--run_nchart` | `true` | Requires `--fastalengths_bin` + `--nchart_script` |

---

## Profiles

| profile | description |
|---|---|
| `singularity` | Singularity containers (recommended on HPC) |
| `docker` | Docker containers |
| `conda` | Conda environments (needed for LAI, Flagger) |
| `local` | Local executor |
| `slurm` | SLURM executor |

Combine profiles: `-profile singularity,slurm`

---

## Tool notes

### LAI (LTR Assembly Index)
LAI runs as three sequential Nextflow processes to avoid memory conflicts from
parallel ltrharvest instances:
1. `LAI_HARVEST` — suffixerator + ltrharvest (~60–120 min for ~130 Mb)
2. `LAI_DIGEST`  — ltrdigest (fast; falls back gracefully if no tRNA HMM)
3. `LAI_RETRIEVER` — LTR_retriever → final LAI score

Requires local installations in PATH or via `params.genometools_bin` / `params.ltr_retriever`.

### kplexity
Uses a custom FASTK-based binary. The FastK directory (containing `FastK`, `Fastrm`, etc.)
must be the same directory as `Kplex`, or on PATH. Set with `--kplexity_bin`.

Typical run: `kplexity_k = '5:151:1'` with `kplexity_h = '1:10000000'` (~8 min per
Arabidopsis-sized assembly with 16 threads).

### Nchart
A lightweight contig-size chart. Uses the `fastalengths` binary from the original
Nchart repo plus an inline R plot. No R packages beyond base R are required.

---

## Output structure

```
results/
└── <sample_id>/
    ├── contigs/          (present if contigs_fasta was provided)
    │   ├── compleasm/    summary.txt
    │   ├── quast/        report.tsv, report.html
    │   ├── merqury/      *.qv, *.completeness.stats
    │   ├── lai/          *.harvest.scn, *.LAI.output
    │   ├── kplexity/     *.kplex.csv
    │   ├── flagger/      *.bed
    │   ├── nucfreq/      *.png
    │   └── nchart/       *.lengths, *.nchart.png
    └── scaffolds/        (present if scaffolds_fasta was provided)
        ├── compleasm/
        ├── quast/
        └── ...           (same structure as contigs/)
```
