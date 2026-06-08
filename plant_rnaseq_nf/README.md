# plant_rnaseq_nf

A **Nextflow DSL2** pipeline that takes **raw Illumina short reads** from two
plant genotypes — a stress-**tolerant** and a stress-**susceptible** line — all
the way to differential expression and a **phenotype-facing functional
summary**. It is the read-processing companion to
[`aml_rnaseq_nf`](../aml_rnaseq_nf): where that pipeline starts from
pre-quantified counts, this one demonstrates the full short-read path —
QC/trim → spliced genome alignment → gene quantification → DE → GO enrichment.

The default target is rice (*Oryza sativa*) under drought or salinity stress
(e.g. IR64 vs. Apo, or IR29 vs. Pokkali), but nothing is rice-specific: point it
at any genome, GTF, GMT, and samplesheet.

## Stages

| #  | Process            | Tool / script                | Output                          |
|----|--------------------|------------------------------|---------------------------------|
| 1  | `SUBSAMPLE`        | `bin/subsample_reads.sh` (seqtk) | downsampled FASTQ           |
| 2  | `QC_TRIM`          | `fastp` + FastQC             | trimmed FASTQ, QC reports       |
| 3  | `HISAT2_BUILD`     | `hisat2-build`              | genome index                    |
| 4  | `ALIGN`            | `hisat2` + `samtools`        | sorted, indexed BAM             |
| 5  | `QUANTIFY`         | `featureCounts`             | `featurecounts.txt`             |
| 6  | `BUILD_MATRIX`     | `bin/build_count_matrix.py`  | `counts_raw.tsv`, `metadata.tsv`|
| 7  | `RUN_DE`           | `bin/run_de.py` (pydeseq2)   | `de_results.tsv`, `norm_counts.tsv` |
| 8  | `ENRICH`           | `bin/run_enrichment.py`      | `enrichment.tsv`                |
| 9  | `MAKE_HEATMAP`     | `bin/make_heatmap.py`        | `heatmap.html`                  |
| 10 | `MAKE_ENRICH_PLT`  | `bin/make_enrichment_plot.py`| `enrichment.html`               |

All outputs land in `results/` via Nextflow `publishDir 'copy'`.

## Run

```bash
cd plant_rnaseq_nf
mamba env create -f envs/plant_rnaseq_env.yml   # or: conda env create -f ...
conda activate plant_rnaseq_env

# Offline demo on a synthetic toy genome + planted DE signal (seconds):
bash run_local.sh --demo
#   or with Nextflow:
nextflow run main.nf -profile conda \
    --samplesheet data/demo/samplesheet.csv --genome data/demo/genome.fa \
    --gtf data/demo/genes.gtf --gmt data/demo/go_sets.gmt
```

For a **real run**, fetch a rice genome + SRA runs once and point the pipeline
at them:

```bash
bash fetch_real_data.sh          # creates data/real/runs.tsv template first time
# ...fill runs.tsv from the SRA Run Selector, then re-run fetch_real_data.sh...
nextflow run main.nf -profile conda \
    --samplesheet data/real/samplesheet.csv --genome data/real/genome.fa \
    --gtf data/real/genes.gtf --gmt /path/to/rice_go.gmt --subsample 5000000
```

`--subsample N` deterministically downsamples each sample to `N` read pairs
(seed-controlled) so a full study runs on a laptop; `--subsample 0` uses all
reads.

## Parameters

Override at the Nextflow CLI (`--param value`) or in `nextflow.config`.

| Parameter      | Default | Meaning                                            |
|----------------|---------|----------------------------------------------------|
| `samplesheet`  | `data/samplesheet.csv` | CSV: `sample_id,genotype,fastq_1,fastq_2` |
| `genome`       | `data/genome.fa`       | reference FASTA                     |
| `gtf`          | `data/genes.gtf`       | gene annotation (`-t exon -g gene_id`) |
| `gmt`          | `data/go_sets.gmt`     | GO/pathway gene sets for enrichment |
| `subsample`    | `0`     | read pairs/sample (0 = all)                        |
| `seed`         | `42`    | subsampling RNG seed                               |
| `fdr`          | `0.05`  | DE padj cutoff (significance / enrichment input)   |
| `lfc`          | `1.0`   | DE \|log2FC\| cutoff                               |
| `topn`         | `40`    | top DE genes shown in the heatmap                  |
| `tolerant` / `susceptible` | `tolerant` / `susceptible` | genotype labels + contrast direction |

## Methods

- **QC / trim.** `fastp` removes adapters and low-quality tails; FastQC reports
  are published for inspection.
- **Alignment.** `HISAT2` (spliced, low-memory vs. STAR) against a genome index,
  sorted and indexed with `samtools`.
- **Quantification.** `featureCounts` over exons (`-g gene_id`), paired-end
  aware (`-p --countReadPairs`).
- **DE.** `pydeseq2` — DESeq2 median-of-ratios size factors and a
  negative-binomial Wald test on `~genotype`, contrasting tolerant vs.
  susceptible. Positive log2FC = higher in the tolerant line.
- **Enrichment.** One-sided hypergeometric over-representation of the
  significant DE genes against the GMT gene sets, with hand-rolled
  Benjamini–Hochberg adjustment (same approach as `aml_rnaseq_nf`). The enriched
  terms are the phenotype link: the biological processes that distinguish the
  two lines.
- **Visualization.** A row-z-scored, hierarchically clustered Plotly heatmap of
  the top DE genes, and a bar chart of the most enriched processes — both in the
  naraen.net "Phalaena Automata" palette.

See [`docs/REPORT.md`](docs/REPORT.md) for an end-to-end demo run.

## Layout

```text
plant_rnaseq_nf/
├── main.nf                  # Nextflow DSL2 workflow (10 processes)
├── nextflow.config          # manifest, params, profiles
├── fetch_real_data.sh       # rice genome + SRA fetch -> data/real/
├── run_local.sh             # mirrors main.nf without Nextflow (--demo offline)
├── bin/
│   ├── subsample_reads.sh
│   ├── make_demo_data.py    # toy genome + planted-DE reads for demo/CI
│   ├── build_count_matrix.py
│   ├── run_de.py            # pydeseq2
│   ├── run_enrichment.py    # hypergeometric ORA + BH FDR
│   ├── make_heatmap.py
│   └── make_enrichment_plot.py
├── envs/plant_rnaseq_env.yml
├── tests/test_pipeline.py   # fast unit tests (no aligner needed)
├── docs/REPORT.md
├── data/                    # inputs (gitignored)
└── results/                 # outputs (gitignored)
```

## Reproducibility

- Deterministic subsampling (single `seed`); all paths project-relative.
- `nextflow.config` sets `PYTHONNOUSERSITE=1` so processes don't leak user
  site-packages.
- Pinned conda env in `envs/plant_rnaseq_env.yml`.
- `run_local.sh --demo` runs the entire DAG offline on a synthetic genome — used
  as the CI smoke test and to keep `main.nf` and `run_local.sh` in sync.

## Notes / next steps

- Subsampling is a documented compute trade-off, not a power-optimized design.
- `STAR` is a drop-in alternative to HISAT2 where RAM allows.
- For real runs, supply a rice GO GMT (g:Profiler / PlantGSEA export); the
  `--gmt` format matches `data/demo/go_sets.gmt`.
