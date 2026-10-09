# RNA-seq Pipeline — miR-150-5p Inhibition in Mouse Embryonic Cortex

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A524.04-brightgreen)](https://www.nextflow.io/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

## Overview

End-to-end RNA-seq differential expression analysis for **PRJNA1533324**,
a prenatal alcohol exposure study examining the effect of miR-150-5p inhibition
in *Mus musculus* (C57BL/6J) embryonic cortex (E18).

**Contrast:** `miR-150-5p inhibition` vs. `negative control inhibition`

---

## Pipeline Architecture

```
SRA (fasterq-dump + seqtk)
         │  1M read pairs / sample
         ▼
    [FastQC] ──────────────────────────► raw QC reports
         │
         ▼
    [Trimmomatic] ── TruSeq3-PE-2 adapters, MINLEN:36
         │
         ▼
    [FastQC] ──────────────────────────► post-trim QC reports
         │
         ▼
    [HISAT2] ── splice-aware align to GRCm39/GENCODE M36
         │
         ▼
    [SAMtools] ── sort + index (BAM)
         │
         ▼
    [featureCounts] ── gene-level counts (reverse-stranded)
         │  counts.txt
         ▼
    [DESeq2 R script] ─────────────────► deseq2_results.csv
                                          volcano_plot.pdf
                                          pca_plot.pdf
                                          normalized_counts.csv
```

---

## Samples

| SRR Accession | Condition | GEO |
|---------------|-----------|-----|
| SRR40794928 | Negative Control | GSM10069104 |
| SRR40794929 | Negative Control | GSM10069103 |
| SRR40794930 | Negative Control | GSM10069102 |
| SRR40794931 | miR-150-5p Inhibition | GSM10069101 |
| SRR40794932 | miR-150-5p Inhibition | GSM10069100 |
| SRR40794933 | miR-150-5p Inhibition | GSM10069099 |

---

## Project Structure

```
rnaseq_pipeline/
├── main.nf                    # Nextflow DSL2 main workflow
├── nextflow.config            # Configuration (local / docker / aws profiles)
├── environment.yml            # Conda environment specification
├── modules/
│   ├── fastqc.nf              # Quality control
│   ├── trimmomatic.nf         # Adapter trimming
│   ├── hisat2.nf              # Alignment + BAM indexing
│   ├── featurecounts.nf       # Gene-level quantification
│   └── deseq2.nf              # Differential expression
├── scripts/
│   └── deseq2_analysis.R      # DESeq2 R script (volcano, PCA, results)
├── data/
│   ├── fastq/                 # Subsampled FASTQ files
│   └── samplesheet.csv        # Sample manifest
├── genome/
│   ├── genome.fa              # GRCm39 genome (GENCODE M36)
│   ├── annotation.gtf         # GENCODE M36 annotation
│   └── hisat2_index/          # HISAT2 genome index
└── results/                   # All pipeline outputs
```

---

## Installation

### Prerequisites
- [Miniconda](https://docs.conda.io/en/latest/miniconda.html)
- [Nextflow](https://www.nextflow.io/docs/latest/install.html) (≥ 24.04)
- Docker Desktop (for `docker` profile)
- AWS CLI + account (for `aws` profile)

### Conda environment

```bash
conda env create -f environment.yml
conda activate rnaseq
```

### Install Nextflow (inside the conda env)

```bash
conda activate rnaseq
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/
```

---

## Usage

### Local — conda (development)

```bash
conda activate rnaseq
nextflow run main.nf -profile local
```

### Local — Docker (reproducible)

```bash
nextflow run main.nf -profile docker
```

### Cloud — AWS Batch

> Set the placeholders in `nextflow.config` first.

```bash
# 1. Upload reference data and FASTQs to S3
aws s3 sync genome/ s3://YOUR-BUCKET-NAME/genome/
aws s3 sync data/   s3://YOUR-BUCKET-NAME/data/
aws s3 cp  scripts/deseq2_analysis.R s3://YOUR-BUCKET-NAME/scripts/

# 2. Run pipeline on AWS Batch
nextflow run main.nf -profile aws
```

### Override parameters

```bash
nextflow run main.nf -profile local \
    --outdir ./my_results \
    --alpha 0.01
```

---

## Outputs

| Path | Description |
|------|-------------|
| `results/fastqc/` | FastQC HTML + ZIP reports |
| `results/trimmomatic/` | Trimmomatic logs |
| `results/hisat2/` | BAM + BAI files per sample |
| `results/featurecounts/counts.txt` | Raw gene count matrix |
| `results/deseq2/deseq2_results.csv` | DE results (log2FC, padj, direction) |
| `results/deseq2/volcano_plot.pdf` | Volcano plot (top 15 DEGs labelled) |
| `results/deseq2/pca_plot.pdf` | PCA of VST-normalised counts |
| `results/deseq2/normalized_counts.csv` | DESeq2 size-factor normalised counts |
| `results/pipeline_report.html` | Nextflow execution report |
| `results/timeline.html` | Nextflow timeline report |

---

## AWS Scaling Architecture

```
Local Machine                         AWS Cloud
─────────────────                     ──────────────────────────────────
Nextflow (orchestrator) ──────────►  AWS Batch (job queue)
                                          │
                              ┌───────────┼───────────┐
                              ▼           ▼           ▼
                          EC2 (FastQC) EC2 (HISAT2) EC2 (DESeq2)
                              └───────────┴───────────┘
                                          │
                                          ▼
                                    Amazon S3 (storage)
                                    ├── genome/
                                    ├── data/fastq/
                                    └── results/
```

**Key advantage:** Switching from local to cloud requires only `-profile aws`.
Nextflow handles job scheduling, data staging from S3, and parallel execution
across multiple EC2 instances automatically.

---

## License

MIT
