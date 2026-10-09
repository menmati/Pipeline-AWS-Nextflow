#!/usr/bin/env nextflow
/*
 * ============================================================
 *  RNA-seq Pipeline – miR-150-5p Inhibition in Mouse Cortex
 * ============================================================
 * Organism : Mus musculus (C57BL/6J), Embryonic day 18
 * Project  : PRJNA1533324
 * Aligner  : HISAT2  |  Quantification : featureCounts
 * DE       : DESeq2  |  Contrast : miR150_inhibition vs negative_control
 * ============================================================
 */
nextflow.enable.dsl = 2

include { FASTQC as FASTQC_RAW     } from './modules/fastqc'
include { FASTQC as FASTQC_TRIMMED } from './modules/fastqc'
include { TRIMMOMATIC              } from './modules/trimmomatic'
include { HISAT2_ALIGN  } from './modules/hisat2'
include { FEATURECOUNTS } from './modules/featurecounts'
include { DESEQ2        } from './modules/deseq2'

// ── Parameter defaults ──────────────────────────────────
params.samplesheet          = "${projectDir}/data/samplesheet.csv"
params.genome               = "${projectDir}/genome/genome.fa"
params.gtf                  = "${projectDir}/genome/annotation.gtf"
params.hisat2_index         = "${projectDir}/genome/hisat2_index"
params.hisat2_index_base    = "GRCm39"
params.trimmomatic_adapters = "TruSeq3-PE-2.fa"
params.outdir               = "${projectDir}/results"
params.contrast             = ["condition", "miR150_inhibition", "negative_control"]
params.r_script             = "${projectDir}/scripts/deseq2_analysis.R"

// ── Workflow ─────────────────────────────────────────────
workflow {

    // 1. Parse samplesheet → channel of [meta, [r1, r2]]
    Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row ->
            def meta = [id: row.sample, condition: row.condition]
            [meta, [file(row.fastq_1), file(row.fastq_2)]]
        }
        .set { reads_ch }

    // 2. QC on raw reads
    FASTQC_RAW(reads_ch)

    // 3. Adapter trimming
    TRIMMOMATIC(reads_ch)

    // 4. QC on trimmed reads
    FASTQC_TRIMMED(
        TRIMMOMATIC.out.trimmed_reads.map { meta, r1, r2 -> [meta, [r1, r2]] }
    )

    // 5. Alignment (trimmed reads)
    HISAT2_ALIGN(
        TRIMMOMATIC.out.trimmed_reads,
        file(params.hisat2_index),
        params.hisat2_index_base
    )

    // 6. Collect all BAMs for a single featureCounts call
    bam_ch = HISAT2_ALIGN.out.bam.map { meta, bam -> bam }.collect()
    bai_ch = HISAT2_ALIGN.out.bai.map { meta, bai -> bai }.collect()

    // 7. Gene-level quantification
    FEATURECOUNTS(bam_ch, bai_ch, file(params.gtf))

    // 8. Differential expression
    DESEQ2(
        FEATURECOUNTS.out.counts,
        file(params.samplesheet),
        params.contrast,
        file(params.r_script)
    )
}
