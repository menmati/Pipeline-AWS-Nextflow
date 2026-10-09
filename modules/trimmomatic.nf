/*
 * Trimmomatic – Adapter trimming and quality filtering for paired-end reads.
 * Uses TruSeq3-PE-2 adapters (standard Illumina NovaSeq).
 */
process TRIMMOMATIC {
    tag { "${meta.id}" }
    label 'process_medium'

    conda 'bioconda::trimmomatic=0.39'
    container 'staphb/trimmomatic:0.39'

    publishDir { "${params.outdir}/trimmomatic/${meta.id}" }, mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}_R1_trimmed.fastq.gz"),
                     path("${meta.id}_R2_trimmed.fastq.gz"),   emit: trimmed_reads
    tuple val(meta), path("${meta.id}_R1_unpaired.fastq.gz"),
                     path("${meta.id}_R2_unpaired.fastq.gz"),  emit: unpaired_reads
    path "${meta.id}_trimmomatic.log",                          emit: log
    path 'versions.yml',                                        emit: versions

    script:
    def (r1, r2) = reads
    """
    trimmomatic PE \
        -threads ${task.cpus} \
        -phred33 \
        ${r1} ${r2} \
        ${meta.id}_R1_trimmed.fastq.gz  ${meta.id}_R1_unpaired.fastq.gz \
        ${meta.id}_R2_trimmed.fastq.gz  ${meta.id}_R2_unpaired.fastq.gz \
        ILLUMINACLIP:${params.trimmomatic_adapters}:2:30:10:2:keepBothReads \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 \
        2> ${meta.id}_trimmomatic.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trimmomatic: \$(trimmomatic -version)
    END_VERSIONS
    """
}
