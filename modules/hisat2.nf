/*
 * HISAT2 – Splice-aware alignment to GRCm39.
 * Pipes directly into SAMtools for sort + index.
 */
process HISAT2_ALIGN {
    tag { "${meta.id}" }
    label 'process_high'

    conda 'bioconda::hisat2=2.2.1 bioconda::samtools=1.21'
    container 'quay.io/biocontainers/mulled-v2-a97e90b3b802d1da3d6958e0867610c718cb5eb1:2388ff74dec9a0e006acdb57ef8a1f7fcba8aeab-0'

    publishDir { "${params.outdir}/hisat2/${meta.id}" }, mode: 'copy'

    input:
    tuple val(meta), path(r1), path(r2)
    path  index_dir
    val   index_base

    output:
    tuple val(meta), path("${meta.id}.bam"),     emit: bam
    tuple val(meta), path("${meta.id}.bam.bai"), emit: bai
    path "${meta.id}_hisat2.log",                emit: log
    path 'versions.yml',                         emit: versions

    script:
    """
    hisat2 \
        -x ${index_dir}/${index_base} \
        -1 ${r1} -2 ${r2} \
        --threads ${task.cpus} \
        --dta \
        --rna-strandness RF \
        --summary-file ${meta.id}_hisat2.log \
    | samtools sort -@ ${task.cpus} -o ${meta.id}.bam -

    samtools index ${meta.id}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hisat2: \$(hisat2 --version | head -1 | sed 's/.*HISAT2 version //')
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """
}
