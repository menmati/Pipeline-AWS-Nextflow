/*
 * featureCounts – Gene-level quantification.
 * All BAMs collected into one call → single count matrix output.
 * -s 2 = reverse-stranded (typical for dUTP/ligation libraries on NovaSeq).
 */
process FEATURECOUNTS {
    label 'process_medium'

    conda 'bioconda::subread=2.0.6'
    container 'quay.io/biocontainers/subread:2.0.6--he4a0461_0'

    publishDir "${params.outdir}/featurecounts", mode: 'copy'

    input:
    path bam_files
    path bai_files
    path gtf

    output:
    path 'counts.txt',         emit: counts
    path 'counts.txt.summary', emit: summary
    path 'versions.yml',       emit: versions

    script:
    """
    featureCounts \
        -T ${task.cpus} \
        -p \
        --countReadPairs \
        -s 2 \
        -t exon \
        -g gene_id \
        -a ${gtf} \
        -o counts.txt \
        ${bam_files}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        subread: \$(featureCounts -v 2>&1 | grep -oP '(?<=v)\\d+\\.\\d+\\.\\d+')
    END_VERSIONS
    """
}
