/*
 * DESeq2 – Differential expression analysis.
 * Calls scripts/deseq2_analysis.R and emits results + plots.
 */
process DESEQ2 {
    label 'process_medium'

    conda 'conda-forge::r-base=4.4.1 bioconda::bioconductor-deseq2=1.46.0 bioconda::bioconductor-apeglm conda-forge::r-ggplot2 conda-forge::r-ggrepel conda-forge::r-pheatmap conda-forge::r-rcolorbrewer conda-forge::r-optparse'
    container 'ghcr.io/rocker-project/rocker/verse:4.4.1'

    publishDir "${params.outdir}/deseq2", mode: 'copy'

    input:
    path counts
    path samplesheet
    val  contrast
    path r_script

    output:
    path 'deseq2_results.csv',    emit: results
    path 'volcano_plot.pdf',      emit: volcano
    path 'pca_plot.pdf',          emit: pca
    path 'normalized_counts.csv', emit: norm_counts
    path 'versions.yml',          emit: versions

    script:
    def factor      = contrast[0]
    def numerator   = contrast[1]
    def denominator = contrast[2]
    """
    Rscript ${r_script} \
        --counts      ${counts}      \
        --metadata    ${samplesheet} \
        --factor      ${factor}      \
        --numerator   ${numerator}   \
        --denominator ${denominator} \
        --outdir      .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | head -1 | sed 's/R version //' | sed 's/ .*//')
        deseq2: \$(Rscript -e "cat(as.character(packageVersion('DESeq2')))")
    END_VERSIONS
    """
}
