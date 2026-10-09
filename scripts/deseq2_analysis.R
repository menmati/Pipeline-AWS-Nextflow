#!/usr/bin/env Rscript
# =============================================================================
# DESeq2 Differential Expression Analysis
# Project  : PRJNA1533324 – miR-150-5p inhibition in mouse embryonic cortex
# Contrast : miR150_inhibition vs. negative_control
# Outputs  : deseq2_results.csv, volcano_plot.pdf, pca_plot.pdf,
#             normalized_counts.csv
# =============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
  library(optparse)
})

# ── Command-line argument parsing ─────────────────────────────
option_list <- list(
  make_option("--counts",      type="character", help="featureCounts output (.txt)"),
  make_option("--metadata",    type="character", help="Samplesheet CSV (sample, condition ...)"),
  make_option("--factor",      type="character", default="condition",
              help="Column in metadata used for grouping [default: condition]"),
  make_option("--numerator",   type="character", help="Numerator level of contrast"),
  make_option("--denominator", type="character", help="Denominator (reference) level"),
  make_option("--outdir",      type="character", default=".", help="Output directory"),
  make_option("--alpha",       type="double",    default=0.05,
              help="Adjusted p-value threshold [default: 0.05]"),
  make_option("--lfc_thresh",  type="double",    default=1.0,
              help="Log2FC threshold for volcano labelling [default: 1.0]")
)
opt <- parse_args(OptionParser(option_list=option_list))
dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)

message("=== DESeq2 Analysis ===")
message(sprintf("Contrast : %s vs %s", opt$numerator, opt$denominator))
message(sprintf("Alpha    : %.3f  |  LFC threshold: %.2f", opt$alpha, opt$lfc_thresh))

# ── 1. Load count matrix ──────────────────────────────────────
message("\n[1/6] Loading count matrix: ", opt$counts)
raw <- read.table(opt$counts, header=TRUE, sep="\t", comment.char="#")

bam_cols  <- grep("\\.bam$", colnames(raw), value=TRUE)
counts_mx <- as.matrix(raw[, bam_cols])
rownames(counts_mx) <- raw$Geneid

# Clean BAM path → sample ID (handles both with and without directory paths)
colnames(counts_mx) <- sub("\\.bam$", "", basename(colnames(counts_mx)))
message(sprintf("  Genes: %d  |  Samples: %d (%s)", nrow(counts_mx), ncol(counts_mx), paste(colnames(counts_mx), collapse=", ")))

# ── 2. Load metadata ──────────────────────────────────────────
message("[2/6] Loading metadata: ", opt$metadata)
meta <- read.csv(opt$metadata, stringsAsFactors=FALSE)
rownames(meta) <- meta$sample

shared    <- intersect(colnames(counts_mx), rownames(meta))
counts_mx <- counts_mx[, shared, drop=FALSE]
meta      <- meta[shared, , drop=FALSE]
message("  Samples matched: ", paste(shared, collapse=", "))

# ── 3. Build DESeqDataSet ─────────────────────────────────────
message("[3/6] Building DESeqDataSet")
meta[[opt$factor]] <- factor(meta[[opt$factor]],
                             levels=c(opt$denominator, opt$numerator))

dds <- DESeqDataSetFromMatrix(
  countData = counts_mx,
  colData   = meta,
  design    = as.formula(paste("~", opt$factor))
)

# Pre-filter: keep genes with >= 10 counts total
keep <- rowSums(counts(dds)) >= 10
dds  <- dds[keep, ]
message(sprintf("  Genes after pre-filtering (>=10 counts): %d", nrow(dds)))

# ── 4. Run DESeq2 ─────────────────────────────────────────────
message("[4/6] Running DESeq2...")
dds        <- DESeq(dds)
res_coef   <- resultsNames(dds)
message("  Result names: ", paste(res_coef, collapse=", "))

res <- results(dds,
               contrast = c(opt$factor, opt$numerator, opt$denominator),
               alpha    = opt$alpha)

# LFC shrinkage with apeglm
coef_name <- res_coef[grep(opt$numerator, res_coef)]
res_shrunk <- lfcShrink(dds, coef=coef_name, type="apeglm", quiet=TRUE)

# Annotate significance
res_df <- as.data.frame(res_shrunk)
res_df$gene_id     <- rownames(res_df)
res_df$significant <- with(res_df,
  !is.na(padj) & padj < opt$alpha & abs(log2FoldChange) >= opt$lfc_thresh)
res_df$direction <- with(res_df, ifelse(
  significant & log2FoldChange >  0, "Up",
  ifelse(significant & log2FoldChange < 0, "Down", "NS")))

up_n   <- sum(res_df$direction == "Up",   na.rm=TRUE)
down_n <- sum(res_df$direction == "Down", na.rm=TRUE)
message(sprintf("  Significant DEGs: %d  (Up: %d | Down: %d)",
                up_n + down_n, up_n, down_n))

out_csv <- file.path(opt$outdir, "deseq2_results.csv")
write.csv(res_df[order(res_df$padj, na.last=TRUE), ], out_csv, row.names=FALSE)
message("  Results written to: ", out_csv)

# ── 5. Normalized counts ──────────────────────────────────────
message("[5/6] Saving normalized counts")
norm_df <- as.data.frame(counts(dds, normalized=TRUE))
norm_df$gene_id <- rownames(norm_df)
write.csv(norm_df, file.path(opt$outdir, "normalized_counts.csv"), row.names=FALSE)

# ── 6. PCA plot ───────────────────────────────────────────────
message("[6/6] Generating plots")
vsd      <- vst(dds, blind=TRUE)
pca_data <- plotPCA(vsd, intgroup=opt$factor, returnData=TRUE)
pct_var  <- round(100 * attr(pca_data, "percentVar"))

pal <- c("miR150_inhibition"="#E63946", "negative_control"="#457B9D")

pca_p <- ggplot(pca_data,
                aes(x=PC1, y=PC2,
                    colour=.data[[opt$factor]],
                    label=name)) +
  geom_point(size=5, alpha=0.9) +
  geom_text_repel(size=3.5, max.overlaps=20) +
  scale_colour_manual(values=pal) +
  labs(
    title    = "PCA of VST-normalised counts",
    subtitle = sprintf("Contrast: %s vs %s", opt$numerator, opt$denominator),
    x        = sprintf("PC1: %d%% variance", pct_var[1]),
    y        = sprintf("PC2: %d%% variance", pct_var[2]),
    colour   = opt$factor
  ) +
  theme_bw(base_size=13) +
  theme(legend.position="bottom")

ggsave(file.path(opt$outdir, "pca_plot.pdf"), pca_p,
       width=7, height=6, device="pdf")

# ── Volcano plot ──────────────────────────────────────────────
top_genes <- head(res_df[res_df$significant, ][order(res_df$padj[res_df$significant]), ], 15)

vol_p <- ggplot(res_df,
                aes(x=log2FoldChange, y=-log10(pvalue), colour=direction)) +
  geom_point(alpha=0.5, size=1.3) +
  geom_point(data=top_genes, size=2.8, alpha=0.9) +
  geom_text_repel(data=top_genes,
                  aes(label=gene_id), size=3, max.overlaps=25,
                  box.padding=0.4) +
  geom_vline(xintercept=c(-opt$lfc_thresh, opt$lfc_thresh),
             linetype="dashed", colour="grey50", linewidth=0.4) +
  geom_hline(yintercept=-log10(opt$alpha),
             linetype="dashed", colour="grey50", linewidth=0.4) +
  scale_colour_manual(values=c("Up"="#E63946", "Down"="#457B9D", "NS"="grey75")) +
  labs(
    title    = "Volcano plot – Differential Expression",
    subtitle = sprintf("%s vs %s  |  padj < %.2f, |log2FC| >= %.1f  |  DEGs: %d",
                       opt$numerator, opt$denominator,
                       opt$alpha, opt$lfc_thresh, up_n + down_n),
    x        = expression(log[2]~"fold change"),
    y        = expression(-log[10]~"p-value"),
    colour   = "Direction"
  ) +
  theme_bw(base_size=13) +
  theme(legend.position="bottom")

ggsave(file.path(opt$outdir, "volcano_plot.pdf"), vol_p,
       width=9, height=7, device="pdf")

# ── Summary ───────────────────────────────────────────────────
message("\n=== Summary ===")
message(sprintf("Contrast         : %s vs %s", opt$numerator, opt$denominator))
message(sprintf("Total genes      : %d", nrow(res_df)))
message(sprintf("Significant DEGs : %d  (Up: %d | Down: %d)",
                up_n + down_n, up_n, down_n))
message(sprintf("Outputs in       : %s", opt$outdir))
message("Done.")
