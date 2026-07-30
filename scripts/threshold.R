#!/usr/bin/env Rscript
# scripts/threshold.R
# Preprocesses data and runs WGCNA soft-thresholding parameter scan.

# 1. Parse Command Line Arguments using optparse
library(optparse)

option_list <- list(
  make_option("--counts", type="character", default=NULL, help="Path to counts CSV/TSV file [required]"),
  make_option("--metadata", type="character", default=NULL, help="Path to metadata CSV/TSV file [required]"),
  make_option("--taxonomy", type="character", default=NULL, help="Path to taxonomy CSV/TSV file [optional]"),
  make_option("--sample_id_col", type="character", default="sample_id", help="Column name for sample IDs in metadata"),
  make_option("--outdir", type="character", default="output", help="Output directory path"),
  make_option("--min_num_samples", type="integer", default=10, help="Minimum sample count threshold for filtering"),
  make_option("--transform", type="character", default="hellinger", help="Data transformation type")
)

opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$counts) || is.null(opt$metadata)) {
  stop("Error: --counts and --metadata are required arguments. Use --help for usage details.")
}

counts_file <- opt$counts
metadata_file <- opt$metadata
taxonomy_file <- opt$taxonomy
sample_id_col <- opt$sample_id_col
outdir <- opt$outdir
min_num_samples <- opt$min_num_samples
transform <- opt$transform


# 2. Load Required Packages Quietly
suppressPackageStartupMessages({
  library(WGCNA)
  library(vegan)
})

cat("Loading input datasets...\n")
# Load counts and metadata
counts <- read.csv(counts_file, row.names = 1, check.names = FALSE)
metadata <- read.csv(metadata_file, row.names = 1, check.names = FALSE)

# If sample_id_col is in metadata columns, set it as rownames
if (sample_id_col %in% colnames(metadata)) {
  rownames(metadata) <- metadata[[sample_id_col]]
}

# 3. Preprocess and Filter Counts Table
cat(sprintf("Preprocessing data (prevalence filtering: min_num_samples = %d)...\n", min_num_samples))
# Ensure taxa are rows and samples are columns
otumat <- as.matrix(counts)
otumat[is.na(otumat)] <- 0

# Keep taxa present in more than min_num_samples samples
vals <- rowSums(otumat > 0)
otumat_filtered <- otumat[vals > min_num_samples, , drop = FALSE]

cat(sprintf("Taxa count before filtering: %d, after filtering: %d\n", nrow(otumat), nrow(otumat_filtered)))

# Transpose counts matrix for WGCNA (samples as rows, taxa as columns)
data_t <- t(otumat_filtered)

# 4. Standardisation
if (transform == "hellinger") {
  cat("Applying Hellinger transformation...\n")
  data_transformed <- decostand(data_t, method = "hellinger")
} else if (transform == "relative") {
  cat("Applying relative abundance transformation...\n")
  data_transformed <- decostand(data_t, method = "total")
} else {
  cat("No transformation applied.\n")
  data_transformed <- data_t
}

# 5. Check good samples and genes
cat("Checking good samples and genes using WGCNA...\n")
gsg <- goodSamplesGenes(data_transformed, verbose = 3)
if (!gsg$allOK) {
  if (sum(!gsg$goodGenes) > 0) {
    cat(sprintf("Removing %d low-quality genes.\n", sum(!gsg$goodGenes)))
  }
  if (sum(!gsg$goodSamples) > 0) {
    cat(sprintf("Removing %d low-quality samples.\n", sum(!gsg$goodSamples)))
  }
  data_transformed <- data_transformed[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
}

# Save preprocessed matrix for downstream commands
prep_path <- file.path(outdir, "preprocessed_data.RData")
save(data_transformed, metadata, file = prep_path)
cat(sprintf("Preprocessed data saved to: %s\n", prep_path))

# 6. Sample Clustering plot
cat("Generating sample outlier clustering dendrogram...\n")
sampleTree <- hclust(dist(data_transformed), method = "average")

pdf_path <- file.path(outdir, "thresholding_diagnostics.pdf")
pdf(pdf_path, width = 12, height = 9)

# Layout for plots: Page 1 has Sample Tree
par(cex = 0.8)
par(mar = c(10, 4, 4, 2))
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 1.8)

# Page 2 has threshold parameter plots
allowWGCNAThreads()
powers <- c(1:10, seq(from = 12, to = 20, by = 2))
cat("Scanning soft thresholding powers...\n")
sft <- pickSoftThreshold(
  data_transformed,
  powerVector = powers,
  verbose = 5,
  networkType = "signed"
)

par(mfrow = c(1, 2))
par(mar = c(5, 5, 4, 2))
cex1 <- 0.9

plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n",
     main = "Scale independence"
)
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red"
)
abline(h = 0.80, col = "red")

plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = "Mean connectivity"
)
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers,
     cex = cex1, col = "red"
)

dev.off()
cat(sprintf("Diagnostic plots saved to: %s\n", pdf_path))
cat("Thresholding analysis completed successfully!\n")
