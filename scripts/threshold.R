#!/usr/bin/env Rscript
# scripts/threshold.R
# Preprocesses data and runs WGCNA soft-thresholding parameter scans for different prevalence thresholds.

# 1. Parse Command Line Arguments using optparse
library(optparse)

option_list <- list(
  make_option("--counts", type="character", default=NULL, help="Path to counts CSV/TSV file [required]"),
  make_option("--metadata", type="character", default=NULL, help="Path to metadata CSV/TSV file [required]"),
  make_option("--taxonomy", type="character", default=NULL, help="Path to taxonomy CSV/TSV file [optional]"),
  make_option("--sample_id_col", type="character", default="sample_id", help="Column name for sample IDs in metadata"),
  make_option("--outdir", type="character", default="output", help="Output directory path"),
  make_option("--min_prevalence_pct", type="double", default=10, help="Minimum prevalence percentage of samples (0-100)"),
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
outdir <- file.path(opt$outdir, "preprocess")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
min_prevalence_pct <- opt$min_prevalence_pct
transform <- opt$transform

# 2. Load Required Packages Quietly
suppressPackageStartupMessages({
  library(WGCNA)
  library(vegan)
})

cat("Loading input datasets...\n")
counts <- read.csv(counts_file, row.names = 1, check.names = FALSE)
metadata <- read.csv(metadata_file, row.names = 1, check.names = FALSE)

if (sample_id_col %in% colnames(metadata)) {
  rownames(metadata) <- metadata[[sample_id_col]]
}

# Helper function to preprocess data for a given prevalence percentage
preprocess_data <- function(counts_df, metadata_df, sample_id, prev_pct, trans_method) {
  otumat <- as.matrix(counts_df)
  otumat[is.na(otumat)] <- 0
  
  # Calculate min_num_samples from prevalence percentage
  n_samples <- ncol(otumat)
  min_samples <- ceiling((prev_pct / 100) * n_samples)
  
  # Filter by prevalence
  vals <- rowSums(otumat > 0)
  otumat_filtered <- otumat[vals >= min_samples, , drop = FALSE]
  
  # Transpose
  data_t <- t(otumat_filtered)
  
  # Transformation
  if (trans_method == "hellinger") {
    data_transformed <- decostand(data_t, method = "hellinger")
  } else if (trans_method == "relative") {
    data_transformed <- decostand(data_t, method = "total")
  } else {
    data_transformed <- data_t
  }
  
  # goodSamplesGenes check
  gsg <- goodSamplesGenes(data_transformed, verbose = 0)
  if (!gsg$allOK) {
    data_transformed <- data_transformed[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
  }
  
  return(list(data = data_transformed, min_samples = min_samples, n_features = nrow(otumat_filtered)))
}

# 3. Preprocess for the main configured prevalence cutoff
cat(sprintf("Preprocessing configured cutoff: %.1f%% prevalence...\n", min_prevalence_pct))
prep_res <- preprocess_data(counts, metadata, sample_id_col, min_prevalence_pct, transform)
data_transformed <- prep_res$data

cat(sprintf("Configured cutoff (%.1f%%) corresponds to presence in >= %d samples.\n", 
            min_prevalence_pct, prep_res$min_samples))
cat(sprintf("Taxa count before filtering: %d, after filtering: %d\n", 
            nrow(counts), prep_res$n_features))

# Save this main preprocessed workspace for network construction
prep_path <- file.path(outdir, "preprocessed_data.RData")
save(data_transformed, metadata, file = prep_path)
cat(sprintf("Preprocessed data saved to: %s\n", prep_path))

# 4. Generate Diagnostic PDF Plots
pdf_path <- file.path(outdir, "thresholding_diagnostics.pdf")
graphics.off()
unlink(pdf_path)
pdf(pdf_path, width = 12, height = 9)

# Subsequent Pages: Scan a range of prevalence levels to test scale independence and sample clustering
prevalence_grid <- unique(sort(c(5, 10, 15, 20, 25, 30, min_prevalence_pct)))

allowWGCNAThreads()
powers <- c(1:10, seq(from = 12, to = 20, by = 2))

for (pct in prevalence_grid) {
  cat(sprintf("\n==================================================\n"))
  cat(sprintf("Scanning Soft Thresholds for Prevalence: %.1f%%\n", pct))
  cat(sprintf("==================================================\n"))
  
  # Preprocess for this grid point
  grid_prep <- preprocess_data(counts, metadata, sample_id_col, pct, transform)
  
  if (grid_prep$n_features < 10) {
    cat(sprintf("Skipping %.1f%%: Too few features remaining (%d)\n", pct, grid_prep$n_features))
    next
  }
  
  cat(sprintf("Prevalence %.1f%% -> >= %d samples, keeping %d features.\n", 
              pct, grid_prep$min_samples, grid_prep$n_features))
  
  # A. Plot Sample Outlier Clustering for this prevalence cutoff
  cat(sprintf("Plotting sample clustering for %.1f%% prevalence...\n", pct))
  grid_sampleTree <- hclust(dist(grid_prep$data), method = "average")
  
  par(mfrow = c(1, 1))
  par(cex = 0.8)
  par(mar = c(10, 4, 4, 2))
  plot(grid_sampleTree, 
       main = sprintf("Sample clustering (Prevalence cutoff: %.1f%%, features: %d)", pct, grid_prep$n_features), 
       sub = "", xlab = "", cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.4)
  
  # B. Run WGCNA soft-threshold parameter scan
  sft <- pickSoftThreshold(
    grid_prep$data,
    powerVector = powers,
    verbose = 5,
    networkType = "signed"
  )
  
  # C. Plot scale independence & mean connectivity results
  par(mfrow = c(1, 2))
  par(mar = c(5, 5, 4, 2))
  cex1 <- 0.9
  
  plot(sft$fitIndices[, 1],
       -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
       xlab = "Soft Threshold (power)",
       ylab = "Scale Free Topology Model Fit, signed R^2",
       type = "n",
       main = sprintf("Scale independence (Prevalence: %.1f%%, features: %d)", pct, grid_prep$n_features),
       cex.main = 1.1
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
       main = sprintf("Mean connectivity (Prevalence: %.1f%%)", pct),
       cex.main = 1.1
  )
  text(sft$fitIndices[, 1],
       sft$fitIndices[, 5],
       labels = powers,
       cex = cex1, col = "red"
  )
}

dev.off()
cat(sprintf("\nDiagnostic plots saved to: %s\n", pdf_path))
cat("Thresholding analysis completed successfully!\n")

