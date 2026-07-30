#!/usr/bin/env Rscript
# scripts/network.R
# Loads preprocessed data and runs WGCNA network construction.

# 1. Parse Command Line Arguments using optparse
library(optparse)

option_list <- list(
  make_option("--counts", type="character", default=NULL, help="Path to counts CSV/TSV file"),
  make_option("--metadata", type="character", default=NULL, help="Path to metadata CSV/TSV file"),
  make_option("--taxonomy", type="character", default=NULL, help="Path to taxonomy CSV/TSV file [optional]"),
  make_option("--sample_id_col", type="character", default="sample_id", help="Column name for sample IDs in metadata"),
  make_option("--outdir", type="character", default="output", help="Output directory path"),
  make_option("--min_num_samples", type="integer", default=10, help="Minimum sample count threshold for filtering"),
  make_option("--transform", type="character", default="hellinger", help="Data transformation type"),
  make_option("--power", type="integer", default=14, help="Soft thresholding power"),
  make_option("--TOMType", type="character", default="signed", help="TOM type (signed/unsigned)"),
  make_option("--networkType", type="character", default="signed", help="Network type (signed/unsigned)"),
  make_option("--mergeCutHeight", type="double", default=0.01, help="Cut height threshold for merging modules"),
  make_option("--maxBlockSize", type="integer", default=1800, help="Maximum block size for network analysis")
)

opt <- parse_args(OptionParser(option_list=option_list))

outdir <- opt$outdir
power <- opt$power
tom_type <- opt$TOMType
net_type <- opt$networkType
merge_cut <- opt$mergeCutHeight
max_block <- opt$maxBlockSize

# 2. Load Packages Quietly
suppressPackageStartupMessages({
  library(WGCNA)
  library(vegan)
})

# 3. Load preprocessed data cache or preprocess on the fly
prep_path <- file.path(outdir, "preprocessed_data.RData")

if (file.exists(prep_path)) {
  cat(sprintf("Loading cached preprocessed data from: %s\n", prep_path))
  load(prep_path)
} else {
  cat("No preprocessed data cache found. Preprocessing raw files...\n")
  if (is.null(opt$counts) || is.null(opt$metadata)) {
    stop("Error: counts and metadata are required since no preprocessed cache was found.")
  }
  
  counts <- read.csv(opt$counts, row.names = 1, check.names = FALSE)
  metadata <- read.csv(opt$metadata, row.names = 1, check.names = FALSE)
  if (opt$sample_id_col %in% colnames(metadata)) {
    rownames(metadata) <- metadata[[opt$sample_id_col]]
  }
  
  otumat <- as.matrix(counts)
  otumat[is.na(otumat)] <- 0
  vals <- rowSums(otumat > 0)
  otumat_filtered <- otumat[vals > opt$min_num_samples, , drop = FALSE]
  
  data_t <- t(otumat_filtered)
  if (opt$transform == "hellinger") {
    data_transformed <- decostand(data_t, method = "hellinger")
  } else if (opt$transform == "relative") {
    data_transformed <- decostand(data_t, method = "total")
  } else {
    data_transformed <- data_t
  }
  
  gsg <- goodSamplesGenes(data_transformed, verbose = 3)
  if (!gsg$allOK) {
    data_transformed <- data_transformed[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
  }
}

# 4. Construct Network Modules
cat(sprintf("Running WGCNA blockwiseModules (power = %d, networkType = %s, TOMType = %s)...\n", 
            power, net_type, tom_type))
allowWGCNAThreads()

net <- blockwiseModules(
  data_transformed,
  power = power,
  TOMType = tom_type,
  networkType = net_type,
  reassignThreshold = 0,
  mergeCutHeight = merge_cut,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = TRUE,
  saveTOMFileBase = file.path(outdir, "network_TOM"),
  verbose = 3,
  maxBlockSize = max_block
)

# 5. Extract and format module assignments
cat("Formatting module membership output...\n")
module_colors <- labels2colors(net$colors)
clusterid <- as.data.frame(net$colors)
clusteridcolors <- as.data.frame(module_colors)
rownames(clusteridcolors) <- rownames(clusterid)

clustermembership <- merge(clusterid, clusteridcolors, by = 'row.names')
rownames(clustermembership) <- clustermembership$Row.names
clustermembership$Row.names <- NULL
colnames(clustermembership) <- c("ModuleLabel", "ModuleColor")

# Write module membership list
membership_csv <- file.path(outdir, "module_membership.csv")
write.csv(clustermembership, file = membership_csv, row.names = TRUE)
cat(sprintf("Module membership saved to: %s\n", membership_csv))

# 6. Extract Module Eigengenes
cat("Extracting module eigengenes...\n")
MAGs_eig <- moduleEigengenes(data_transformed, module_colors)$eigengenes
pops <- orderMEs(MAGs_eig)
rownames(pops) <- rownames(data_transformed)

# Write module eigengenes table
eigengenes_csv <- file.path(outdir, "module_eigengenes.csv")
write.csv(pops, file = eigengenes_csv, row.names = TRUE)
cat(sprintf("Module eigengenes saved to: %s\n", eigengenes_csv))

# 7. Save WGCNA outputs to RData for downstream scripts (correlate, pls)
network_rdata <- file.path(outdir, "network_construction.RData")
geneTree <- net$dendrograms[[1]]
moduleColors <- module_colors
moduleLabels <- net$colors

save(pops, moduleLabels, moduleColors, geneTree, file = network_rdata)
cat(sprintf("Network construction workspace saved to: %s\n", network_rdata))
cat("Network construction completed successfully!\n")
