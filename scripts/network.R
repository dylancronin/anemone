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
  make_option("--min_prevalence_pct", type="double", default=10, help="Minimum prevalence percentage of samples (0-100)"),
  make_option("--min_num_samples", type="integer", default=NULL, help="Minimum absolute sample count (strict inequality >)"),
  make_option("--transform", type="character", default="hellinger", help="Data transformation type"),
  make_option("--power", type="integer", default=14, help="Soft thresholding power"),
  make_option("--TOMType", type="character", default="signed", help="TOM type (signed/unsigned)"),
  make_option("--networkType", type="character", default="signed", help="Network type (signed/unsigned)"),
  make_option("--mergeCutHeight", type="double", default=0.01, help="Cut height threshold for merging modules"),
  make_option("--maxBlockSize", type="integer", default=1800, help="Maximum block size for network analysis")
)

opt <- parse_args(OptionParser(option_list=option_list))

outdir <- file.path(opt$outdir, "network")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
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
prep_path <- file.path(opt$outdir, "preprocess", "preprocessed_data.RData")

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
  
  # Calculate min_num_samples from prevalence percentage or use absolute threshold
  if (!is.null(opt$min_num_samples)) {
    min_samples <- opt$min_num_samples
    vals <- rowSums(otumat > 0)
    otumat_filtered <- otumat[vals > min_samples, , drop = FALSE]
  } else {
    n_samples <- ncol(otumat)
    min_samples <- ceiling((opt$min_prevalence_pct / 100) * n_samples)
    vals <- rowSums(otumat > 0)
    otumat_filtered <- otumat[vals >= min_samples, , drop = FALSE]
  }
  
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

# 6. Extract Module Eigengenes & Quantitative Module Membership (MM)
cat("Extracting module eigengenes...\n")
MAGs_eig <- moduleEigengenes(data_transformed, module_colors)$eigengenes
pops <- orderMEs(MAGs_eig)
rownames(pops) <- rownames(data_transformed)

# Write module eigengenes table
eigengenes_csv <- file.path(outdir, "module_eigengenes.csv")
write.csv(pops, file = eigengenes_csv, row.names = TRUE)
cat(sprintf("Module eigengenes saved to: %s\n", eigengenes_csv))

# Calculate quantitative Module Membership (MM) and Student p-values
cat("Calculating quantitative module membership (MM)...\n")
nSamples <- nrow(data_transformed)
modNames <- substring(colnames(pops), 3) # Strip ME prefix
geneModuleMembership <- as.data.frame(cor(data_transformed, pops, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
colnames(geneModuleMembership) <- paste0("MM", modNames)
colnames(MMPvalue) <- paste0("p.MM", modNames)

# Write full module membership matrix CSV
gene_mm_csv <- file.path(outdir, "gene_module_membership.csv")
full_mm_df <- cbind(Feature = rownames(geneModuleMembership), geneModuleMembership, MMPvalue)
write.csv(full_mm_df, file = gene_mm_csv, row.names = FALSE)
cat(sprintf("Full gene module membership saved to: %s\n", gene_mm_csv))

# Annotate module_membership.csv with MM_ownModule and p.MM_ownModule
clustermembership$MM_ownModule <- NA_real_
clustermembership$p.MM_ownModule <- NA_real_

for (i in 1:nrow(clustermembership)) {
  feat <- rownames(clustermembership)[i]
  col_color <- clustermembership$ModuleColor[i]
  mm_col <- paste0("MM", col_color)
  pmm_col <- paste0("p.MM", col_color)
  if (mm_col %in% colnames(geneModuleMembership)) {
    clustermembership$MM_ownModule[i] <- geneModuleMembership[feat, mm_col]
    clustermembership$p.MM_ownModule[i] <- MMPvalue[feat, pmm_col]
  }
}

membership_csv <- file.path(outdir, "module_membership.csv")
write.csv(clustermembership, file = membership_csv, row.names = TRUE)
cat(sprintf("Module membership saved to: %s\n", membership_csv))

# Plot the unmerged module eigengene clustering dendrogram to help choose mergeCutHeight
cat("Generating unmerged module eigengene clustering dendrogram...\n")
unmerged_colors <- labels2colors(net$unmergedColors)
unmerged_eig <- moduleEigengenes(data_transformed, unmerged_colors)$eigengenes
unmerged_pops <- orderMEs(unmerged_eig)
rownames(unmerged_pops) <- rownames(data_transformed)

unmerged_pdf <- file.path(outdir, "unmerged_module_eigengene_clustering.pdf")
graphics.off() # Close any dangling graphics devices
unlink(unmerged_pdf) # Force delete existing PDF if present
pdf(unmerged_pdf, width = 8, height = 6)
unmerged_dissim <- 1 - cor(unmerged_pops)
unmerged_meTree <- hclust(as.dist(unmerged_dissim), method = "average")

par(mar = c(5, 5, 4, 2))
plot(unmerged_meTree, main = "Unmerged Eigengene dendrogram (Before Merging)", xlab = "", sub = "", cex = 0.8)
abline(h = seq(0, 2, by = 0.05), col = "gray80", lty = "dashed", lwd = 0.5)
dev.off()

# Plot the merged module eigengene clustering dendrogram to verify final state
cat("Generating merged module eigengene clustering dendrogram...\n")
merged_pdf <- file.path(outdir, "module_eigengene_clustering.pdf")
unlink(merged_pdf) # Force delete existing PDF if present
pdf(merged_pdf, width = 8, height = 6)

# Calculate dissimilarity and cluster on merged eigengenes
dissim <- 1 - cor(pops)
meTree <- hclust(as.dist(dissim), method = "average")

# Plot with native y-axis to ensure perfect alignment
par(mar = c(5, 5, 4, 2))
plot(meTree, main = "Merged Eigengene dendrogram (Final State)", xlab = "", sub = "", cex = 0.8)

# Add horizontal dashed gridlines for easy reading
abline(h = seq(0, 2, by = 0.05), col = "gray80", lty = "dashed", lwd = 0.5)

dev.off()

# 7. Save WGCNA outputs to RData for downstream scripts (correlate, pls)
network_rdata <- file.path(outdir, "network_construction.RData")
geneTree <- net$dendrograms[[1]]
moduleColors <- module_colors
moduleLabels <- net$colors

save(pops, moduleLabels, moduleColors, geneTree, geneModuleMembership, MMPvalue, file = network_rdata)
cat(sprintf("Network construction workspace saved to: %s\n", network_rdata))
cat("Network construction completed successfully!\n")

