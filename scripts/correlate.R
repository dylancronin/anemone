#!/usr/bin/env Rscript
# scripts/correlate.R
# Computes correlation between WGCNA module eigengenes and metadata traits.

# 1. Parse Command Line Arguments using optparse
library(optparse)

option_list <- list(
  make_option("--metadata", type="character", default=NULL, help="Path to metadata CSV/TSV file [required]"),
  make_option("--sample_id_col", type="character", default="sample_id", help="Column name for sample IDs in metadata"),
  make_option("--outdir", type="character", default="output", help="Output directory path"),
  make_option("--traits", type="character", default="", help="Comma-separated list of metadata traits to correlate")
)

opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$metadata)) {
  stop("Error: --metadata is a required argument.")
}

metadata_file <- opt$metadata
sample_id_col <- opt$sample_id_col
outdir <- file.path(opt$outdir, "correlation")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
traits_str <- opt$traits

# 2. Load Packages Quietly
suppressPackageStartupMessages({
  library(WGCNA)
  library(RColorBrewer)
})

# 3. Load Network Construction Results
network_rdata <- file.path(opt$outdir, "network", "network_construction.RData")
if (!file.exists(network_rdata)) {
  stop(sprintf("Error: Network construction data not found at %s. Run 'network' step first.", network_rdata))
}
cat(sprintf("Loading network workspace from: %s\n", network_rdata))
load(network_rdata) # Loads 'pops' (eigengenes)

# 4. Load and Process Metadata
cat("Loading and matching metadata...\n")
metadata <- read.csv(metadata_file, row.names = 1, check.names = FALSE)
if (sample_id_col %in% colnames(metadata)) {
  rownames(metadata) <- metadata[[sample_id_col]]
}

# Keep only samples present in module eigengenes (pops)
samples_intersect <- intersect(rownames(pops), rownames(metadata))
if (length(samples_intersect) == 0) {
  stop("Error: No overlapping samples found between module eigengenes and metadata.")
}

metadata_sub <- metadata[samples_intersect, , drop = FALSE]
pops_sub <- pops[samples_intersect, , drop = FALSE]

# 5. Extract target traits
if (traits_str != "") {
  target_traits <- unlist(strsplit(traits_str, ","))
  missing_traits <- setdiff(target_traits, colnames(metadata_sub))
  if (length(missing_traits) > 0) {
    warning(sprintf("The following traits are missing from metadata and will be skipped: %s", 
                    paste(missing_traits, collapse=", ")))
  }
  target_traits <- intersect(target_traits, colnames(metadata_sub))
} else {
  # Default to all numeric columns in metadata
  cat("No traits specified. Selecting all numeric columns in metadata...\n")
  numeric_cols <- sapply(metadata_sub, is.numeric)
  target_traits <- names(numeric_cols)[numeric_cols]
}

if (length(target_traits) == 0) {
  stop("Error: No numeric metadata traits available for correlation analysis.")
}

metadata_traits <- metadata_sub[, target_traits, drop = FALSE]

# Convert factor/character columns to numeric if needed
for (col in colnames(metadata_traits)) {
  if (!is.numeric(metadata_traits[[col]])) {
    factor_vals <- as.factor(metadata_traits[[col]])
    levels_map <- levels(factor_vals)
    metadata_traits[[col]] <- as.numeric(factor_vals)
    
    map_details <- paste(sprintf("'%s'=%d", levels_map, 1:length(levels_map)), collapse=", ")
    warning(sprintf("Metadata column '%s' is categorical (factor/character) and has been auto-encoded to numeric using levels: [%s]", 
                    col, map_details), call. = FALSE)
  }
}

# Remove columns with zero variance
var_traits <- apply(metadata_traits, 2, var, na.rm = TRUE)
zero_var_cols <- names(var_traits)[is.na(var_traits) | var_traits == 0]
if (length(zero_var_cols) > 0) {
  cat(sprintf("Removing zero-variance trait columns: %s\n", paste(zero_var_cols, collapse=", ")))
  metadata_traits <- metadata_traits[, !colnames(metadata_traits) %in% zero_var_cols, drop = FALSE]
}

if (ncol(metadata_traits) == 0) {
  stop("Error: No valid numeric metadata traits remaining after variance filtering.")
}

# 6. Compute Pearson Correlation and Student p-values
cat("Calculating module-trait correlations...\n")
nSamples <- nrow(pops_sub)
moduleTraitCor <- cor(pops_sub, metadata_traits, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)


# 8. Write Long-form Output Table with FDR BH correction
p_vector <- as.vector(moduleTraitPvalue)
p_adj_values <- p.adjust(p_vector, method = "BH", n = sum(!is.na(p_vector)))
modules_rep <- rownames(moduleTraitCor)
traits_rep <- colnames(moduleTraitCor)

long_form <- data.frame(
  Module = rep(modules_rep, times = length(traits_rep)),
  Metadata_Parameter = rep(traits_rep, each = length(modules_rep)),
  Correlation = as.vector(moduleTraitCor),
  P_value = as.vector(moduleTraitPvalue),
  p.adj = p_adj_values,
  stringsAsFactors = FALSE
)
long_csv <- file.path(outdir, "longform_module_trait_table.csv")
write.csv(long_form, file = long_csv, row.names = FALSE)
cat(sprintf("Long-form table with BH adjustment saved to: %s\n", long_csv))

# 9. Plot Eigengene-Trait Relationships Heatmap
cat("Generating heatmap plot...\n")
# Calculate FDR adjusted p-values (q-values) matrix for display
p_adj_vector <- p.adjust(as.vector(moduleTraitPvalue), method = "BH")
moduleTraitPadj <- matrix(p_adj_vector, nrow = nrow(moduleTraitPvalue), ncol = ncol(moduleTraitPvalue))

textMatrix <- paste0(
  signif(moduleTraitCor, 2),
  "\np = ", signif(moduleTraitPvalue, 1),
  "\nq = ", signif(moduleTraitPadj, 1)
)
dim(textMatrix) <- dim(moduleTraitCor)

pdf_path <- file.path(outdir, "module_trait_relationships.pdf")
# Dynamically scale height and width based on number of traits and modules
plot_width <- max(10, ncol(moduleTraitCor) * 0.5 + 4)
plot_height <- max(8, nrow(moduleTraitCor) * 0.4 + 4)

graphics.off()
unlink(pdf_path)
pdf(pdf_path, width = plot_width, height = plot_height)
par(mar = c(10, 14, 4, 5))

labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = colnames(metadata_traits),
  yLabels = colnames(pops_sub),
  ySymbols = substring(colnames(pops_sub), 3), # Strips the "ME" prefix from MEblue -> blue
  colorLabels = FALSE,
  colors = brewer.pal(11, "RdBu"),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 0.6,
  cex.lab = 1.0,
  zlim = c(-1, 1),
  main = "Module-trait relationships"
)
dev.off()
cat(sprintf("Heatmap plot saved to: %s\n", pdf_path))
cat("Correlation analysis completed successfully!\n")
