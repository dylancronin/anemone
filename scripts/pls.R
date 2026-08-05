#!/usr/bin/env Rscript
# scripts/pls.R
# Fits Partial Least Squares (PLS) regression for target Module-Trait pairs and calculates VIP scores.

# 1. Parse Command Line Arguments using optparse
library(optparse)

option_list <- list(
  make_option("--counts", type="character", default=NULL, help="Path to counts CSV/TSV file [required]"),
  make_option("--metadata", type="character", default=NULL, help="Path to metadata CSV/TSV file [required]"),
  make_option("--taxonomy", type="character", default=NULL, help="Path to taxonomy CSV/TSV file [optional]"),
  make_option("--sample_id_col", type="character", default="sample_id", help="Column name for sample IDs in metadata"),
  make_option("--outdir", type="character", default="output", help="Output directory path"),
  make_option("--min_prevalence_pct", type="double", default=10, help="Minimum prevalence percentage of samples (0-100)"),
  make_option("--min_num_samples", type="integer", default=NULL, help="Minimum absolute sample count (strict inequality >)"),
  make_option("--transform", type="character", default="hellinger", help="Data transformation type"),
  make_option("--module", type="character", default=NULL, help="Target WGCNA module name [required]"),
  make_option("--parameter", type="character", default=NULL, help="Target metadata trait parameter [required]"),
  make_option("--power", type="integer", default=12, help="Soft-threshold power [ignored]"),
  make_option("--min_r2", type="double", default=0.30, help="Minimum R^2 threshold to proceed to VIP calculation")
)

opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$counts) || is.null(opt$metadata) || is.null(opt$module) || is.null(opt$parameter)) {
  stop("Error: --counts, --metadata, --module, and --parameter are all required arguments.")
}

counts_file <- opt$counts
metadata_file <- opt$metadata
taxonomy_file <- opt$taxonomy
sample_id_col <- opt$sample_id_col
outdir <- file.path(opt$outdir, "pls_vip")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
min_prevalence_pct <- opt$min_prevalence_pct
min_num_samples <- opt$min_num_samples
transform_type <- opt$transform
target_module <- opt$module
target_param <- opt$parameter
min_r2 <- opt$min_r2

# 2. Load Packages Quietly
suppressPackageStartupMessages({
  library(pls)
  library(vegan)
  library(WGCNA)
})

# 3. Load preprocessed data cache or preprocess on the fly
prep_path <- file.path(opt$outdir, "preprocess", "preprocessed_data.RData")

if (file.exists(prep_path)) {
  cat(sprintf("Loading cached preprocessed data from: %s\n", prep_path))
  load(prep_path) # Loads data_transformed, counts, metadata, otumat_filtered
} else {
  cat("No preprocessed data cache found. Preprocessing raw files...\n")
  counts <- read.csv(counts_file, row.names = 1, check.names = FALSE)
  metadata <- read.csv(metadata_file, row.names = 1, check.names = FALSE)
  if (sample_id_col %in% colnames(metadata)) {
    rownames(metadata) <- metadata[[sample_id_col]]
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
    min_samples <- ceiling((min_prevalence_pct / 100) * n_samples)
    vals <- rowSums(otumat > 0)
    otumat_filtered <- otumat[vals >= min_samples, , drop = FALSE]
  }
  
  data_t <- t(otumat_filtered)
  if (transform_type == "hellinger") {
    data_transformed <- decostand(data_t, method = "hellinger")
  } else if (transform_type == "relative") {
    data_transformed <- decostand(data_t, method = "total")
  } else {
    data_transformed <- data_t
  }
}

# 4. Load Module Membership to find features in the target module
membership_file <- file.path(opt$outdir, "network", "module_membership.csv")
if (!file.exists(membership_file)) {
  stop(sprintf("Error: Module membership file not found at %s. Run 'network' step first.", membership_file))
}
membership <- read.csv(membership_file, row.names = 1, check.names = FALSE)

# Filter features in target module
module_features <- rownames(membership)[membership$ModuleColor == target_module]
if (length(module_features) == 0) {
  stop(sprintf("Error: No features found belonging to module '%s'. Check your module name.", target_module))
}
cat(sprintf("Found %d features belonging to module '%s'.\n", length(module_features), target_module))

# Load network construction workspace for module eigengenes (pops)
network_rdata <- file.path(opt$outdir, "network", "network_construction.RData")
if (file.exists(network_rdata)) {
  load(network_rdata)
} else {
  stop(sprintf("Error: Network workspace file not found at %s. Run 'network' step first.", network_rdata))
}

# Subset data_transformed to only target module features
X_data <- data_transformed[, colnames(data_transformed) %in% module_features, drop = FALSE]
if (ncol(X_data) == 0) {
  stop("Error: No overlapping features between preprocessed data and module membership.")
}

# 5. Align with target metadata trait
if (!target_param %in% colnames(metadata)) {
  stop(sprintf("Error: Target parameter '%s' not found in metadata.", target_param))
}

samples_intersect <- intersect(rownames(X_data), rownames(metadata))
if (length(samples_intersect) == 0) {
  stop("Error: No overlapping samples found between data and metadata.")
}

X <- as.matrix(X_data[samples_intersect, , drop = FALSE])
y_raw <- metadata[samples_intersect, target_param]
if (!is.numeric(y_raw)) {
  factor_vals <- as.factor(y_raw)
  levels_map <- levels(factor_vals)
  y <- as.numeric(factor_vals)
  
  map_details <- paste(sprintf("'%s'=%d", levels_map, 1:length(levels_map)), collapse=", ")
  warning(sprintf("Target parameter '%s' is categorical and has been auto-encoded to numeric using levels: [%s]", 
                  target_param, map_details), call. = FALSE)
} else {
  y <- as.numeric(y_raw)
}

# Remove samples with NA in y
valid_samples <- !is.na(y)
if (sum(valid_samples) < 5) {
  stop("Error: Too few samples (< 5) remain with non-NA values for the target parameter.")
}
X <- X[valid_samples, , drop = FALSE]
y <- y[valid_samples]

cat(sprintf("Fitting PLS model predicting '%s' using %d features across %d samples...\n", 
            target_param, ncol(X), nrow(X)))

# 6. Fit PLS Regression with Cross-Validation
max_comps <- min(5, ncol(X), nrow(X) - 2)
if (max_comps < 1) {
  stop("Error: Data dimensions are too small to fit PLS.")
}

pls_model <- plsr(y ~ X, ncomp = max_comps, validation = "LOO", method = "oscorespls")

# R^2 Cutoff Check (matching th_r2 <- 0.3 in original refactored.Rmd)
r2_res <- R2(pls_model)
r2_vals <- as.vector(r2_res$val)
max_r2 <- max(r2_vals[-1], na.rm = TRUE)

if (max_r2 < min_r2) {
  cat(sprintf("Skipping PLS-VIP: max R^2 = %.4f is below the threshold of %.2f. No good correlation.\n", 
              max_r2, min_r2))
  quit(save = "no", status = 0)
}

# Select optimal number of components based on LOO MSEP
cv_res <- RMSEP(pls_model, estimate = "CV")
# Find index of minimum LOO error (excluding 0 components)
best_comp <- which.min(cv_res$val[1, 1, -1])
if (length(best_comp) == 0 || best_comp == 0) {
  best_comp <- 1
}
cat(sprintf("Optimal number of PLS components selected by LOO: %d\n", best_comp))

# Extract final predictions and fit statistics
y_pred <- predict(pls_model, ncomp = best_comp, newdata = X)
r2_val <- cor(y, y_pred)^2
rmse_val <- sqrt(mean((y - y_pred)^2))

# 7. Local VIP Score Calculation Function
# Computes Variable Importance in Projection (VIP) scores for a PLS model.
calculate_vip <- function(pls_obj, ncomp) {
  W <- pls_obj$loading.weights[, 1:ncomp, drop = FALSE]
  q <- pls_obj$Yloadings[, 1:ncomp, drop = FALSE]
  TT <- pls_obj$scores[, 1:ncomp, drop = FALSE]
  p <- nrow(W)
  
  ssy <- numeric(ncomp)
  for (a in 1:ncomp) {
    ssy[a] <- sum(q[, a]^2) * sum(TT[, a]^2)
  }
  
  vip <- numeric(p)
  for (j in 1:p) {
    val <- 0
    for (a in 1:ncomp) {
      w_norm <- W[, a] / sqrt(sum(W[, a]^2))
      val <- val + ssy[a] * (w_norm[j]^2)
    }
    vip[j] <- sqrt(p * val / sum(ssy))
  }
  names(vip) <- rownames(W)
  return(vip)
}

vip_scores <- calculate_vip(pls_model, best_comp)
coef_vals <- coef(pls_model, ncomp = best_comp)[, 1, 1]

# Compute Module Membership (MM) and Trait Significance (GS) statistics for features
nSamples_full <- nrow(data_transformed)
modNames <- substring(colnames(pops), 3)
geneModuleMembership <- as.data.frame(cor(data_transformed, pops, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples_full))
colnames(geneModuleMembership) <- paste0("MM", modNames)
colnames(MMPvalue) <- paste0("p.MM", modNames)

# Trait significance for features against target parameter y
samples_valid_all <- intersect(rownames(data_transformed), rownames(X))
geneTraitSignificance <- as.data.frame(cor(data_transformed[samples_valid_all, , drop = FALSE], y, use = "p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), length(y)))
colnames(geneTraitSignificance) <- paste0("GS.", target_param)
colnames(GSPvalue) <- paste0("p.GS.", target_param)

# 8. Load Taxonomy mapping if available
tax_df <- NULL
if (!is.null(taxonomy_file) && file.exists(taxonomy_file)) {
  cat(sprintf("Loading taxonomy map from: %s\n", taxonomy_file))
  tax_df <- read.csv(taxonomy_file, row.names = 1, check.names = FALSE)
}

# 9. Format Ranking Table
mm_target_col <- paste0("MM", target_module)
pmm_target_col <- paste0("p.MM", target_module)

feats_target <- names(vip_scores)
mm_vals <- geneModuleMembership[feats_target, mm_target_col]
pmm_vals <- MMPvalue[feats_target, pmm_target_col]
gs_vals <- geneTraitSignificance[feats_target, 1]
pgs_vals <- GSPvalue[feats_target, 1]

ranking_df <- data.frame(
  Feature = feats_target,
  VIP = vip_scores,
  Coefficient = coef_vals,
  Correlation = gs_vals,
  ModuleMembership = mm_vals,
  MM_Pvalue = pmm_vals,
  TraitSignificance = gs_vals,
  GS_Pvalue = pgs_vals,
  stringsAsFactors = FALSE
)

# Merge Taxonomy columns if available
if (!is.null(tax_df)) {
  intersect_feats <- intersect(ranking_df$Feature, rownames(tax_df))
  tax_sub <- tax_df[ranking_df$Feature, , drop = FALSE]
  ranking_df <- cbind(ranking_df, tax_sub)
}

# Sort by VIP score descending
ranking_df <- ranking_df[order(-ranking_df$VIP), ]

# Save ranking CSV table
csv_out <- file.path(outdir, sprintf("pls_%s_%s_vip_rankings.csv", target_module, target_param))
write.csv(ranking_df, file = csv_out, row.names = FALSE)
cat(sprintf("Feature rankings saved to: %s\n", csv_out))

# 10. Plot Diagnostic Reports (Page 1: MM vs GS, Page 2: Eigengene Adjacency/Dendrogram, Page 3: PLS LOO & VIP)
plots_pdf <- file.path(outdir, sprintf("pls_%s_%s_plots.pdf", target_module, target_param))
graphics.off()
unlink(plots_pdf)
pdf(plots_pdf, width = 12, height = 8)

# Page 1: Module Membership vs Trait Significance Scatter Plot
if (mm_target_col %in% colnames(geneModuleMembership)) {
  mod_genes_in_data <- intersect(module_features, rownames(geneModuleMembership))
  par(mfrow = c(1, 1))
  par(mar = c(6, 6, 4, 4))
  verboseScatterplot(abs(geneModuleMembership[mod_genes_in_data, mm_target_col]),
                     abs(geneTraitSignificance[mod_genes_in_data, 1]),
                     xlab = paste("Module Membership -", target_module),
                     ylab = paste("Correlation to", target_param),
                     lmFnc = lm,
                     abline = TRUE, abline.color = 1, abline.lty = 1,
                     corFnc = "cor", corOptions = "use = 'p'",
                     main = sprintf("Module Membership vs. Trait Significance\n(Module: %s | Trait: %s)", target_module, target_param),
                     cex.main = 1.4, cex.lab = 1.3, cex.axis = 1.1, pch = 21, bg = target_module, col = "black")
}

# Page 2: Eigengene Dendrogram & Adjacency Heatmap with Trait
y_df <- as.data.frame(y)
rownames(y_df) <- rownames(X)
colnames(y_df) <- target_param

pops_matched <- pops[rownames(y_df), , drop = FALSE]
MET <- orderMEs(cbind(pops_matched, y_df))

par(mfrow = c(1, 2))
par(mar = c(5, 5, 4, 2))
plotEigengeneNetworks(MET, paste("Eigengene dendrogram with", target_param), marDendro = c(0, 4, 2, 0), plotHeatmaps = FALSE)
plotEigengeneNetworks(MET, paste("Eigengene adjacency heatmap with", target_param), marHeatmap = c(3, 4, 2, 2), plotDendrograms = FALSE, xLabelsAngle = 90)

# Page 3: PLS LOO Measured-vs-Predicted & Top 15 VIP Barplot
par(mfrow = c(1, 2))
par(mar = c(5, 5, 4, 2))

# Plot A: Measured vs Predicted Scatter plot
plot(y, y_pred, 
     xlab = sprintf("Measured %s", target_param), 
     ylab = sprintf("Predicted %s", target_param),
     main = sprintf("PLS LOO Performance\n(comp = %d, R^2 = %.3f, RMSE = %.3f)", 
                    best_comp, r2_val, rmse_val),
     pch = 19, col = target_module, cex = 1.2,
     cex.lab = 1.2, cex.axis = 1.1)
abline(0, 1, col = "gray50", lty = "dashed", lwd = 1.5)

# Plot B: Barplot of Features with VIP > 1.0
top_vip <- ranking_df[ranking_df$VIP > 1.0, , drop = FALSE]
if (nrow(top_vip) == 0) {
  # Fallback to top 5 if no features exceed 1.0 to ensure a plot is rendered
  top_vip <- head(ranking_df, 5)
}
bar_labels <- top_vip$Feature
if ("Genus" %in% colnames(top_vip)) {
  bar_labels <- paste0(top_vip$Genus, " (", substring(top_vip$Feature, 1, 8), "...)")
}

# Adjust label font size dynamically based on number of VIP > 1 features
cex_font <- ifelse(nrow(top_vip) > 25, 0.5, ifelse(nrow(top_vip) > 15, 0.65, 0.8))

par(mar = c(5, 12, 4, 2))
barplot(rev(top_vip$VIP), 
        names.arg = rev(bar_labels), 
        horiz = TRUE, 
        las = 1, 
        col = target_module,
        xlab = "VIP Score", 
        main = sprintf("Predictors with VIP > 1 in Module: %s\n(Predicting %s)", target_module, target_param),
        cex.names = cex_font,
        cex.lab = 1.1)
abline(v = 1.0, col = "red", lty = "dashed", lwd = 1.5)

dev.off()
cat(sprintf("Diagnostic plots saved to: %s\n", plots_pdf))
cat("PLS analysis completed successfully!\n")

