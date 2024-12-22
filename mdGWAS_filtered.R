# Load required libraries
library(data.table)
library(ggplot2)
library(gridExtra)
library(DCA)
library(impute)

# Set paths
input_dir <- "data/input/"
output_dir <- "data/output/"
results_dir <- "results/"

# ======================
# MODULE 1: Preprocess Metabolites
# ======================
preprocess_metabolites <- function(input_path, output_path) {
  # Load raw data
  metabolites <- fread(input_path)
  
  # Impute missing values
  imputation <- impute.knn(t(metabolites[,-1]), k = 5)
  metabolites_imputed <- data.frame(
    studyid = metabolites[, 1],
    t(imputation$data)
  )
  
  # Save preprocessed data
  write.csv(metabolites_imputed, output_path, row.names = FALSE)
}

# ======================
# MODULE 2: Run DCA
# ======================
run_dca <- function(data, n_factors, output_file) {
  # Run DCA analysis
  dca_results <- DCA::dca_witho(
    t(data),
    n.fac = n_factors,
    normalization = "standardize"
  )
  
  # Save results
  saveRDS(dca_results, file = output_file)
}

# ======================
# MODULE 3: DCA Visualization
# ======================
plot_variance_explained <- function(dca_results, output_file) {
  variance_df <- data.frame(
    Factor = 1:length(dca_results$ss.proj),
    Variance = dca_results$ss.proj * 100
  )
  
  plot <- ggplot(variance_df, aes(x = Factor, y = Variance)) +
    geom_bar(stat = "identity", fill = "#69b3a2") +
    labs(title = "Variance Explained by Factors", x = "Factors", y = "Variance (%)") +
    theme_minimal()
  
  ggsave(output_file, plot)
}

# ======================
# MODULE 4: Factor Selection
# ======================
select_factors <- function(dca_results, top_n, output_file) {
  # Select top factors based on explained variance
  factor_importance <- data.frame(
    Factor = 1:length(dca_results$ss.proj),
    Variance = dca_results$ss.proj * 100
  )
  top_factors <- factor_importance[order(-factor_importance$Variance), ][1:top_n, ]
  
  # Save selected factors
  write.csv(top_factors, output_file, row.names = FALSE)
}

# ======================
# MODULE 5: Compare DCA Results (10 vs 20 Factors)
# ======================
compare_dca_results <- function(dca_10, dca_20, output_file) {
  # Compare explained variance
  comparison <- data.frame(
    Factor = 1:20,
    Variance_10 = c(dca_10$ss.proj * 100, rep(NA, 10)),
    Variance_20 = dca_20$ss.proj * 100
  )
  
  # Save comparison table
  write.csv(comparison, output_file, row.names = FALSE)
}

# ======================
# MODULE 6: Correlation Analysis
# ======================
metabolite_correlation <- function(data, output_file) {
  # Compute correlation matrix
  cor_matrix <- cor(data)
  
  # Save correlation matrix
  write.csv(cor_matrix, output_file, row.names = TRUE)
}

plot_correlation_heatmap <- function(cor_matrix, output_file) {
  heatmap_df <- as.data.frame(as.table(cor_matrix))
  
  plot <- ggplot(data = heatmap_df, aes(x = Var1, y = Var2, fill = Freq)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(
      low = "blue", high = "red", mid = "white",
      midpoint = 0, limit = c(-1, 1), name = "Correlation"
    ) +
    theme_minimal() +
    labs(title = "Metabolite Correlation Heatmap", x = "Metabolites", y = "Metabolites") +
    theme(axis.text.x = element_blank(), axis.text.y = element_blank())
  
  ggsave(output_file, plot)
}

# ======================
# MODULE 7: Conditional Correlation Analysis
# ======================
conditional_correlation <- function(metabolites, genotypes, output_file) {
  # Compute conditional correlation
  conditional_corr <- cor(metabolites, genotypes)
  
  # Save results
  write.csv(conditional_corr, output_file, row.names = TRUE)
}

# ======================
# MODULE 8: Compare Related vs Unrelated Samples
# ======================
compare_relateds_unrelateds <- function(factors_related, factors_unrelated, output_file) {
  # Compare correlations
  correlations <- cor(factors_related, factors_unrelated)
  
  # Save results
  write.csv(correlations, output_file, row.names = TRUE)
}

# ======================
# MODULE 9: Compare Single vs Double Normalization
# ======================
compare_single_double <- function(data, output_file) {
  double_inverse_normalize <- function(data) {
    data <- apply(data, 2, function(x) {
      residual <- residuals(lm(x ~ 1))
      return(qnorm((rank(residual) - 0.5) / length(residual)))
    })
    return(data)
  }
  
  # Double normalization
  data_double <- double_inverse_normalize(data)
  
  # Save results
  write.csv(data_double, output_file, row.names = TRUE)
}

# ======================
# MODULE 10: EPACTS Analysis
# ======================
run_epacts_pipeline <- function(trait, chr, vcf_dir, out_dir, pheno_file) {
  # Construct EPACTS command
  epacts_cmd <- paste(
    "epacts single",
    "--vcf", paste0(vcf_dir, "/chr", chr, ".vcf.gz"),
    "--pheno", pheno_file,
    "--test q.emmax",
    "--out", paste0(out_dir, "/", trait, "_chr", chr),
    sep = " "
  )
  
  # Run the command
  system(epacts_cmd)
  
  # Run the makefile for EPACTS
  system(paste("make -f", paste0(out_dir, "/", trait, "_chr", chr, ".Makefile")))
}

analyze_epacts <- function(epacts_results_file, output_file) {
  # Load EPACTS results
  epacts_results <- fread(epacts_results_file)
  
  # Extract significant SNPs
  significant_snps <- epacts_results[epacts_results$pvalue < 5e-8, ]
  
  # Save significant results
  write.csv(significant_snps, output_file, row.names = FALSE)
}

# ======================
# EXECUTE PIPELINE
# ======================

# Example usage of all modules

# Preprocess metabolites
preprocess_metabolites(
  input_path = paste0(input_dir, "metabolites_raw.csv"),
  output_path = paste0(output_dir, "metabolites_preprocessed.csv")
)

# Run DCA
metabolites <- fread(paste0(output_dir, "metabolites_preprocessed.csv"))
run_dca(
  data = metabolites[,-1],
  n_factors = 10,
  output_file = paste0(output_dir, "dca_results_10_factors.rds")
)

# Visualize DCA results
dca_results <- readRDS(paste0(output_dir, "dca_results_10_factors.rds"))
plot_variance_explained(
  dca_results = dca_results,
  output_file = paste0(results_dir, "plots/variance_explained.pdf")
)

# Select top factors
select_factors(
  dca_results = dca_results,
  top_n = 5,
  output_file = paste0(results_dir, "tables/selected_factors.csv")
)

# Correlation analysis
metabolite_correlation(
  data = metabolites[,-1],
  output_file = paste0(results_dir, "tables/metabolite_correlations.csv")
)

# Run EPACTS
run_epacts_pipeline(
  trait = "trait_name",
  chr = 1,
  vcf_dir = paste0(input_dir, "vcf"),
  out_dir = paste0(output_dir, "epacts"),
  pheno_file = paste0(input_dir, "phenotype.ped")
)

# Analyze EPACTS results
analyze_epacts(
  epacts_results_file = paste0(output_dir, "epacts/trait_name_chr1.epacts.gz"),
  output_file = paste0(results_dir, "tables/significant_snps.csv")
)