# ======================
# Load Required Packages
# ======================
library(data.table)
library(ggplot2)
library(DCA)
library(coin)
library(igraph)
library(locfdr)

# ======================
# Utility Functions
# ======================
normrow <- function(x) {
  # Standardize each row
  x <- sweep(x, 1, rowMeans(x))
  x <- sweep(x, 1, sqrt(rowSums(x^2)/ncol(x)), "/")
  return(x)
}

# ======================
# DCA Core Implementation
# ======================
dca_witho <- function(array, top.pairs.prop = 0.95, max.pairs = 1e+06, n.fac = 10, 
                     sumabsv = sqrt(max.pairs)/10, normalization = "standardize", 
                     method = "PCA") {
  # Data standardization
  if (normalization == "standardize") {
    array <- normrow(array)
  }
  
  # Calculate correlations
  n.pairs <- min(choose(nrow(array), 2) * (1 - top.pairs.prop), max.pairs)
  ccc.abs <- cor(abs(t(array)))
  abs.ccc <- abs(cor(t(array)))
  ccc.diff <- ccc.abs - abs.ccc
  rm(ccc.abs, abs.ccc); gc()
  
  # Select significant pairs
  pos <- which(ccc.diff > quantile(ccc.diff, 1 - 2 * n.pairs/nrow(array)/nrow(array)), 
               arr.ind = T)
  pos <- pos[which(pos[, 1] > pos[, 2]), ]
  rm(ccc.diff); gc()
  
  # Calculate dynamic correlations
  b <- array[pos[, 1], ] * array[pos[, 2], ]
  for (i in 1:nrow(b)) b[i, ] <- coin:::normal_trafo(b[i, ])
  
  # PCA analysis
  ccc <- t(b) %*% b
  e <- eigen(ccc)
  fac <- e$vec[, 1:n.fac]
  
  # Calculate explained variance
  proj.len.unrotated <- apply((b %*% e$vec)^2, 2, sum)
  ss.proj.unrotated.full <- proj.len.unrotated/sum(proj.len.unrotated)
  
  proj.len.unrotated <- apply((b %*% fac)^2, 2, sum)
  ss.proj.unrotated <- proj.len.unrotated/sum(proj.len.unrotated)
  
  # Factor rotation
  if (nrow(b) > 1e+05) 
    b <- b[sample(nrow(b), 1e+05, replace = FALSE), ]
  proj <- b %*% fac
  rot <- diag(ncol(fac))
  try(rot <- varimax(proj)$rot)
  fac2 <- fac %*% rot
  
  proj <- b %*% fac2
  proj.len <- apply(proj^2, 2, sum)
  o <- order(-proj.len)
  fac2 <- fac2[, o]
  ss.proj <- proj.len[o]/sum(proj.len)
  
  # Return results
  return(list(
    fac = fac, 
    rotated = fac2, 
    ss.proj = ss.proj,
    ss.proj.unrotated = ss.proj.unrotated,
    ss.proj.unrotated.full = ss.proj.unrotated.full,
    order = o,
    metabolite_pairs = pos
  ))
}

# ======================
# Local FDR Calculation
# ======================
calculate_local_fdr <- function(la_scores, cutoff = 0.2) {
  # Fit mixture distribution
  null_scores <- la_scores[la_scores < quantile(la_scores, 0.95)]
  null_mean <- mean(null_scores)
  null_sd <- sd(null_scores)
  
  # Calculate densities
  f0 <- dnorm(la_scores, null_mean, null_sd)
  f <- density(la_scores)
  
  # Calculate local FDR
  pi0 <- min(1, length(null_scores)/length(la_scores))
  fdr <- pi0 * f0 / approx(f$x, f$y, xout = la_scores)$y
  
  # Return significant pairs
  return(la_scores[fdr < cutoff])
}

# ======================
# Metabolic Pathway Analysis
# ======================
analyze_metabolic_pathways <- function(metabolite_pairs, pathway_db) {
  # Hypergeometric test
  enrich_test <- function(pathway_genes, metabolite_set, background) {
    overlap <- intersect(pathway_genes, metabolite_set)
    return(phyper(length(overlap)-1, 
                 length(pathway_genes),
                 length(background)-length(pathway_genes),
                 length(metabolite_set),
                 lower.tail=FALSE))
  }
  
  # Analyze each pathway
  results <- data.frame(
    pathway = names(pathway_db),
    genes = sapply(pathway_db, function(x) paste(x$genes, collapse=";")),
    p_value = sapply(pathway_db, function(x) 
      enrich_test(x$genes, metabolite_pairs, unique(unlist(pathway_db))))
  )
  
  # FDR correction
  results$q_value <- p.adjust(results$p_value, method = "BH")
  
  return(results)
}

# ======================
# GWAS Analysis
# ======================
run_epacts <- function(phenotype, genotype, covariates = NULL) {
  # Build EPACTS command
  cmd <- paste(
    "epacts single",
    "--vcf", genotype,
    "--pheno", phenotype,
    "--test q.emmax",
    "--out", "results"
  )
  
  if(!is.null(covariates)) {
    cmd <- paste(cmd, "--cov", covariates)
  }
  
  # Run EPACTS
  system(cmd)
  
  # Read results
  results <- fread("results.epacts")
  
  return(results)
}

# ======================
# Gene Set Enrichment Analysis
# ======================
perform_gsea <- function(gene_list, gene_sets) {
  # Calculate enrichment score
  calc_es <- function(gene_set, ranked_genes) {
    hits <- cumsum(ranked_genes %in% gene_set)
    misses <- cumsum(!(ranked_genes %in% gene_set))
    es <- max(hits/length(gene_set) - misses/(length(ranked_genes)-length(gene_set)))
    return(es)
  }
  
  # GSEA analysis
  results <- lapply(gene_sets, function(gs) {
    es <- calc_es(gs$genes, gene_list)
    perms <- replicate(1000, calc_es(gs$genes, sample(gene_list)))
    p_value <- mean(perms >= es)
    return(list(es=es, p_value=p_value))
  })
  
  return(results)
}

# ======================
# Result Visualization
# ======================
plot_pathway_network <- function(pathway_results, threshold=0.05) {
  # Build network
  sig_paths <- pathway_results$q_value < threshold
  g <- graph_from_data_frame(
    d = data.frame(
      from = pathway_results$pathway1[sig_paths],
      to = pathway_results$pathway2[sig_paths],
      weight = -log10(pathway_results$q_value[sig_paths])
    ),
    directed = FALSE
  )
  
  # Plot network
  plot(g,
       layout = layout_with_fr,
       vertex.size = degree(g),
       edge.width = E(g)$weight,
       main = "Metabolic Pathway Network")
}

plot_manhattan <- function(gwas_results, threshold = 5e-8) {
  # Manhattan plot
  manhattan <- ggplot(gwas_results, 
                     aes(x = POS, y = -log10(P), color = factor(CHR))) +
    geom_point() +
    geom_hline(yintercept = -log10(threshold), linetype = "dashed") +
    facet_grid(. ~ CHR, scales = "free_x", space = "free_x") +
    theme_bw() +
    theme(axis.text.x = element_blank())
  
  print(manhattan)
}

# ======================
# Main Analysis Pipeline
# ======================
run_mdgwas <- function(metabolite_file, genotype_file, pathway_db, 
                      n_factors = 10, fdr_cutoff = 0.2) {
  # 1. Load data
  message("Loading data...")
  metabolites <- fread(metabolite_file)
  
  # 2. DCA analysis
  message("Running DCA...")
  dca_results <- dca_witho(
    t(as.matrix(metabolites[,-1])), 
    n.fac = n_factors
  )
  
  # 3. Local FDR filtering
  message("Calculating local FDR...")
  sig_pairs <- calculate_local_fdr(
    dca_results$metabolite_pairs, 
    cutoff = fdr_cutoff
  )
  
  # 4. GWAS analysis
  message("Running GWAS...")
  gwas_results <- lapply(1:n_factors, function(i) {
    run_epacts(
      phenotype = dca_results$rotated[,i],
      genotype = genotype_file
    )
  })
  
  # 5. Pathway analysis
  message("Analyzing pathways...")
  pathway_results <- analyze_metabolic_pathways(
    sig_pairs,
    pathway_db
  )
  
  # 6. Functional annotation
  message("Performing GSEA...")
  gsea_results <- perform_gsea(
    unique(unlist(lapply(gwas_results, function(x) x$SNP))),
    pathway_db
  )
  
  # 7. Generate visualizations
  message("Generating plots...")
  plot_pathway_network(pathway_results)
  lapply(gwas_results, plot_manhattan)
  
  # 8. Return results
  return(list(
    dca = dca_results,
    significant_pairs = sig_pairs,
    gwas = gwas_results,
    pathways = pathway_results,
    gsea = gsea_results
  ))
}

# ======================
# Usage Example
# ======================
# Define input/output paths
input_dir <- "data/input/"
output_dir <- "data/output/"

# Run analysis
results <- run_mdgwas(
  metabolite_file = paste0(input_dir, "metabolites.txt"),
  genotype_file = paste0(input_dir, "genotypes.vcf"),
  pathway_db = readRDS(paste0(input_dir, "pathway_db.rds")),
  n_factors = 10,
  fdr_cutoff = 0.2
)

# Save results
saveRDS(results, file = paste0(output_dir, "mdgwas_results.rds"))