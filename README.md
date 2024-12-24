# mdGWAS: Metabolic Dynamics GWAS Analysis Pipeline

This repository contains the implementation of mdGWAS (Metabolic Dynamics GWAS), a framework for identifying SNPs associated with changing correlation patterns of metabolite pairs by integrating dynamic correlation analysis with GWAS.

## Overview

The mdGWAS pipeline consists of the following key components in a single R script (`mdGWAS_pipeline.R`):

1. **Dynamic Correlation Analysis (DCA)**
2. **Local False Discovery Rate (FDR) filtering**
3. **GWAS analysis**
4. **Metabolic pathway analysis**
5. **Functional annotation**
6. **Result visualization**

## Installation

### Prerequisites
- R (>= 3.6.0)
- Required R packages:
  ```r
  install.packages(c("data.table", "ggplot2", "DCA", "coin", "igraph", "locfdr"))
  ```

- **EPACTS** software installed and configured.

### Usage

#### Basic Usage

1. Source the pipeline script:

   ```r
   source("mdGWAS_pipeline.R")
   ```

2. Run the complete analysis:

   ```r
   results <- run_mdgwas(
     metabolite_file = "data/metabolites.txt",
     genotype_file = "data/genotypes.vcf",
     pathway_db = "data/pathway_db.rds",
     n_factors = 10,
     fdr_cutoff = 0.2
   )
   ```

## Pipeline Components

### DCA Core Implementation

```r
dca_witho() # Performs Dynamic Correlation Analysis
```

### Local FDR Calculation

```r
calculate_local_fdr() # Calculates local false discovery rates
```

### Pathway Analysis

```r
analyze_metabolic_pathways() # Performs pathway enrichment analysis
```

### GWAS Analysis

```r
run_epacts() # Executes EPACTS GWAS analysis
```

### Visualization Functions

```r
plot_pathway_network() # Visualizes pathway networks
plot_manhattan() # Creates Manhattan plots
```

## Output

The pipeline generates:

- DCA results
- GWAS summary statistics
- Pathway analysis results
- Visualization plots:
  - Manhattan plots
  - Pathway networks
  - Correlation heatmaps

## Example Results

### Load example results

```r
results <- readRDS("results/mdgwas_results.rds")
```

### Access different components:

```r
dca_results <- results$dca
gwas_results <- results$gwas
pathway_results <- results$pathways
```

## Acknowledgments

This implementation is based on the methodology described in: "Dissecting Genetic Regulation of Metabolic Dynamics"

