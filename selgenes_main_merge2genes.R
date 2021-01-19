#!/usr/bin/env Rscript

# Libraries and imports.
library(data.table)
source("selgenes_functions.R")

# Merge two genes bidirectionally.
mergeTwoGenes("input_genes/2020-06-20_SMN/SMN1.tsv", "SMN1", "SMN2", "SMN12")
mergeTwoGenes("input_genes/2020-06-20_SMN/SMN2.tsv", "SMN1", "SMN2", "SMN21")
input_genes = c("input_genes/2020-06-20_SMN/SMN12.tsv", "input_genes/2020-06-20_SMN/SMN21.tsv")

# Prepare cluster for parallel execution.
cl = cluster()

# For each input file.
for (file in input_genes) {
  selgenes_file(file)
}

# Stop cluster.
stopCluster(cl)

# Merge correlators of two genes (without repeated genes).
mergeTwoResults("output_genes/2020-06-20_SMN", "SMN12", "SMN21", "SMN")


