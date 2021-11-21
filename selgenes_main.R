#!/usr/bin/env Rscript

# Libraries and imports.
library(data.table)
source("selgenes_functions.R")

# Get input file names.
input_genes = list.files("input_genes/2021-01-18_asm3/", full.names = T)

# Prepare cluster for parallel execution.
cl = cluster()

# For each input file.
for (file in input_genes) {
  selgenes_file(file)
  # To allow other gene types different to the default (protein_coding), use: 
  # e.g. to use both protein_coding and lncRNA gene types
  # selgenes_file(file, genetypes_filter = c("protein_coding", "lncRNA")) 
}

# Stop cluster.
stopCluster(cl)
