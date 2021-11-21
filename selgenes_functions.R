library(compiler)
library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(doParallel)
library(psych)


# Prepare cluster for parallel execution.
cluster <- function () {
  if (Sys.info()['sysname'] == "Windows") {
    cl<-makeCluster(spec = as.numeric(detectCores()), type="PSOCK")
    clusterExport(cl, c('data.table', 'metric_SFCS') )
  } else {
    cl<-makeCluster(spec = as.numeric(detectCores()), type="FORK")
  }
  registerDoParallel(cl)
  return(cl)
}

selgenes_file <- function (input_file, genetypes_filename = c("data/biomart_human_type.tsv", "data/biomart_celegans_type.tsv"), genetypes_filter = c("protein_coding")) {
  
  # Extract, transform and load gene fold changes from experiments.
  data = etl_experiments(input_file, genetypes_filename, genetypes_filter)
  xgene = data$xgene
  
  # Select genes by using all available methods.
  res = selgenes(data, xgene)
  
  # Write results into the output file.
  dir.create(dirname(sub("input", "output", input_file)))
  output_path_basename = sub("input", "output", sub(".tsv", "", input_file))
  fwrite(data.table(res$direct), paste0(output_path_basename, "_selgenes_direct.csv"), row.names = F, col.names = F)
  fwrite(data.table(res$inverse), paste0(output_path_basename, "_selgenes_inverse.csv"), row.names = F, col.names = F)
  fwrite(data.table(res$direct_scores), paste0(output_path_basename, "_selgenes_direct_scores.csv"), row.names = F, col.names = T)
  fwrite(data.table(res$inverse_scores), paste0(output_path_basename, "_selgenes_inverse_scores.csv"), row.names = F, col.names = T)
  fwrite(data.table(res$scores), paste0(output_path_basename, "_all_scores.csv"), row.names = F, col.names = T)
  
  return(res)
}

etl_experiments <- function (filename, genetypes_filename, genetypes_filter) {
  # Read fold change values from a experiments file.
  data = fread(filename, header=T, blank.lines.skip=T)
  
  # Ignore the three first rows. 
  data = data[V1 != ""]  ## Antes: data[4:.N] ## Corregido el 2021-01-19 (Gualberto)
  
  # Convert data columns to numeric.
  data[,3:ncol(data)] = data[, lapply(.SD, as.numeric), .SDcols = 3:ncol(data)]
  
  # Determine the xgene based upon file name.
  xgene = tools::file_path_sans_ext(basename(filename))
  if (!startsWith(xgene, "ENSG")) {
    xgene = data[V2 == xgene, V1]  
  }
  if (length(xgene) == 0)
    stop(paste("The gene name", tools::file_path_sans_ext(basename(filename)), "is not included in the experiments file."))
  
  # Rename the column of the gene name.
  setnames(data, "V1", "Gene")
  data[, V2 := NULL]
  
  # Filter genes with type "protein_coding".
  i = 1; ok = FALSE
  while (!ok & i <= length(genetypes_filename)) {
    browser()
    types = fread(genetypes_filename[i])
    types = types[, `Gene type`:=as.factor(`Gene type`)]
    coding_genenames = types[`Gene type` %in% genetypes_filter]$`Gene stable ID`  ## Gualberto 2021-11-21 (allowing any multiple gene types as filter)
    data_filtered = data[Gene %in% coding_genenames]
    if (nrow(data_filtered) > 0) ok = TRUE
    i = i + 1
  }
  if (ok) data = data_filtered
  else stop(paste("No protein-coding genes were found in", filename))

  # Determination of duplicated genes.
  dup = data[duplicated(data[,.(Gene)]),Gene]
  genes = unique(data[,.(Gene)])
  
  # Averaging fold change values for repeated genes.
  data = data[, lapply(.SD, mean), by=Gene, .SDcols = 3:ncol(data)]
  
  # Transpose the data table (experiments -x- genes).
  datat = as.data.table(t(data[,-1]))
  setnames(datat, names(datat), t(genes))
  
  return(list(exp2genes = datat, genes2exp = data, duplicated = dup, different_genes = genes, xgene = xgene))
  
}

mergeTwoGenes <- cmpfun(function (file, gene1, gene2, output_genename) {
  # Read file
  data = fread(file, header=T)
  
  # Sum fold changes of two genes (for the same experiments).
  g1 = data[V2 == gene1]
  g2 = data[V2 == gene2]
  s = as.numeric(g1[,3:(ncol(g1))]) + as.numeric(g2[,3:ncol(g2)])
  
  # Replace new fold changes into the dataset for the first gene.
  for (i in 3:ncol(g1)) {
    data[V2==gene1, i] = s[i-2]
  }
  
  # Remove the row corresponding to the second gene.
  data = data[V2!=gene2]
  
  # Relabel the first gene.
  data[V2==gene1, V2 := output_genename]
  
  # Write results.
  fwrite(data, paste0(dirname(file), "/", output_genename, ".tsv"), quote = F, sep = "\t")
  
})

mergeTwoResults <- cmpfun(function (basepath, gene1, gene2, outgene) {
  
  # Step 1. Merge direct scores.
  f1 = fread(paste0(basepath, "/", gene1, "_selgenes_direct_scores.csv"))
  f2 = fread(paste0(basepath, "/", gene2, "_selgenes_direct_scores.csv"))
  direct = rbindlist(list(f1,f2))
  direct = direct %>% group_by(gene) %>% summarise_all(mean) %>% arrange(ppv)
  fwrite(direct, paste0(basepath, "/", outgene, "_selgenes_direct_scores.csv"))
  
  # Step 2. Merge inverse scores.
  f1 = fread(paste0(basepath, "/", gene1, "_selgenes_inverse_scores.csv"))
  f2 = fread(paste0(basepath, "/", gene2, "_selgenes_inverse_scores.csv"))
  inverse = rbindlist(list(f1,f2))
  inverse = inverse %>% group_by(gene) %>% summarise_all(mean) %>% arrange(npv)
  fwrite(inverse, paste0(basepath, "/", outgene, "_selgenes_inverse_scores.csv"))
  
  # Step 3. Store both direct and inverse gene names separately.
  fwrite(direct %>% select(gene), paste0(basepath, "/", outgene, "_selgenes_direct.csv"))
  fwrite(inverse %>% select(gene), paste0(basepath, "/", outgene, "_selgenes_inverse.csv"))
  
})

selgenes <- function(data, xgene) {
  res = selgenes_method_SFCS(data, xgene)
  return(res)
}

metric_SFCS <- cmpfun(function (fc_xgene, fc_other) {
  r = fc_xgene * fc_other
  p = r[r[[1]] > 0,.N] / r[,.N]
  n = r[r[[1]] < 0,.N] / r[,.N]
  z = 1 - p - n
  if (length(which(r[[1]] < 0)) > 0) {
    ppen = fc_other[which(r[[1]] < 0), sum(abs(.SD))]
  } else {
    ppen = 0
  }
  if (length(which(r[[1]] > 0)) > 0) {
    npen = fc_other[which(r[[1]] > 0), sum(abs(.SD))]
  } else {
    npen = 0
  }
  return(data.table(p, n, z, ppen, npen))
})

empiric_pvalue <- function (distr) {
  
}

selgenes_method_SFCS <- function(data_pack, xgene, method = "pvalue", sd_factor = 2) {
  
  # Extract only the used variables from the data_pack.
  data = data_pack$genes2exp
  datat = data_pack$exp2genes

  # Assess the score for each gene.
  score <- foreach(col=1:ncol(datat), .combine=rbind) %dopar% {   ## dopar
    metric_SFCS(datat[,xgene,with=F], datat[,col,with=F])
  }
  
  # Build the result table.
  ppv = 1 - pnorm(scale(score$p))
  npv = 1 - pnorm(scale(score$n))
  ppenpv = 1 - pnorm(scale(score$ppen))
  npenpv = 1 - pnorm(scale(score$npen))
  genes_table = data.table(gene = names(datat), score = score, ppv = ppv, npv = npv, ppenpv = ppenpv, npenpv = npenpv)
  setnames(genes_table, c("ppv.V1","npv.V1","ppenpv.V1","npenpv.V1"), c("ppv","npv","ppenpv","npenpv"))
  genes_table = genes_table[gene != xgene]    # [score.p != 1]
  
  # Selection of relevant related genes to xgene.
  if (method == "pvalue") {
    # (by p-value).
    sel_genes_direct = genes_table[ppv <= 0.05]   #  & ppenpv <= 0.05
    q3_ppen = quantile(sel_genes_direct$score.ppen)[4]
    sel_genes_direct = sel_genes_direct[score.ppen < q3_ppen]
    
    sel_genes_inverse = genes_table[npv <= 0.05]   #  & npenpv <= 0.05
    q3_npen = quantile(sel_genes_inverse$score.npen)[4]
    sel_genes_inverse = sel_genes_inverse[score.npen < q3_npen]
  } else {
    # (by mean + sd * sd_factor).
    threshold_direct = mean(score$p) + sd_factor * sd(score$p)
    threshold_inverse = mean(score$n) + sd_factor * sd(score$n)
    sel_genes_direct = genes_table[score.p >= threshold_direct]
    sel_genes_inverse = genes_table[score.n >= threshold_inverse]
  }

  # Sort by their scores.
  sel_genes_direct = sel_genes_direct[order(-score.p)]
  sel_genes_inverse = sel_genes_inverse[order(-score.n)]
  
  return(list(direct = sel_genes_direct[,gene], inverse = sel_genes_inverse[,gene], scores = genes_table, 
              direct_scores = sel_genes_direct, inverse_scores = sel_genes_inverse))
  
}
