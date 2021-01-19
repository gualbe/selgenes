# Test de normalidad

library(data.table)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(stringr)

base_dir = "output_genes/2020-10-05_35genesPaper/"
genes = list.files(path = base_dir, pattern = "*_all_scores.csv")
lgraf = list()
pvalues = c()
difppv = c()

for (gen in genes) {
  # Leer valores de Pi.
  x = fread(paste0(base_dir, gen))
  
  # Calcular p-valor muestral.
  x[, ppvm := as.numeric(.N)]
  for (i in seq(1,x[,.N])) {
    v = x[i, score.p]
    p = x[score.p >= v, .N]
    x[i, ppvm := p / ppvm]
  }
  x[, difppv := abs(ppv - ppvm)]
  fwrite(x, paste0(base_dir, str_replace(gen, ".csv", ""), "_ppvm.csv"))
  difppv = c(difppv, x[ppv <= 0.05, mean(difppv)])
  # difppv = rbindlist(difppv, list(data.table(gen = str_replace(gen, "_all_scores.csv", ""), difppvm05 = x[ppv <= 0.05, mean(difppv)])))
  
  # Test de normalidad muestreando cada 'h' valores de la secuencia (limite test = 5000).
  h = ceiling(length(x[[1]]) / 5000)
  s = x[order(score.p), score.p][seq(1, length(x[[1]]), h)]
  test = shapiro.test(s)
  pvalues = c(pvalues, test$p.value)

  # Agregar grafico de densidad de Pi.
  xlab = paste0(str_replace(gen, "_all_scores.csv", ""), " (", test$p.value, ")")
  lgraf[[length(lgraf)+1]] <- ggdensity(x$score.p, xlab = xlab)
}

grid.arrange(grobs = lgraf)
cat(pvalues)



# ---------------------------------------

# x = fread("output_genes/2020-10-05_35genesPaper/ABL2_all_scores.csv")
# g1 = ggdensity(x$score.p)
# g2 = ggqqplot(x$score.p)

# No se permite mas de 5000 muestras para el test.
# shapiro.test(x$score.p)

# Idea 1: Muestrear aleatoriamente 5000 valores.
# s = sample(x$score.p, size = 5000)
# shapiro.test(s)
# ggdensity(s)
# ggqqplot(s)

# Idea 2: Muestrear cada 4 valores de la secuencia.
# s = x[order(score.p), score.p][seq(1, length(x[[1]]), 4)]
# c = shapiro.test(s)
# ggdensity(s)
# ggqqplot(s)

# Idea 3: Hacer un espejo en los valores por debajo de la media.
# max <- which.max(density(x$score.p)$y)
# m = density(x$score.p)$x[which.max(density(x$score.p)$y)]
# m = 0.5 # x[, mode(score.p)]
# xder = x[score.p >= 0.6, score.p]
# xizq = 1 - xder
# s = c(xizq, xder)
# ggdensity(xder)
# ggqqplot(s)
