library(ggplot2)

source("r_scripts/gene_expression/02_gene_correlation_matrix.R")

gene_intensity <- gene_na_omit %>% select(c("EH_ID", "delta", ends_with("_N"))) %>%
  pivot_longer(cols = 3:23, names_to = "Gene", values_to = "expression")

ggplot(gene_na_omit, aes(ends_with("_N")))


function(x) {}
ggplot(gene_intensity, (aes()))