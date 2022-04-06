library(ggplot2)
library(tidyr)
library(dplyr)

source("r_scripts/gene_expression/02_gene_correlation_matrix.R")

### Remove unecssary files from the environment
rm(annotation_df, basics_gene, gene, heatmap_data, gene)

# Plotting weight against each dpi of infection, dividing into different species 

#switch off all dev devices
while (!is.null(dev.list()))  dev.off() 

jpeg("output_data/gene_expression/04_weight_dpi.jpg", width = 1400, height = 1000)

Challenge %>% 
  filter(infection == "challenge", !dpi %in% "0") %>%
  drop_na(weight_dpi0, relative_weight) %>%
  group_by("EH_ID") %>%
  ggplot(aes(x = dpi, y = relative_weight, color = challenge_infection)) +
  geom_jitter() +
  stat_smooth()

dev.off()

# Now plot the gene expression agains infection intensity
gene_expr_delta <- gene_na_omit %>%
  pivot_longer(cols = 8:28, names_to = "Gene", values_to = "gene_expression") %>%
  na.omit(expression) %>%
  ggplot(aes(x = delta, y = gene_expression, color = challenge_infection))

jpeg("output_data/gene_expression/05_gene_expression_intensity", width = 800, height = 600)

gene_expr_delta +
  geom_jitter() +
  facet_wrap(~ Gene, scales = "free") +
  theme_bw()

dev.off()

ggplot(data = gene_intensity, aes(x = c))


ggplot(gene_na_omit, aes(ends_with("_N")))


function(x) {}
ggplot(gene_intensity, (aes()))