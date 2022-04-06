library(ggplot2)
library(tidyr)
library(dplyr)

source("r_scripts/gene_expression/02_gene_correlation_matrix.R")

### Remove unecssary files from the environment
rm(annotation_df, basics_gene, gene, heatmap_data)

# Plotting weight against each dpi of infection, dividing into different species 

#switch off all dev devices
while (!is.null(dev.list()))  dev.off() 


jpeg("output_data/gene_expression/04_weight_dpi_primary.jpg", width = 800, 
     height = 600)

# Primary: 
Challenge %>% 
  filter(infection == "primary", !dpi %in% "0", death == "primary") %>%
  drop_na(weight_dpi0, relative_weight) %>%
  group_by("EH_ID") %>%
  ggplot(aes(x = dpi, y = relative_weight, color = primary_infection)) +
  geom_jitter() +
  stat_smooth() +
  labs(x = "Days Post Infection", y = "Relative weight to first day",
       title = "Weight changes during the course of the primary infection")

dev.off()

jpeg("output_data/gene_expression/05_weight_dpi_challenge.jpg", width = 800, 
     height = 600)

Challenge %>% 
  filter(infection == "challenge", !dpi %in% "0", death == "challenge") %>%
  drop_na(weight_dpi0, relative_weight) %>%
  group_by("EH_ID") %>%
  ggplot(aes(x = dpi, y = relative_weight, color = challenge_infection)) +
  geom_jitter() +
  stat_smooth() +
  labs(x = "Days Post Infection", y = "Relative weight to first day",
       title = "Weight changes during the course of infection")

dev.off()


## OPG against dpi 
# Primary infections: 

jpeg("output_data/gene_expression/06_oocysts_primary_dpi.jpg", 
     width = 800, height = 600)

Challenge %>% 
  drop_na(OOC) %>%
  group_by("EH_ID") %>%
  filter(infection == "primary", death == "primary") %>%
  ggplot(aes(x = dpi, y = OOC, color = primary_infection)) +
  geom_point(position = position_jitterdodge()) +
  stat_smooth() +
  labs(x = "Days Post Infection", y = "Oocysts per gram",
       title = "Oocyst shedding in primary infections during the 
       course of infection")

dev.off()

# Challenge infections:

jpeg("output_data/gene_expression/07_oocysts_challenge_dpi.jpg", width = 800, height = 600)

Challenge %>% 
  drop_na(OOC) %>%
  group_by("EH_ID") %>%
  filter(infection == "challenge") %>%
  ggplot(aes(x = dpi, y = OOC, color = challenge_infection)) +
  geom_point(position = position_jitterdodge()) +
  stat_smooth() +
  labs(x = "Days Post Infection", y = "Oocysts per gram",
                     title = "Oocyst shedding in challenge infections during the 
       course of infection")

dev.off()

## Plot the weight loss against the infection intensity

jpeg("output_data/gene_expression/08_weight_intensity_primary.jpg", width = 800, height = 600)

# primary 
Challenge %>% 
  drop_na(delta, max_WL) %>%
  group_by("EH_ID") %>%
  filter(infection == "primary", Eim_MC == "TRUE", death == "primary") %>%
  ggplot(aes(x = delta, y = max_WL, color = primary_infection)) +
  geom_jitter() +
  labs(x = "Maximum weight loss of each mouse", y = "Delta Ct, Infection intensity",
       title = "Maximum Weight loss for each mouse and infection intensity, 
       primary infections")

dev.off()

jpeg("output_data/gene_expression/09_weight_intensity_primary.jpg", width = 800, height = 600)

# challenge 
Challenge %>% 
  drop_na(delta, max_WL) %>%
  group_by("EH_ID") %>%
  filter(infection == "challenge", Eim_MC == "TRUE", death == "challenge") %>%
  ggplot(aes(x = max_WL, y = delta, color = challenge_infection)) +
  geom_jitter() +
  labs(x = "Maximum weight loss of each mouse", y = "Delta Ct, Infection intensity",
       title = "Maximum Weight loss for each mouse and infection intensity, 
       challenge infections") + 
  stat_smooth()

dev.off()


# Now plot the gene expression agains infection intensity
# Use the data frame that has been already cleaned for nas 
gene_expr_delta <- gene_na_omit %>%
  pivot_longer(cols = 8:28, names_to = "Gene", values_to = "gene_expression") %>%
  na.omit(expression) %>%
  ggplot(aes(x = delta, y = gene_expression, color = challenge_infection))

jpeg("output_data/gene_expression/10_gene_expression_intensity", width = 800, height = 600)

gene_expr_delta +
  geom_jitter() +
  facet_wrap(~ Gene, scales = "free") +
  theme_bw()

dev.off()
