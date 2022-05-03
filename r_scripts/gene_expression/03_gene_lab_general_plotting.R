library(ggplot2)
library(tidyr)
library(dplyr)

source("r_scripts/gene_expression/02_gene_correlation_matrix.R")

### Remove unecssary files from the environment
rm(annotation_df, basics_gene, gene, heatmap_data)

# Plotting weight against each dpi of infection, dividing into different species 

#switch off all dev devices
while (!is.null(dev.list()))  dev.off() 


png("output_data/gene_expression/04.01_weight_dpi_primary.png", res =  100)

# Primary:
Challenge$dpi <- as.double(Challenge$dpi)
Challenge %>%
    mutate(weight_dpi1 = case_when(
        dpi == 1 ~ Challenge$weight,
        TRUE ~ as.integer(NA)))

typeof(Challenge$weight)
    mutate(relative_weight = weight/weight)
  
  ggplot(aes(x = dpi, y = relative_weight, color = primary_infection)) +
  geom_jitter() +
  stat_smooth() +
  labs(x = "Days Post Infection", y = "Relative weight to first day",
       title = "Weight changes during the course of the primary infection") +
    theme_bw()

dev.off()

# Primary /mouse stains: 
png("output_data/gene_expression/04.02_weight_dpi_primary_mouse_strain.png", width = 800, 
     height = 600)

Challenge %>% 
    filter(infection == "primary", !dpi %in% "0", death == "primary") %>%
    drop_na(weight_dpi0, relative_weight) %>%
    group_by("EH_ID") %>%
    ggplot(aes(x = dpi, y = relative_weight, color = primary_infection)) +
    geom_jitter() +
    stat_smooth() +
    facet_wrap(~ mouse_strain) +
    labs(x = "Days Post Infection", y = "Relative weight to first day",
         title = "Weight changes during the course of the primary infection") +
    theme_bw()

dev.off()

# Challenge:
png("output_data/gene_expression/05.01_weight_dpi_challenge.png", width = 800, 
     height = 600)

Challenge %>% 
  filter(infection == "challenge", !dpi %in% "0", death == "challenge") %>%
  drop_na(weight_dpi0, relative_weight) %>%
  group_by("EH_ID") %>%
  ggplot(aes(x = dpi, y = relative_weight, color = challenge_infection)) +
  geom_jitter() +
  stat_smooth() +
  labs(x = "Days Post Infection", y = "Relative weight to first day",
       title = "Weight changes during the course of infection") +
    theme_bw()

dev.off()

# Challenge /mouse_stains:
png("output_data/gene_expression/05.02_weight_dpi_challenge_mouse_strains.png", width = 800, 
     height = 600)

Challenge %>% 
    filter(infection == "challenge", !dpi %in% "0", death == "challenge") %>%
    drop_na(weight_dpi0, relative_weight) %>%
    group_by("EH_ID") %>%
    ggplot(aes(x = dpi, y = relative_weight, color = challenge_infection)) +
    geom_jitter() +
    stat_smooth() +
    labs(x = "Days Post Infection", y = "Relative weight to first day",
         title = "Weight changes during the course of infection") +
    theme_bw()  +
    facet_wrap(~ mouse_strain) 

dev.off()

## OPG against dpi 
# Primary infections: 

png("output_data/gene_expression/06.01_oocysts_primary_dpi.png", 
     res =  100)

Challenge %>% 
  drop_na(OOC) %>%
  group_by("EH_ID") %>%
  filter(infection == "primary", death == "primary") %>%
  ggplot(aes(x = dpi, y = OOC, color = primary_infection)) +
  geom_point(position = position_jitterdodge()) +
  stat_smooth() +
    scale_y_log10() +
  labs(x = "Days Post Infection", y = "Oocysts per gram",
       title = "Oocyst shedding in primary infections during the 
       course of infection") +
    theme_bw()

dev.off()

## OPG against dpi 
# Primary infections: mouse_strain

png("output_data/gene_expression/06.02_oocysts_primary_dpi.png", 
     res =  100)

Challenge %>% 
    drop_na(OOC) %>%
    group_by("EH_ID") %>%
    filter(infection == "primary", death == "primary") %>%
    ggplot(aes(x = dpi, y = OOC, color = primary_infection)) +
    geom_point(position = position_jitterdodge()) +
    scale_y_log10() +
    
    labs(x = "Days Post Infection", y = "Oocysts per gram",
         title = "Oocyst shedding in primary infections during the 
       course of infection") +
    theme_bw() +
    facet_wrap(~ mouse_strain, scales = "free")

dev.off()

# Challenge infections:

png("output_data/gene_expression/07.01_oocysts_challenge_dpi.png", res = 50)

Challenge %>% 
  drop_na(OOC) %>%
  group_by("EH_ID") %>%
  filter(infection == "challenge") %>%
  ggplot(aes(x = dpi, y = OOC, color = challenge_infection)) +
  geom_point(position = position_jitterdodge()) +
  stat_smooth() +
  labs(x = "Days Post Infection", y = "Oocysts per gram",
                     title = "Oocyst shedding in challenge infections during the 
       course of infection") +
    theme_bw()

dev.off()

# Challenge infections: mouse_strain

png("output_data/gene_expression/07.01_oocysts_challenge_dpi_mouse_strain.png", res =  100)

Challenge %>% 
    drop_na(OOC) %>%
    group_by("EH_ID") %>%
    filter(infection == "challenge") %>%
    ggplot(aes(x = dpi, y = OOC, color = challenge_infection)) +
    geom_point(position = position_jitterdodge()) +
    stat_smooth() +
    labs(x = "Days Post Infection", y = "Oocysts per gram",
         title = "Oocyst shedding in challenge infections during the 
       course of infection") +
    theme_bw() +
    facet_wrap(~ mouse_strain, scales = "free")

dev.off()

## Plot the weight loss against the infection intensity

png("output_data/gene_expression/08.01_weight_intensity_primary.png", res =  100)

# primary 
Challenge %>% 
  drop_na(delta, max_WL) %>%
  group_by("EH_ID") %>%
  filter(infection == "primary", Eim_MC == "TRUE", death == "primary") %>%
  ggplot(aes(x = delta, y = max_WL, color = primary_infection)) +
  geom_jitter() +
  labs(x = "Delta Ct, Infection intensity", y = "Maximum weight loss of each mouse",
       title = "Maximum Weight loss for each mouse and infection intensity, 
       primary infections") +
    theme_bw()

dev.off()

## Plot the weight loss against the infection intensity_mouse strain

png("output_data/gene_expression/08.02_weight_intensity_primary_mouse_strain.png", res =  100)

# primary 
Challenge %>% 
    drop_na(delta, max_WL) %>%
    group_by("EH_ID") %>%
    filter(infection == "primary", Eim_MC == "TRUE", death == "primary") %>%
    ggplot(aes(x = delta, y = max_WL, color = primary_infection)) +
    geom_jitter() +
    labs(x = "Delta Ct, Infection intensity", y = "Maximum weight loss of each mouse",
         title = "Maximum Weight loss for each mouse and infection intensity, 
       primary infections") +
    theme_bw() +
    facet_wrap(~ mouse_strain, scales = "free")

dev.off()



# challenge 
png("output_data/gene_expression/09.01_weight_intensity_challenge.png", res =  100)

Challenge %>% 
  drop_na(delta, max_WL) %>%
  group_by("EH_ID") %>%
  filter(infection == "challenge", Eim_MC == "TRUE", death == "challenge") %>%
  ggplot(aes(x = delta, y = max_WL, color = challenge_infection)) +
  geom_jitter() +
  labs(x = "Delta Ct, Infection intensity", y = "Maximum weight loss of each mouse",
       title = "Maximum Weight loss for each mouse and infection intensity, 
       challenge infections") + 
  stat_smooth() +
    theme_bw()

dev.off()

# challenge: mouse_strain
png("output_data/gene_expression/09.02_weight_intensity_challenge_mouse_strain.png", res =  100)

Challenge %>% 
    drop_na(delta, max_WL) %>%
    group_by("EH_ID") %>%
    filter(infection == "challenge", Eim_MC == "TRUE", death == "challenge") %>%
    ggplot(aes(x = delta, y = max_WL, color = challenge_infection)) +
    geom_jitter() +
    labs(x = "Delta Ct, Infection intensity", y = "Maximum weight loss of each mouse",
         title = "Maximum Weight loss for each mouse and infection intensity, 
       challenge infections") + 
    stat_smooth() +
    theme_bw() +
    facet_wrap(~ mouse_strain, scales = "free")

dev.off()


# Now plot the gene expression agains infection intensity
# Use the data frame that has been already cleaned for nas 
gene_expr_delta <- gene_na_omit %>%
  pivot_longer(cols = 8:28, names_to = "Gene", values_to = "gene_expression") %>%
  na.omit(expression) %>%
    group_by(EH_ID) %>%
  ggplot(aes(x = delta, y = gene_expression, color = challenge_infection)) 

png("output_data/gene_expression/10.01_gene_expression_intensity.png", 
     res =  100)

gene_expr_delta +
  geom_jitter() +
  facet_wrap(~ Gene, scales = "free") +
  theme_light() +
    labs(x = "Delta Ct, Infection intensity", y = "Gene expression",
         title = "Gene expression in response to infection intensity") +
    theme_bw()

dev.off()

png("output_data/gene_expression/10.02_gene_expression_intensity_mouse_strains.png", 
     res =  100)

gene_expr_delta +
    geom_jitter() +
    facet_wrap(~ Gene + mouse_strain, scales = "free") +
    theme_light() +
    labs(x = "Delta Ct, Infection intensity", y = "Gene expression",
         title = "Gene expression in response to infection intensity") +
    theme_bw()

dev.off()

png("output_data/gene_expression/11_gene_expression_eimeria_boxplot.png", 
     res =  50)

gene_na_omit %>%
    group_by(EH_ID) %>%
  pivot_longer(cols = 8:28, names_to = "Gene", values_to = "gene_expression") %>%
  na.omit(expression) %>%
  ggplot(aes(x = challenge_infection, y = gene_expression, color = challenge_infection)) + 
  geom_boxplot() +
  geom_jitter() +
  facet_wrap(~ Gene) +
  theme_bw() +
    labs(x = "Infection groups, E64 = E. ferrisi, E88 = E.falciformis, 
         UNI = Uninfected", y = "Gene expression",
         title = "Gene expression in response to infection group") 

dev.off()





