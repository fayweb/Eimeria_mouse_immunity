
# positive regulation of regulatory T cell differentiation

#################### CASP1, PRF1, CXCR3, IL6, MUC5AC
##  "positive regulation of interleukin-1 beta production
## positive regulation of inflammatory reslibrary(ggplot2)
library(ggrepel)
library(dplyr)
library(viridis)  # for color-blind-friendly palette
library(scales)
library(cowplot)
library(ggthemes)

options(ggrepel.max.overlaps = Inf)

# read the lab data with pca vectors
lab <- read.csv("output_data/lab_pca")

# change the labels pc1 and pc2 to PC1 / PC2
lab <- lab %>%
  dplyr::rename(PC1 = "pc1", PC2 = "pc2")

# read the variance explained by each gene for the pca 
vpg <- read.csv("output_data/variance_contr_gene_lab")

# Change the first column of the variance contribution of variables to the gene
# names
vpg <- vpg %>%
  dplyr::rename(Variable = vars, PC1 = Dim.1, PC2 = Dim.2)

# add cos2 to lab
lab <- lab %>% mutate(cos2 = lab$PC1^2 + lab$PC2^2)

# Define color palette
color_palette <- c("E_ferrisi" = "#66C2A5", "uninfected" = "#8DA0CB", "E_falciformis" = "#FC8D62")

# PCA graph of individuals
ggplot(lab, aes(x = PC1, y = PC2, color = infection, shape = infection)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray50") + 
  geom_vline(xintercept = 0, linetype = "dotted", color = "gray50") +
  geom_point(size = 3, alpha = 0.8) +
  labs(x = "PC1", y = "PC2", title = "PCA graph of individuals",
       colour = "Current infection", shape ="Current infection") +
  theme_minimal() +
  theme(plot.title = element_text(size = 24, face = "bold"),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = "right") +
  scale_color_manual(values = color_palette) +
  scale_shape_manual(values = c("E_ferrisi" = 17, "uninfected" = 16, "E_falciformis" = 18)) +
  guides(color = guide_legend(override.aes = list(size = 4)))


####### PCA graph of variables

# Add cos2 variable to the dataframe
vpg$cos2 <- with(vpg, PC1^2 + PC2^2)
vpg_labels <- vpg

# Add specific labels for variables
vpg_labels$Additional_Label <- NA
vpg_labels$Additional_Label[vpg_labels$Variable %in% 
                              c("IFNy", "IL.6", "CASP1", "IDO1", "TNF")] <- 
  "*"

# Plotting the factor map with labels
ggplot(vpg_labels, aes(x = PC1, y = PC2)) +
  geom_segment(aes(xend = 0, yend = 0), color = "gray50") +
  geom_point(size = 3) +
  geom_point(data= vpg_labels %>%
               dplyr::filter(Pro_infl == "positive regulation of inflammatory response"),
             pch = 24,
             size=5, 
             colour = "red") +
  geom_label_repel(aes(label = Variable), size = 3, box.padding = 0.5, max.overlaps = 20) +
  coord_equal() +
  xlab("PC1") +
  ylab("PC2") +
  ggtitle("PCA Plot of Variables") +
  theme_minimal() +
  theme(legend.position = "right") +
  theme(plot.title = element_text(size = 18))

# the red triangles signify the variables that are 
#into positive regulation of inflammatory response

####################
# Define color palette
color_palette <- c("E_ferrisi" = "#66C2A5", "uninfected" = "#8DA0CB", "E_falciformis" = "#FC8D62")

# PCA: PC1 predicting weight loss
ggplot(lab, aes(x = PC1, y = WL_max, color = infection, shape = infection)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "black", size = 0.5) +
  labs(x = "PC1", y = "Maximum weight loss", title = "First principal component predicting maximum weight loss during an infection") +
  scale_color_manual(values = color_palette, guide = FALSE) +
  scale_shape_manual(values = c("E_ferrisi" = 16, "uninfected" = 17, "E_falciformis" = 18)) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.key = element_blank()
  )


# PCA: PC1 predicting weight loss with single regression line
ggplot(lab, aes(x = PC1, y = WL_max, color = infection, shape = infection)) +
  geom_point(colour = infection, size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "black", size = 0.5) +
  labs(x = "PC1", y = "Maximum weight loss", title = "First principal component predicting maximum weight loss during an infection") +
  scale_color_manual(values = color_palette, guide = FALSE) +
  scale_shape_manual(values = c("E_ferrisi" = 16, "uninfected" = 17, "E_falciformis" = 18)) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.key = element_blank()
  )

# pc1 predicting weight loss
ggplot(lab, aes(x = PC1, y = WL_max)) +
    geom_point(aes(color = infection, alpha = 0.8, 
                   shape = infection, size = 2)) +
    geom_smooth(method = "lm",formula = y ~ x, 
                se = T, color = "black",
                size = 0.5) +
    labs(x = "PC1", y = "Maximum weight loss", title = 
    "First principal compenonent predicting maximum weight loss during an infection") +
    theme_bw()





# for the model with pc2 and infection
ggplot(lab, aes(x = pc2, y = WL_max)) +
    geom_point(aes(color = infection, alpha = 0.8, 
                   shape = infection, size = 2)) +
    geom_smooth(method = "lm",formula = y ~ x, se = T, 
                color = "black",
                size = 0.5) +
    labs(x = "PC1", y = "Maximum Weight Loss") +
    theme_bw()


### Residual Plots

# calculate residuals for the model with pc1 and pc2
lab$residuals_pc1_pc2 <- resid(model_1_pc1_pc2)

ggplot(lab, aes(x = pc1, y = residuals_pc1_pc2)) +
    geom_point() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(x = "PC1", y = "Residuals") +
    theme_bw()


ggplot(lab, aes(x = pc2, y = residuals_pc1_pc2)) +
    geom_point() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(x = "PC2", y = "Residuals") +
    theme_bw()


ggplot(lab, aes(x = infection, y = residuals_pc1_pc2)) +
    geom_point() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(x = "Infection group", y = "Residuals") +
    theme_bw()


### 3D plots



# First, make sure infection is a factor
lab$infection <- as.factor(lab$infection)

# Then, define the color for each level of infection
color_mapping <- c("E_falciformis" = "salmon", 
                   "E_ferrisi" = "green", 
                   "uninfected" = "blue")

# Now create the scatter plot using this color mapping
ggplot(lab, aes(x = pc1, y = WL_max, color = infection)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE, aes(color = infection)) +
    scale_color_manual(values = color_mapping) +
    labs(x = "PC1", y = "Maximum Weight Loss") +
    theme_bw()



# Now create the scatter plot using this color mapping
ggplot(lab, aes(x = pc2, y = WL_max, color = infection)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE, aes(color = infection)) +
    scale_color_manual(values = color_mapping) +
    labs(x = "PC2", y = "Maximum Weight Loss") +
    theme_bw()






# Create a new color column in your dataframe by mapping 'infection' to your colors
lab$color <- color_mapping[lab$infection]

# 3D scatter plot
scatterplot3d(lab$pc1, lab$pc2, lab$WL_max, pch = 16, color = lab$color,
              xlab = "PC1", ylab = "PC2", zlab = "Maximum Weight Loss")


