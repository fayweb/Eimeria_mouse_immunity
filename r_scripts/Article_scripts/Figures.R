library(ggplot2)

# read the lab data with pca vectors
lab <- read.csv("output_data/lab_pca")

# read the variance explained by each gene for the pca 
vpg <- read.csv("output_data/variance_contr_gene_lab")

# pc1 predicting weight loss
ggplot(lab, aes(x = pc1, y = WL_max)) +
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


