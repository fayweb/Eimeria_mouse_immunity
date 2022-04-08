library(dplyr)
library(tidyr)
library(lattice)
library(ggplot2)
library(GGally)
library(scatterplot3d)

source("r_scripts/gene_expression/03_gene_lab_general_plotting.R")

### Multivariate Probability Distributions in R

# exploring the data frame
glimpse(gene_na_omit) 
dim(gene_na_omit)

# calculate the mean of columns
colMeans(gene_na_omit[6:28], na.rm = TRUE)

# mean by infection type
by(gene_na_omit[6:28], gene_na_omit$challenge_infection, colMeans, na.rm = TRUE)

# Calculate the variance-covariance matrix of the genes
var.gene <- var(gene_na_omit[8:28], na.rm = TRUE)

# Round the matrix values to two decimal places 
round(var.gene, 2)

# data frame including only thhe genes 
GN <- gene_na_omit[8:28]

# Calculate the correlation matrix of the genes
cor.gene <- cor(GN %>% na.omit())

# Round the matrix to two decimal places 
round(cor.gene, 2)

# Plot the correlations 
# Visualizing correlates
# Method = eelipse: high correlation = thin circle 
# low correlation = thick circle
corrplot(cor.gene, method = "ellipse", tl.col = "black")

# Scatter plot matrix using the base R plot function
#selecting genes with high correlation to eachother
graphics.off()
pairs(GN[13:15]) #if you run into error then grow the viewing window

# Scatter plot matrix with lattice  
splom( ~ GN[13:15]) #doesn't work out nicely

# Scatter plot matrix colored by groups
splom( ~ gene_na_omit[, 20:22], pch = 16, col = gene_na_omit$challenge_infection)
#not working

#switch off all dev devices
while (!is.null(dev.list()))  dev.off() 

# Produce a matrix of plots 
ggpairs(data = gene_na_omit, columns = 20:22, aes(color = challenge_infection))


# Plot the gene variables]
#label the challenge infections with colors
colors <- c("red","green","blue")
dfcolor <- as.data.frame(cbind(unique(gene_na_omit$challenge_infection), colors)) %>%
    rename(challenge_infection = V1)

gene_na_omit <- gene_na_omit %>% left_join(dfcolor, by = "challenge_infection")

scatterplot3d(gene_na_omit[, c(20, 21, 22)], color = gene_na_omit$colors)

# Play around with multivariate normal distribution



