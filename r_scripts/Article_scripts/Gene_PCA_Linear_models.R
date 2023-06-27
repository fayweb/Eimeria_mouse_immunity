library(tidyverse)
library(dplyr)
library(stringr)
library(FactoMineR)
library(reshape2)
library(corrplot)
library(factoextra)
library(lmtest)
library(ggpubr)
library(janitor)
library(pheatmap)
library(visdat)
library(scatterplot3d)
library(clusterProfiler) # gene enrichment analysis
library(org.Mmu.eg.db) # gene ids identifiers


# Load the normalized and imputed data set
hm <- read.csv("output_data/2.imputed_MICE_data_set.csv")


# Vectors for selecting the relevant immune genes
Gene_lab   <- c("IFNy", "CXCR3", "IL.6", "IL.13",
                "IL1RN","CASP1", "CXCL9", "IDO1", "IRGM1", "MPO", 
                "MUC2", "MUC5AC", "MYD88", "NCR1", "PRF1", "RETNLB", "SOCS1", 
                "TICAM1", "TNF") # "IL.12", "IRG6")

# Data cleaning and preparation

## Select the laboratory data. As we have duplicate data for each mice with 
## the same gene expression values, we select the mice once, using the mesenterial
## lymphnode data as a separator ("mLN")
lab <- hm %>%
    dplyr::filter(origin == "Lab", Position == "mLN") 

## check for duplicates
lab <- unique(lab)

## select the mice labels and the immune genes
gene <- lab %>%
    dplyr::select(c(Mouse_ID, all_of(Gene_lab)))

## duplicates?
genes <- unique(gene)

## remove the mice labes
genes <- genes[, -1]

#remove rows with only nas
genes <- genes[,colSums(is.na(genes))<nrow(genes)]

#remove colums with only nas 
genes <- genes[rowSums(is.na(genes)) != ncol(genes), ]

#select same rows in the first table
gene <- gene[row.names(genes), ]

# we need to change the  in challenge infections to a factor
lab$infection <- as.factor(lab$infection)
lab$MC.Eimeria <- as.factor(lab$MC.Eimeria)

# Here I create a new column, where we get the actual infection status
# According to the melting curve for eimeria 
lab <- lab %>%
    dplyr::mutate(current_infection = case_when(
        Parasite_challenge == "E_ferrisi" & MC.Eimeria == "TRUE" ~ "E_ferrisi",
        Parasite_challenge == "E_ferrisi" & MC.Eimeria == "FALSE" ~ "uninfected",
        Parasite_challenge == "E_falciformis" & MC.Eimeria == "TRUE" ~ "E_falciformis",
        Parasite_challenge == "E_falciformis" & MC.Eimeria == "FALSE" ~ "uninfected",
        Parasite_challenge == "uninfected" & MC.Eimeria == "TRUE" ~ "infected_eimeria",
        Parasite_challenge == "uninfected" & MC.Eimeria == "FALSE" ~ "uninfected",
        TRUE ~ ""
    ))

# current falciformis
lab <- lab %>%
    dplyr::mutate(infection = case_when(
        Parasite_challenge == "E_ferrisi" & MC.Eimeria == "TRUE" ~ "E_ferrisi",
        Parasite_challenge == "E_ferrisi" & MC.Eimeria == "FALSE" ~ "uninfected",
        Parasite_challenge == "E_falciformis" & MC.Eimeria == "TRUE" ~ "E_falciformis",
        Parasite_challenge == "E_falciformis" & MC.Eimeria == "FALSE" ~ "uninfected",
        Parasite_challenge == "uninfected" & MC.Eimeria == "TRUE" ~ "E_falciformis",
        Parasite_challenge == "uninfected" & MC.Eimeria == "FALSE" ~ "uninfected",
        TRUE ~ ""
    ))

# PCA
## we can now run a normal pca on the complete data set
res.pca <- PCA(genes)


## How much does each dimension contribute to variance?
    
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 70))


fviz_pca_var(res.pca, col.var = "cos2",
             gradient.cols = c("#DB6212", "#CC8733", "#5f25e6", "#073DA8"),
             repel = TRUE, title = "")

fviz_pca_ind(res.pca, col.ind = "cos2", 
                  gradient.cols = c("#DB6212", "#CC8733", "#5f25e6", "#073DA8"), 
                  repel = TRUE, title = "")

## Description of the dimensions
## We get a correlation between each variable and the first dimension
dimdesc(res.pca)


## Adding the PC Eigenvectors to the data set. 

mouse_id <- gene %>%
  dplyr::select(Mouse_ID)

mouse_id$pc1 <- res.pca$ind$coord[, 1] # indexing the first column

mouse_id$pc2 <- res.pca$ind$coord[, 2]  # indexing the second column

lab <- lab %>% 
  left_join(mouse_id, by = "Mouse_ID")



## We also need to extract the data for the variable contributions to each of 
## the pc axes.
pca.vars <- res.pca$var$coord %>% data.frame

pca.vars$vars <- rownames(pca.vars)

pca.vars.m <- melt(pca.vars, id.vars = "vars")

source("r_scripts/functions/circle_fun.R")

circ <- circleFun(c(0,0),2,npoints = 500)


#It’s possible to use the function corrplot() [corrplot package] to highlight 
#the most contributing variables for each dimension:
var.contrib <- as.data.frame(res.pca$var$contrib)
corrplot(var.contrib, is.corr=FALSE) 


pca_var <- as.data.frame(pca.vars)

## save the variance contribution of each gene 
##save the normalized data 
write.csv(pca_var, "output_data/variance_contr_gene_lab", row.names = TRUE)

### Contributions to the first dimension

# Contributions of variables to PC1
fviz_contrib(res.pca, choice = "var", axes = 1, top = 18, 
             title = "Contribution of immune genes to the first dimension of the PCA")

# res.pca$var$contrib


### Contributions to the second dimension

## Contributions of variables to PC2
fviz_contrib(res.pca, choice = "var", axes = 2, top = 18, 
             title = "Contribution of immune genes to the second dimension of the PCA")

fviz_contrib(res.pca, choice = "var", axes = 1:2, top = 18)

# Total contribution on PC1 and PC2
fviz_contrib(res.pca, choice = "ind", axes = 1:2)

#select same rows in the first table
lab <- lab[row.names(genes), ]

fviz_pca_biplot(res.pca, 
                col.ind = lab$infection, palette = "jco", 
                addEllipses = TRUE, label = "var",
                col.var = "black", repel = TRUE,
                legend.title = "Infection groups",
                title = "") 




################## Linear models: Predicting weight loss with the PCA eigenvectors

# predicting weight loss with the pc1 and pc2
model_1_pc1_pc2 <- lm(WL_max ~ pc1 + pc2, data = lab)
summary(model_1_pc1_pc2)
AIC(model_1_pc1_pc2)

### use the ggefects package
# 

# predicting weight loss with pc1
model_1_pc1 <- lm(WL_max ~ pc1, data = lab)
summary(model_1_pc1)
AIC(model_1_pc1)


# Here the base of comparison is E_ferrisi. I should change it to compare to the
# uninfected mice
model_2_pc1_pc2_challenge <- lm(WL_max ~ pc1 + pc2 + infection, data = lab)
summary(model_2_pc1_pc2_challenge)
AIC(model_2_pc1_pc2_challenge)


# homozygous / heterozygous infection


# correlations between infection and the immune responses
plot(model_2_pc1_pc2_challenge)

# covariance matrix of the fixed effects of the model


model_3_infection_hybrid_status <- lm(WL_max ~ pc1 + pc2 + hybrid_status, 
                 data = lab)
summary(model_3_infection_hybrid_status)
AIC(model_3_infection_hybrid_status)

### information of previous infections, some mice are actually reacting for the 
# first time, as the



# Compare
llr_test <- anova(model_1_pc1_pc2, model_2_pc1_pc2_challenge)
print(llr_test)

# model_2_pc1_pc2_challenge <- lm(WL_max ~ pc1 + pc2 + infection, data = lab)

weight_no_pc1 <- lm(WL_max ~ pc2 + infection, data = lab)
weight_no_pc2 <- lm(WL_max ~ pc1  + infection, data = lab)
weight_no_infection <- lm(WL_max ~ pc1 + pc2, data = lab)
lrtest(model_2_pc1_pc2_challenge, weight_no_pc1)
lrtest(model_2_pc1_pc2_challenge, weight_no_pc2)
lrtest(model_2_pc1_pc2_challenge, weight_no_infection)
lrtest(weight_no_pc1, weight_no_pc2)


## Visualizing the regression models
### scatter plot with regression lines 

# for the model with pc1 and pc2
ggplot(lab, aes(x = pc1, y = WL_max)) +
  geom_point(aes(color = infection)) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  labs(x = "PC1", y = "Maximum Weight Loss") +
  theme_bw()


# for the model with pc2 and infection
ggplot(lab, aes(x = pc2, y = WL_max)) +
  geom_point(aes(color = infection)) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  labs(x = "PC2", y = "Maximum Weight Loss") +
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



# Heatmap
### repeating the heatmap on the now imputed data

 # turn the data frame into a matrix and transpose it. We want to have each cell 
 # type as a row name 
 gene <- t(as.matrix(gene))
 
 # turn the first row into column names
 gene %>%
     row_to_names(row_number = 1) -> heatmap_data
 
 heatmap_data <- as.data.frame(heatmap_data)
 
 table(rowSums(is.na(heatmap_data)) == nrow(heatmap_data))

 
# turn the columns to numeric other wise the heatmap function will not work
 heatmap_data[] <- lapply(heatmap_data, function(x) as.numeric(as.character(x)))

 # remove columns with only NAs 
 heatmap_data <- Filter(function(x)!all(is.na(x)), heatmap_data) 
 
 #remove rows with only Nas
 heatmap_data <-  heatmap_data[, colSums(is.na(heatmap_data)) != 
                                   nrow(heatmap_data)]
 
  
#Prepare the annotation data frame
annotation_df <- as_tibble(lab) %>%
    dplyr::select(c("Mouse_ID",  "WL_max", "infection")) 
  
annotation_df <- unique(annotation_df) 

annotation_df <- as.data.frame(annotation_df)




### Prepare the annotation columns for the heatmap
rownames(annotation_df) <- annotation_df$Mouse_ID


# Match the row names to the heatmap data frame
rownames(annotation_df) <- colnames(heatmap_data)

#remove the unecessary column
annotation_df <- annotation_df %>% dplyr::select(-Mouse_ID, )


# Define colors for each parasite
parasite_colors <- c("E_falciformis" = "coral2",
                     "E_ferrisi" = "chartreuse4",
                     "uninfected" = "cornflowerblue")

# Generate the heat map
pheatmap(heatmap_data, annotation_col = annotation_df,  scale = "row",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         annotation_colors = list(infection = parasite_colors)) # use annotation_colors


####### Gene enrichment analysis
#create a new vector n to match the genes with gene ids in the package
# Enrichr

# Get the available keytypes in the database
keytypes <- keytypes(org.Mmu.eg.db)

# View the list of keytypes
print(keytypes)

# creating vector according to gene names NIH
#https://www.ncbi.nlm.nih.gov/gene/?term=Mus+musculus+Prf1
gene_ids <- c("15978", #IFNG, 
               "12766", #CXCR3"
               "16193", #IL6
               "16163", #IL13",
               "16181", #IL1RN
               "12362", #CASP1
               "17329", #CXCL9
               "15930", #IDO1
               "15944", #IRGM1
               "17523", #MPO"
               "17831", #MUC2
               "17833", #MUC5AC
               "17874", #MYD88
               "17086", #NCR1
               "18646", #PRF1
               "57263", #RETNLB
               "12703", #SOCS1
               "106759", #TICAM1
               "21926") #TNF


gene_symbols <- c("IFNG", "CXCR3", "IL6", "IL13", "IL1RN", "CASP1", "CXCL9", "IDO1",
              "IRGM1", "MPO", "MUC2", "MUC5AC", "MYD88", "NCR1", "PRF1", "RETNLB",
              "SOCS1", "TICAM1", "TNF")



# Perform gene ontology enrichment analysis
enrich_result <- enrichGO(gene = gene_symbols,
                          OrgDb = org.Mmu.eg.db,
                          ont = "BP",
                          keyType = "SYMBOL")


# View the enrichment result
enrich_result

           
# Extract the relevant columns from the enrichment table
enriched_terms <- enrich_result$Description
p_values <- enrich_result$p.adjust
gene_ratio <- enrich_result$GeneRatio

# Sort the enriched terms based on p-values
sorted_terms <- enriched_terms[order(p_values)]

# Create the bar plot
barplot(-log10(p_values), names.arg = sorted_terms, horiz = FALSE,
        xlab = "Enriched GO Terms", ylab = "-log10(p-value)",
        main = "Gene Ontology Enrichment Analysis",
        col = "steelblue", border = "black")

## Now go on to select the interest groupings seen on the pca
# IL.13

# TICAM1

# NCR1, SOCS1, IRGM1, MUC2
gene_symbols <- c("NCR1", "SOCS1", "IRGM1", "MUC2")

# Perform gene ontology enrichment analysis
enrich_result <- enrichGO(gene = gene_symbols,
                          OrgDb = org.Mmu.eg.db,
                          ont = "BP",
                          keyType = "SYMBOL")

# View the enrichment result
enrich_result


# Extract the relevant columns from the enrichment table
enriched_terms <- enrich_result$Description
p_values <- enrich_result$p.adjust
gene_ratio <- enrich_result$GeneRatio

# Sort the enriched terms based on p-values
sorted_terms <- enriched_terms[order(p_values)]

# Create the bar plot
barplot(-log10(p_values), names.arg = sorted_terms, horiz = FALSE,
        xlab = "Enriched GO Terms", ylab = "-log10(p-value)",
        main = "Gene Ontology Enrichment Analysis",
        col = "steelblue", border = "black")

# Extract the enriched terms and gene ratios
enriched_terms <- enrich_result$Description
gene_ratio <- enrich_result$GeneRatio
p_value <- enrich_result$pvalue

# Create a data frame for the heatmap data
data_df <- data.frame(Gene_Ratio = as.numeric(as.factor(gene_ratio[complete_cases])),
                      P_Value = as.numeric(p_value[complete_cases]),
                      stringsAsFactors = FALSE)

rownames(data_df) <- enriched_terms

# Convert the data frame to a matrix
data_matrix <- as.matrix(data_df)

# Sort the data frame by P_Value in ascending order
sorted_df <- data_df[order(data_df$P_Value), ]

# Select the first 15 rows
sorted_df[1:15, ]

#negative regulation of alpha-beta T cell differentiation                        
#positive regulation of CD4-positive, alpha-beta T cell differentiation


# CASP1, PRF1, CXCR3, IL6, MUC5AC

# ILRN

# MPO, IFNG, CXCL9, TNF



# save the lab data frame for figures
write.csv(lab, "output_data/lab_pca", row.names = FALSE)
