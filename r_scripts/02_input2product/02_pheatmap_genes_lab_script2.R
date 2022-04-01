library(pheatmap)
library(tidyverse)
library(matrixStats)
library(tidyr)
library(janitor)
library(tibble)

source("r_scripts/01_raw2input/02_gene_lab_script.R")

### Import the data 

### Import the data 
gene_lab <- df_gene

## vectors for selecting columns

### Prepare the annotation data frame 

### Prepare the annotation data frame for the heatmap
annotation_df <- gene_lab %>% 
  select(c("EH_ID", "primary_infection", "challenge_infection"))

#missing from the annotation columns 
#"max_weight_loss"
### Data tidying for the heatmap function
gene <- gene_lab  %>% 
  select(c(EH_ID, Gene.Exp.cols.norm ))

# turn the data frame into a matrix and transpose it. We want to have each cell 
# type as a row name 
gene <- t(as.matrix(gene))

#switch the matrix back to a data frame format
gene <- as.data.frame(gene)

# turn the first row into column names
gene %>%
  row_to_names(row_number = 1) -> gene

# Now further prepare the data frame for plotting by removing the first row
## and convert the column to row names with the cells 
gene[-1, ] -> heatmap_data


table(rowSums(is.na(heatmap_data))==nrow(heatmap_data))

#hhh
# How do our data look like?
glimpse(heatmap_data) #alll columns appear to be numeric

# turn the columns to numeric other wise the heatmap function will not work
heatmap_data[] <- lapply(heatmap_data, function(x) as.numeric(as.character(x)))

annotation_df <- as.data.frame(annotation_df)

### Prepare the annotation columns for the heatmap
rownames(annotation_df) <- annotation_df$EH_ID

# Match the row names to the heatmap data frame
rownames(annotation_df) <- colnames(gene)

#remove the unecessary column
annotation_df <- annotation_df %>% select(-EH_ID, )

heatmap_data2 <- heatmap_data[rowSums(is.na(heatmap_data))<41,
                              colSums(is.na(heatmap_data))<6]


#plot the heatmap
jpeg("data_products/02_output_files/02_plots/Pheatmap_gene_lab.jpg", width = 800, height = 1000)

pheatmap(heatmap_data2, annotation_col = annotation_df, scale = "row")


#close the jpeg file
dev.off()

