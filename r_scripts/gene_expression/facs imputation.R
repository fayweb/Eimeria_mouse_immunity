################
###############
#######
########
### repeating the heatmap on the now imputed data
```{r}

### Prepare the annotation data frame for the heatmap
#annotation_df <-  %>%
#  dplyr::select(c("EH_ID", "Parasite_challenge", "infection_history"))

### Data tidying for the heatmap function
gene_imp <- g %>% left_join(imputed_gene, by = c("pc1", "pc2"))

#remove all columns of the non-imputed data
gene_imp = gene_imp[,!grepl(".x$",names(gene_imp))]

#remove the suffix y
gene_imp <- gene_imp %>% rename_with(~str_remove(., '.y'))

gene <- gene_imp %>% dplyr::select(c("EH_ID", "CXCR3_bio", "IL.6",
                                     "IL.10", "IL.13", "IL.10", "IL.13", "IL1RN", "CASP1", "CXCL9", 
                                     "IDO1", "IRGM1", "MPO", "MUC2", "MUC5AC", "MYD88", 
                                     "NCR1", "PRF1", "RETNLB", "SOCS1", "TICAM1", "TNF"))

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


table(rowSums(is.na(heatmap_data)) == nrow(heatmap_data))


# turn the columns to numeric other wise the heatmap function will not work
heatmap_data[] <- lapply(heatmap_data, function(x) as.numeric(as.character(x)))

# remove columns with only NAs 
heatmap_data <- Filter(function(x)!all(is.na(x)), heatmap_data) 

#remove rows with only Nas
heatmap_data <-  heatmap_data[, colSums(is.na(heatmap_data)) != nrow(heatmap_data)]

rownames(annotation_df) <- colnames(heatmap_data)
```
Heatmap on gene expression data: 
  ```{r, echo = FALSE}
pheatmap(heatmap_data, annotation_col = annotation_df, scale = "row")
```
