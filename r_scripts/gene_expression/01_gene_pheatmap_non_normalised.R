library(pheatmap)
library(tidyverse)
library(matrixStats)
library(tidyr)
library(janitor)
library(tibble)


### Import the data
Challenge <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data_products/Challenge_infections.csv")

Genes <- c("IFNy", "CXCR3_bio", "IL.6",
  "IL.10", "IL.13", "IL.10", "IL.13", "IL1RN", "CASP1", "CXCL9", 
  "IDO1", "IRGM1", "MPO", "MUC2", "MUC5AC", "MYD88", 
  "NCR1", "PRF1", "RETNLB", "SOCS1", "TICAM1", "TNF")

basics_gene <- as_tibble(Challenge) %>%
  dplyr::filter(infection == "challenge", dpi == "8") %>%
  dplyr::group_by(EH_ID, infection) %>%
  dplyr::select(c("EH_ID", "primary_infection", "challenge_infection", "infection_history",
           "mouse_strain", "max_WL", Genes, "delta"))

Challenge <- Challenge %>%
    dplyr::mutate(Parasite_primary = case_when(
        primary_infection == "E64" ~ "Eimeria ferrisi",
        primary_infection == "E88" ~ "Eimeria falciformis",
        primary_infection == "Eflab" ~ "Eimeria falciformis",
        primary_infection == "E139" ~ "Eimeria ferrisi",
        primary_infection == "UNI" ~ "uninfected",
        TRUE ~ ""))


Challenge <- Challenge %>%
    dplyr::mutate(Parasite_challenge = case_when(    
        challenge_infection == "E64" ~ "Eimeria ferrisi",
        challenge_infection == "E88" ~ "Eimeria falciformis",
        challenge_infection == "Eflab" ~ "E. falciformis",
        challenge_infection == "E139" ~ "Eimeria ferrisi",
        challenge_infection == "UNI" ~ "uninfected",
        TRUE ~ ""))

# remove duplicated and change the format to data frame 
basics_gene <- unique(basics_gene) %>% as.data.frame(basics_gene) %>% select(-infection)
 
gene <- basics_gene %>% dplyr::select(c(EH_ID, Genes))

 
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
 

 
 # How do our data look like?
 glimpse(heatmap_data) #alll columns appear to be numeric
 
 # turn the columns to numeric other wise the heatmap function will not work
 heatmap_data[] <- lapply(heatmap_data, function(x) as.numeric(as.character(x)))

 # remove columns with only NAs 
 heatmap_data <- Filter(function(x)!all(is.na(x)), heatmap_data) 
 
 #remove rows with only Nas
 heatmap_data <-  heatmap_data[, colSums(is.na(heatmap_data)) != nrow(heatmap_data)]

 rowSums(is.na(heatmap_data))
 

 
### Prepare the annotation data frame for the heatmap
 
## First select the common mice between your two data frames
 #now we have different row names in the two data frames because I removed the Nas
 #join the heatmap data to the basic data frame and get only the common columns
 gene_na_omit <- basics_gene %>% 
   select(c(EH_ID, primary_infection, challenge_infection, infection_history, mouse_strain, 
            max_WL, delta)) %>%
     inner_join((t(heatmap_data) %>% as.data.frame() %>% tibble::rownames_to_column("EH_ID")), 
                by = "EH_ID")
 
 gene_na_omit <- gene_na_omit %>%
     dplyr::mutate(Parasite_challenge = case_when(    
         challenge_infection == "E64" ~ "Eimeria ferrisi",
         challenge_infection == "E88" ~ "Eimeria falciformis",
         challenge_infection == "Eflab" ~ "E. falciformis",
         challenge_infection == "E139" ~ "Eimeria ferrisi",
         challenge_infection == "UNI" ~ "uninfected",
         TRUE ~ ""))
 
annotation_df <- gene_na_omit %>%
  select(c("EH_ID", "Parasite_challenge", "infection_history", "max_WL"))


annotation_df <- as.data.frame(annotation_df)

### Prepare the annotation columns for the heatmap
rownames(annotation_df) <- annotation_df$EH_ID



# Match the row names to the heatmap data frame
rownames(annotation_df) <- colnames(heatmap_data)

#remove the unecessary column
annotation_df <- annotation_df %>% select(-EH_ID, )




#plot the heatmap
jpeg("output_data/gene_expression/01_Pheatmap_gene_lab_non_normalised.jpg", width = 1400, height = 1000)

pheatmap(heatmap_data, annotation_col = annotation_df, scale = "row")

#close the jpeg file
dev.off()

pdf("output_data/gene_expression/01_Pheatmap_gene_lab_non_normalised.pdf", width = 14, height = 10)

pheatmap(heatmap_data, annotation_col = annotation_df, scale = "row")


#close the jpeg file
dev.off()


#switch off all dev devices
while (!is.null(dev.list()))  dev.off()
