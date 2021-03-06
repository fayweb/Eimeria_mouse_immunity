library(tidyverse)
library(tidyr)
library(janitor)
library(tibble)
library(pheatmap)
library(matrixStats)
library(tibble)

### Import the data
Challenge <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data_products/Challenge_infections.csv")

## vectors for selecting columns
CellCount.cols <- c("Position", "CD4", "Treg", "Div_Treg", "Treg17", "Th1", 
                    "Div_Th1", "Th17", "Div_Th17", "CD8", "Act_CD8", 
                    "Div_Act_CD8", "IFNy_CD4", "IFNy_CD8","Treg_prop", 
                    "IL17A_CD4")  

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

### Select the necessary columns (as shown in the vectors above)
#df_facs <- Challenge %>% select(basics, weight_loss, CellCount.cols)



### Select the measurements from the mesenterial lymphnodes
# The spleen measurements are not so plenty and not comparable to the mln 
# measurements
Challenge <- Challenge %>% filter(Position == "mLN")


### Drop the columns that contain nas in the column CD4 of the facs columns
FACS <- Challenge %>% drop_na("CD4")

### Prepare the annotation data frame for the heatmap
annotation_df <- FACS %>% 
  filter(infection == "challenge") %>%
  select(c("EH_ID", "Parasite_challenge", "infection_history"))

### Data tidying for the heatmap function
FACS <- FACS  %>% 
  select(c(EH_ID, CellCount.cols))

# turn the data frame into a matrix and transpose it. We want to have each cell 
# type as a row name 
FACS <- t(as.matrix(FACS))

#switch the matrix back to a data frame format
FACS <- as.data.frame(FACS)

# turn the first row into column names
FACS %>%
  row_to_names(row_number = 1) -> FACS

# Now further prepare the data frame for plotting by removing the first row
## and convert the column to row names with the cells 
FACS[-1, ] -> heatmap_data

# How do our data look like?
glimpse(heatmap_data) #alll columns appear to be characters

# turn the columns to numeric other wise the heatmap function will not work
heatmap_data[] <- lapply(heatmap_data, function(x) as.numeric(as.character(x)))

### Prepare the annotation columns for the heatmap
rownames(annotation_df) <- annotation_df$EH_ID

# Match the row names to the heatmap data frame
rownames(annotation_df) <- colnames(FACS)

#remove the unecessary column
annotation_df <- annotation_df %>% dplyr::select(-EH_ID, )


#plot the heatmap


pdf("Eimeria_mouse_immunity/output_data/facs/01_Pheatmap_facs_labr.pdf", width = 14, height = 10)
heatmap_data %>% 
  pheatmap(annotation_col = annotation_df, scale = "row")
#close the jpeg file
dev.off()

while (!is.null(dev.list()))  dev.off()
