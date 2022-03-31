library(pheatmap)
library(tidyverse)
library(matrixStats)
library(tidyr)
library(janitor)
library(tibble)

### Import the data 

### Import the data 
facs_lab <- read.csv("https://raw.githubusercontent.com/fayweb/Eimeria_mouse_immunity/main/data_products/01_intermediate_files/01_lab_facs_mLN")

# Here is a vector to select the cell count columns
CellCount.cols <- c("CD4", "Treg", "Div_Treg", "Treg17", "Th1", "Div_Th1", "Th17", 
                    "Div_Th17", "CD8", "Act_CD8", "Div_Act_CD8", "IFNy_CD4", 
                    "IFNy_CD8","Treg_prop", "IL17A_CD4")

### Prepare the annotation data frame 

### Prepare the annotation data frame for the heatmap
annotation_df <- facs_lab %>% 
  select(c("EH_ID", "primary_infection", "challenge_infection", "infection_history",
           "mouse_strain", "max_weight_loss"))

### Data tidying for the heatmap function
facs <- facs_lab  %>% 
  select(c(EH_ID, CellCount.cols))

# turn the data frame into a matrix and transpose it. We want to have each cell 
# type as a row name 
facs <- t(as.matrix(facs))

#switch the matrix back to a data frame format
facs <- as.data.frame(facs)

 # turn the first row into column names
facs %>%
  row_to_names(row_number = 1) -> facs

# Now further prepare the data frame for plotting by removing the first row
 ## and convert the column to row names with the cells 
facs[-1, ] -> heatmap_data

 # How do our data look like?
  glimpse(heatmap_data) #alll columns appear to be numeric

# turn the columns to numeric other wise the heatmap function will not work
heatmap_data[] <- lapply(heatmap_data, function(x) as.numeric(as.character(x)))

### Prepare the annotation columns for the heatmap
rownames(annotation_df) <- annotation_df$EH_ID

# Match the row names to the heatmap data frame
rownames(annotation_df) <- colnames(facs)

#remove the unecessary column
annotation_df <- annotation_df %>% select(-EH_ID, )


#plot the heatmap
jpeg("data_products/02_output_files/02_plots/Pheatmap_facs_lab.jpg", width = 800, height = 1000)
heatmap_data %>% 
   pheatmap(annotation_col = annotation_df, scale = "row")
#close the jpeg file
dev.off()
