library(pheatmap)
library(tidyverse)
library(matrixStats)
library(tidyr)
library(janitor)
library(tibble)

### Import the data 

facs_lab <- read.csv("https://raw.githubusercontent.com/fayweb/Eimeria_mouse_immunity/main/data_products/01_intermediate_files/01_lab_facs_mLN")


### Prepare the annotation data frame 

annotation_df <- fac %>% select(c("EH_ID", "primary_infection", 
                                     "challenge_infection", "infection_history",
                                     "mouse_strain", "relative_weight"))
