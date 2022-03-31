library(tidyverse)
library(tidyr)
library(janitor)
library(tibble)


### Import the data
Challenge <- read.csv("https://raw.githubusercontent.com/fayweb/Eimeria_mouse_immunity/main/data/01_raw_files/Challenge_infections.csv")

## vectors for selecting columns
basics <- c("EH_ID", "mouse_strain", "experiment", "primary_infection", 
            "challenge_infection", "labels", "dpi", "infection", 
            "infection_history","death")

weight_loss <- c("weight", "weight_dpi0", "relative_weight")

CellCount.cols <- c("Position", "CD4", "Treg", "Div_Treg", "Treg17", "Th1", 
                    "Div_Th1", "Th17", "Div_Th17", "CD8", "Act_CD8", 
                    "Div_Act_CD8", "IFNy_CD4", "IFNy_CD8","Treg_prop", 
                    "IL17A_CD4", "batch")  

### Select the necessary columns (as shown in the vectors above)
df_facs <- Challenge %>% select(basics, weight_loss, CellCount.cols)

### Select the measurements from the mesenterial lymphnodes
# The spleen measurements are not so plenty and not comparable to the mln 
# measurements
df_facs_mln <- df_facs %>% filter(Position == "mLN")

### Drop the columns that contain nas in some of the facs columns
df_facs_mln <- df_facs_mln %>% drop_na("CD4")

### See the structure of the data frame
glimpse(df_facs_mln) #Cell count columns values are indeed numeric

write.csv(df_facs_mln, "data_products/01_intermediate_files/01_lab_facs_mLN", row.names=FALSE)



