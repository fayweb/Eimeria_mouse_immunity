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

### Create a new column showing maximum weight loss over every dpi for each
## mouse
# select the necessary columns
Cha_dpi_chal <- Challenge %>% filter(infection == "challenge") %>%
  select(c(EH_ID, dpi, weight, weight_dpi0, death))

Cha_dpi_prim  <- Challenge %>% filter(death == "prim_11") %>%
  select(c(EH_ID, dpi, weight, weight_dpi0, death))

Cha_dpi <- rbind(Cha_dpi_chal, Cha_dpi_prim)

# weight: the weight of the mouse at this dpi
# weight_dpi0: the weight at the day of infection

# create a new column showing weight loss
Cha_dpi <- Cha_dpi %>%
  mutate(weight_change = ((weight_dpi0 - weight)/ weight_dpi0) * 100 ) %>%
  #remove now the columns weight_dpi0 and weight
  select(-c(weight_dpi0, weight))

Cha_dpi <- unique(Cha_dpi)

# change the format to wide
Cha_dpi_wide <- Cha_dpi %>% pivot_wider(names_from = dpi, values_from = weight_change)

Cha_dpi_wide$max_weight_loss <- apply(Cha_dpi_wide[2:8], MARGIN =  1, FUN = min, na.rm = FALSE)

### Select the measurements from the mesenterial lymphnodes
# The spleen measurements are not so plenty and not comparable to the mln 
# measurements
df_facs_mln <- df_facs %>% filter(Position == "mLN")


### Drop the columns that contain nas in some of the facs columns
df_facs_mln <- df_facs_mln %>% drop_na("CD4")

### See the structure of the data frame
glimpse(df_facs_mln) #Cell count columns values are indeed numeric

### join the maximum weight loss
Cha_weight_loss <- Cha_dpi_wide %>% select(EH_ID, max_weight_loss)

df_facs_mln <- df_facs_mln %>% left_join(Cha_weight_loss, by = "EH_ID")

# for some mice the weight was not noted on day 0
#how should I evaluate this?


#make github work
write.csv(df_facs_mln, "data_products/01_intermediate_files/01_lab_facs_mLN", row.names=FALSE)



