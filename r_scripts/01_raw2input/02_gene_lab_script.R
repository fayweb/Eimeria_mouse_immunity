library(tidyverse)
library(tidyr)
library(janitor)
library(tibble)


### Import the data
Challenge <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data_products/Challenge_infections.csv")

## vectors for selecting columns
basics <- c("EH_ID", "mouse_strain", "experiment", "primary_infection", 
            "challenge_infection", "labels", "dpi", "infection", 
            "infection_history","death")

weight_loss <- c("weight", "weight_dpi0", "relative_weight")

Gene.Exp.cols   <- c("IFNy", "IL.12", "IRG6", "CXCR3", "IL.6", "IL.10", 
                     "IL.13", "IL.10", "IL.13", "IL1RN","CXCR3", "CASP1", "CXCL9", 
                     "IDO1", "IRGM1", "MPO", "MUC2", "MUC5AC", "MYD88", 
                     "NCR1", "PRF1", "RETNLB", "SOCS1", "TICAM1", "TNF")

#somehow gene "GBP2" doesn't exist in challenge, hence removed from vector
### Select the necessary columns (as shown in the vectors above)
df_gene <- Challenge %>% select(basics, weight_loss, Gene.Exp.cols)

### Create a new column showing maximum weight loss over every dpi for each
## mouse
# select the necessary columns
Cha_dpi_chal <- Challenge %>% filter(infection == "challenge") %>%
  select(c(EH_ID, dpi, weight, weight_dpi0, death))

Cha_dpi_prim  <- Challenge %>% filter(death == "prim_11") %>%
  select(c(EH_ID, dpi, weight, weight_dpi0, death))

Cha_dpi <- rbind(Cha_dpi_chal, Cha_dpi_prim)
rm(Cha_dpi_chal)
rm(Cha_dpi_prim)

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
rm(Cha_dpi)

Cha_dpi_wide$max_weight_loss <- apply(Cha_dpi_wide[2:8], MARGIN =  1, FUN = min, na.rm = FALSE)

### Select the measurements from the mesenterial lymphnodes
# The spleen measurements are not so plenty and not comparable to the mln 
# measurements
#df_facs_mln <- df_facs %>% filter(Position == "mLN")


### Drop the columns that contain nas in some of the facs columns
df_facs_mln <- df_facs_mln %>% drop_na("CD4")

### See the structure of the data frame
glimpse(df_facs_mln) #Cell count columns values are indeed numeric

### join the maximum weight loss
Cha_weight_loss <- Cha_dpi_wide %>% select(EH_ID, max_weight_loss)

df_facs_mln <- df_facs_mln %>% left_join(Cha_weight_loss, by = "EH_ID")

# for some mice the weight was not noted on day 0
#how should I evaluate this?

write.csv(df_facs_mln, "data_products/01_intermediate_files/01_lab_facs_mLN", row.names=FALSE)
