library(tidyverse)
library(tidyr)
library(janitor)
library(tibble)

### Go back into the data prep (producing Challenge.csv)
# Give the columns names that allow selection:
# E.g. GE_xyz for gene-expression NGE_xyz for normalized gene-expression (xyz is
# a gene-name), FACS_xyz etc...

### clean the naming conventions in the analysis repository! 

## use summarize data per mouse at the beginning of the analysis here (e.g. 
# using the example code in the repo or better mutate instead of summarize
# to keep the other values)

## get one dataset to work with, try to not litter the environment

## in addition to pheatmap:
# plot gene expression against intensity and ESpecies in challenge infection
# (different genes e.g. in different panels/facets)

## cor/cor.test for each gene expression with intensity, 
## even lm (linear models?)

## generally you want the data in a state at which you can do any of this
## in as few lines of code as possible

### Import the data
Challenge <- read.csv("https://raw.githubusercontent.com/derele/Eimeria_Lab/master/data_products/Challenge_infections.csv")

cor(Challenge[, grepl("_N", colnames(Challenge))], use="pairwise.complete.obs")

## vectors for selecting columns
basics <- c("EH_ID", "mouse_strain", "experiment", "primary_infection", 
            "challenge_infection", "labels", "dpi", "infection", 
            "infection_history","death")

weight_loss <- c("weight", "weight_dpi0", "relative_weight")

Gene.Exp.cols   <- c("IFNy", "IL.12", "IRG6", "CXCR3", "IL.6", "IL.10", 
                     "IL.13", "IL.10", "IL.13", "IL1RN","CXCR3", "CASP1", "CXCL9", 
                     "IDO1", "IRGM1", "MPO", "MUC2", "MUC5AC", "MYD88", 
                     "NCR1", "PRF1", "RETNLB", "SOCS1", "TICAM1", "TNF")

Gene.Exp.cols.norm   <-paste0(Gene.Exp.cols, "_N")

Gene.Exp.cols.norm <- Gene.Exp.cols.norm[!Gene.Exp.cols.norm%in%c("IL.12_N", "IRG6_N")]

## this is only a figure, so ideally we would not produce any objects but just
## create the figure
##Challenge %>% ### just start with the data
##  filter(dpi==death) ## select only the dissection dpi 

### this falls into the same category 



##  select(ends_with("_N")) %>% ### select the columns you want to plot
##  filter(if_any(everything(), ~ !is.na(.)))%>% ## remove rows only NA
##  summarise_each(funs(sum(is.na(.)))) ### remove columns only NA
##.... 
##... 
##pheatmap() ## do the pheatmap

#somehow gene "GBP2" doesn't exist in challenge, hence removed from vector
### Select the necessary columns (as shown in the vectors above)
df_gene <- Challenge %>% select(basics, weight_loss, Gene.Exp.cols, Gene.Exp.cols.norm)

# where are the columns IL.

### Create a new column showing maximum weight loss over every dpi for each
## mouse
# select the necessary columns
Cha_dpi_chal <- Challenge %>% filter(infection == "challenge") %>%
  select(c(EH_ID, dpi, weight, weight_dpi0, death, infection, experiment))

Cha_dpi_prim  <- Challenge %>% filter(death == "prim_11") %>%
  select(c(EH_ID, dpi, weight, weight_dpi0, death, infection, experiment))

# combine the two data frames
Cha_dpi <- rbind(Cha_dpi_chal, Cha_dpi_prim)

# remove the unecessary data frames
rm(Cha_dpi_chal)
rm(Cha_dpi_prim)

# weight: the weight of the mouse at this dpi
# weight_dpi0: the weight at the day of infection

# create a new column showing weight loss
Cha_dpi <- Cha_dpi %>%
  mutate(weight_change = ((weight_dpi0 - weight)/ weight_dpi0) * 100 ) %>%
  #remove now the columns weight_dpi0 and weight
  select(-c(weight_dpi0, weight))

# remove duplicates
Cha_dpi <- unique(Cha_dpi)

# change the format to wide
Cha_dpi_wide <- Cha_dpi %>% pivot_wider(names_from = dpi, values_from = weight_change)
rm(Cha_dpi)

Cha_dpi_wide$max_weight_loss <- apply(Cha_dpi_wide[2:8], MARGIN =  1, FUN = min, na.rm = FALSE)

### join the maximum weight loss
Cha_weight_loss <- Cha_dpi_wide %>% select(EH_ID, max_weight_loss, death, 
                                           infection, experiment)
## Add a column for dpi 
Cha_weight_loss1 <- Cha_weight_loss %>% 
  filter(death == "chal_8") %>%
  mutate(dpi = "8")

Cha_weight_loss2 <- Cha_weight_loss %>% 
  filter(death == "prim_11") %>%
  mutate(dpi = "11")

Cha_weight_loss <- rbind(Cha_weight_loss1, Cha_weight_loss2)

Cha_weight_loss$dpi <- as.integer(Cha_weight_loss$dpi)

rm(Cha_weight_loss1)
rm(Cha_weight_loss2)
rm(Cha_dpi_wide)

# now merge the gene expression data
df_gene <- Cha_weight_loss %>% 
  left_join(df_gene, by = intersect(colnames(Cha_weight_loss), colnames(df_gene)))

df_gene <- unique(df_gene)

# for some mice the weight was not noted on day 0
#how should I evaluate this?
