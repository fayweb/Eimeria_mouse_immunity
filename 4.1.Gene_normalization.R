## ----setup, include=FALSE--------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## --------------------------------------------------------------------------------
library(mice)
library(tidyr)
library(tidyverse)
library(VIM)
library(fitdistrplus)
library(fitur)
library(visdat)
library(DESeq2)


## --------------------------------------------------------------------------------
hm <- read.csv("output_data/1.MICE_cleaned_data.csv")


## --------------------------------------------------------------------------------
# Vectors for selecting genes
#Lab genes
# The measurements of IL.12 and IRG6 are done with an other assay and will 
#ignore for now
Gene_lab   <- c("IFNy", "CXCR3", "IL.6", "IL.13", "IL.10",
                "IL1RN","CASP1", "CXCL9", "IDO1", "IRGM1", "MPO", 
                "MUC2", "MUC5AC", "MYD88", "NCR1", "PRF1", "RETNLB", "SOCS1", 
                "TICAM1", "TNF") #"IL.12", "IRG6")

Genes_wild   <- c("IFNy", "CXCR3", "IL.6", "IL.13", "IL.10", 
                  "IL1RN","CASP1", "CXCL9", "IDO1", "IRGM1", "MPO", 
                  "MUC2", "MUC5AC", "MYD88", "NCR1", "PRF1", "RETNLB", "SOCS1", 
                  "TICAM1", "TNF") #, "IL.12", "IRG6")

Facs_lab <- c("CD4", "Treg", "Div_Treg", "Treg17", "Th1", 
                    "Div_Th1", "Th17", "Div_Th17", "CD8", "Act_CD8", 
                    "Div_Act_CD8", "IFNy_CD4", "IFNy_CD8","Treg_prop", 
                    "IL17A_CD4")  

Facs_wild <- c( "Treg", "CD4", "Treg17", "Th1", "Th17", "CD8",
                     "Act_CD8", "IFNy_CD4", "IL17A_CD4", "IFNy_CD8")


## ----imputing_mice---------------------------------------------------------------
hm$Mouse_ID <- str_replace(hm$Mouse_ID, "_", "")

field <- hm %>%
  dplyr::filter(origin == "Field") 

field <- unique(field)

genes_mouse_field <- field %>%
  dplyr::select(c(Mouse_ID, all_of(Genes_wild), GAPDH)) 

genes_field <- genes_mouse_field  %>%
  dplyr::select(-Mouse_ID)
#remove rows with only nas
genes_field <- genes_field[,colSums(is.na(genes_field))<nrow(genes_field)]
#remove colums with only nas 
genes_field <- genes_field[rowSums(is.na(genes_field)) != ncol(genes_field), ]
genes_mouse_field <- genes_mouse_field[row.names(genes_field), ]

##select same rows in the first table
field <- field[row.names(genes_field), ]


###############lab
#select the genes and lab muce
lab <- hm %>%
  dplyr::filter(origin == "Lab", Position == "mLN") #selecting for mln to avoid
# duplicates
lab <- unique(lab)
gene_lab_mouse <- lab %>%
  dplyr::select(c(Mouse_ID, all_of(Gene_lab), PPIB)) 

gene_lab_mouse <- unique(gene_lab_mouse)

genes_lab <- gene_lab_mouse[, -1]

#remove rows with only nas
genes_lab <- genes_lab[,colSums(is.na(genes_lab))<nrow(genes_lab)]

#remove colums with only nas 
genes_lab <- genes_lab[rowSums(is.na(genes_lab)) != ncol(genes_lab), ]

genes_lab <- unique(genes_lab)

#select same rows in the first table
gene_lab_mouse <- gene_lab_mouse[row.names(genes_lab), ]

##select same rows in the first table
lab <- lab[row.names(genes_lab), ]





## --------------------------------------------------------------------------------
#glimpse(hm_selection_g)

#dplyr::select(-Mouse_ID)
# looking at patterns of nas)
#pattern_na <-as.data.frame(md.pattern(field_genes))
sapply(field %>%
         dplyr::select(c(all_of(Genes_wild), "PPIB", "GAPDH")), 
                      function(x) sum(is.na(x)))

sapply(lab %>%
         dplyr::select(c(all_of(Gene_lab), "PPIB", "GAPDH")), 
                      function(x) sum(is.na(x)))
         


## --------------------------------------------------------------------------------


####################### field ##########################
# select first the field samples 

df <- genes_mouse_field

############### IL.13 
# dct
df <- df %>%
  mutate(IL.13_dct = IL.13 - GAPDH)
# mean of dct
dct_mean <- mean(df$IL.13_dct, na.rm = TRUE)
#fold gene expression
 df <- df %>%
   mutate(IL.13_N = 2^ - (IL.13_dct - dct_mean)) %>%
   mutate(IL.13_N = round(IL.13_N, digits = 2))
 
############# I tried writing a nice function but I failed so now I am
 # going t repeat this many many times
 
 ############### IFNy
 # dct
df <- df %>%
  mutate(IFNy_dct = IFNy - GAPDH)

# mean of dct
dct_mean <- mean(df$IFNy_dct, na.rm = TRUE)
 
#fold gene expression
 df <- df %>%
   mutate(IFNy_N = 2^ - (IFNy_dct - dct_mean)) %>%
   mutate(IFNy_N = round(IFNy_N, digits = 2))
 

 ############### CXCR3
 df <- df %>%
  mutate(CXCR3_dct = CXCR3 - GAPDH)

# mean of dct
dct_mean <- mean(df$CXCR3_dct, na.rm = TRUE)
 
#fold gene expression
 df <- df %>%
   mutate(CXCR3_N = 2^ - (CXCR3_dct - dct_mean)) %>%
   mutate(CXCR3_N = round(CXCR3_N, digits = 2))
 
 ############### IL.6
 df <- df %>%
  mutate(IL.6_dct = IL.6 - GAPDH)

# mean of dct
dct_mean <- mean(df$IL.6_dct, na.rm = TRUE)
 
#fold gene expression
 df <- df %>%
   mutate(IL.6_N = 2^ - (IL.6_dct - dct_mean)) %>%
   mutate(IL.6_N = round(IL.6_N, digits = 2))

 ##############  IL1RN
 df <- df %>%
  mutate(IL1RN_dct = IL1RN - GAPDH)

# mean of dct
dct_mean <- mean(df$IL1RN_dct, na.rm = TRUE)
 
#fold gene expression
 df <- df %>%
   mutate(IL1RN_N = 2^ - (IL1RN_dct - dct_mean)) %>%
   mutate(IL1RN_N = round(IL1RN_N, digits = 2))
 
 
 ##############  CASP1
 df <- df %>%
  mutate(CASP1_dct = CASP1 - GAPDH)

# mean of dct
dct_mean <- mean(df$CASP1_dct, na.rm = TRUE)
 
#fold gene expression
 df <- df %>%
   mutate(CASP1_N = 2^ - (CASP1_dct - dct_mean)) %>%
   mutate(CASP1_N = round(CASP1_N, digits = 2))
 
 
 ##############  CXCL9
 df <- df %>%
  mutate(CXCL9_dct = CXCL9 - GAPDH)

# mean of dct
dct_mean <- mean(df$CXCL9_dct, na.rm = TRUE)
 
#fold gene expression
 df <- df %>%
   mutate(CXCL9_N = 2^ - (CXCL9_dct - dct_mean)) %>%
   mutate(CXCL9_N = round(CXCL9_N, digits = 2))
 
  ##############  IDO1
 df <- df %>%
  mutate(IDO1_dct = IDO1 - GAPDH)
# mean of dct
dct_mean <- mean(df$IDO1_dct, na.rm = TRUE)
#fold gene expression
 df <- df %>%
   mutate(IDO1_N = 2^ - (IDO1_dct - dct_mean)) %>%
   mutate(IDO1_N = round(IDO1_N, digits = 2))
 
  ##############  IRGM1
 df <- df %>%
  mutate(IRGM1_dct = IRGM1 - GAPDH)
# mean of dct
dct_mean <- mean(df$IRGM1_dct, na.rm = TRUE)
#fold gene expression
 df <- df %>%
   mutate(IRGM1_N = 2^ - (IRGM1_dct - dct_mean)) %>%
   mutate(IRGM1_N = round(IRGM1_N, digits = 2))
 
  ##############  MPO
 df <- df %>%
  mutate(MPO_dct = MPO - GAPDH)
# mean of dct
dct_mean <- mean(df$MPO_dct, na.rm = TRUE)
#fold gene expression
 df <- df %>%
   mutate(MPO_N = 2^ - (MPO_dct - dct_mean)) %>%
   mutate(MPO_N = round(MPO_N, digits = 2))
 
 
  ##############  MUC2
 df <- df %>%
  mutate(MUC2_dct = MUC2 - GAPDH)
# mean of dct
dct_mean <- mean(df$MUC2_dct, na.rm = TRUE)
#fold gene expression
 df <- df %>%
   mutate(MUC2_N = 2^ - (MUC2_dct - dct_mean)) %>%
   mutate(MUC2_N = round(MUC2_N, digits = 2))
 
  ##############  MUC5AC
 df <- df %>%
  mutate(MUC5AC_dct = MUC5AC - GAPDH)
# mean of dct
dct_mean <- mean(df$MUC5AC_dct, na.rm = TRUE)
#fold gene expression
 df <- df %>%
   mutate(MUC5AC_N = 2^ - (MUC5AC_dct - dct_mean)) %>%
   mutate(MUC5AC_N = round(MUC5AC_N, digits = 2))
 
  ##############  MYD88
 df <- df %>%
  mutate(MYD88_dct = MYD88 - GAPDH)
# mean of dct
dct_mean <- mean(df$MYD88_dct, na.rm = TRUE)
#fold gene expression
 df <- df %>%
   mutate(MYD88_N = 2^ - (MYD88_dct - dct_mean)) %>%
   mutate(MYD88_N = round(MYD88_N, digits = 2))
 
  ##############  NCR1
 df <- df %>%
  mutate(NCR1_dct = NCR1 - GAPDH)
# mean of dct
dct_mean <- mean(df$NCR1_dct, na.rm = TRUE)
#fold gene expression
 df <- df %>%
   mutate(NCR1_N = 2^ - (NCR1_dct - dct_mean)) %>%
   mutate(NCR1_N = round(NCR1_N, digits = 2))
 
 
  ##############  PRF1
 df <- df %>%
  mutate(PRF1_dct = PRF1 - GAPDH)
# mean of dct
dct_mean <- mean(df$PRF1_dct, na.rm = TRUE)
#fold gene expression
 df <- df %>%
   mutate(PRF1_N = 2^ - (PRF1_dct - dct_mean)) %>%
   mutate(PRF1_N = round(PRF1_N, digits = 2))
 
  ##############  RETNLB
 df <- df %>%
  mutate(RETNLB_dct = RETNLB - GAPDH)
# mean of dct
dct_mean <- mean(df$RETNLB_dct, na.rm = TRUE)
#fold gene expression
 df <- df %>%
   mutate(RETNLB_N = 2^ - (RETNLB_dct - dct_mean)) %>%
   mutate(RETNLB_N = round(RETNLB_N, digits = 2))
 
 
 
  ##############  SOCS1
 df <- df %>%
  mutate(SOCS1_dct = SOCS1 - GAPDH)
# mean of dct
dct_mean <- mean(df$SOCS1_dct, na.rm = TRUE)
#fold gene expression
 df <- df %>%
   mutate(SOCS1_N = 2^ - (SOCS1_dct - dct_mean)) %>%
   mutate(SOCS1_N = round(SOCS1_N, digits = 2))
 
  ##############  TICAM1
 df <- df %>%
  mutate(TICAM1_dct = TICAM1 - GAPDH)
# mean of dct
dct_mean <- mean(df$TICAM1_dct, na.rm = TRUE)
#fold gene expression
 df <- df %>%
   mutate(TICAM1_N = 2^ - (TICAM1_dct - dct_mean)) %>%
   mutate(TICAM1_N = round(TICAM1_N, digits = 2))
 
 ##############  TNF
 df <- df %>%
  mutate(TNF_dct = TNF - GAPDH)
# mean of dct
dct_mean <- mean(df$TNF_dct, na.rm = TRUE)
#fold gene expression
 df <- df %>%
   mutate(TNF_N = 2^ - (TNF_dct - dct_mean)) %>%
   mutate(TNF_N = round(TNF_N, digits = 2))
  
 
 df -> df_field



## --------------------------------------------------------------------------------
 ################################# lab 
 # select first the field samples 


df_lab <- gene_lab_mouse

  ############### IFNy
  df_lab <- df_lab %>%
  mutate(IFNy_dct = IFNy - PPIB)
# mean of dct
dct_mean <- mean(df_lab$IFNy_dct, na.rm = TRUE)
#fold gene expression
 df_lab <- df_lab %>%
   mutate(IFNy_N = 2^ - (IFNy_dct - dct_mean)) %>%
   mutate(IFNy_N = round(IFNy_N, digits = 2))
 
 
 
 ############### CXCR3
  df_lab <- df_lab %>%
  mutate(CXCR3_dct = CXCR3 - PPIB)
# mean of dct
dct_mean <- mean(df_lab$CXCR3_dct, na.rm = TRUE)
#fold gene expression
 df_lab <- df_lab %>%
   mutate(CXCR3_N = 2^ - (CXCR3_dct - dct_mean)) %>%
   mutate(CXCR3_N = round(CXCR3_N, digits = 2))
 
 ############### IL.6
   df_lab <- df_lab %>%
  mutate(IL.6_dct = IL.6 - PPIB)
# mean of dct
dct_mean <- mean(df_lab$IL.6_dct, na.rm = TRUE)
#fold gene expression
 df_lab <- df_lab %>%
   mutate(IL.6_N = 2^ - (IL.6_dct - dct_mean)) %>%
   mutate(IL.6_N = round(IL.6_N, digits = 2))
 
 
 ############## IL.13
   df_lab <- df_lab %>%
  mutate(IL.13_dct = IL.13 - PPIB)
# mean of dct
dct_mean <- mean(df_lab$IL.13_dct, na.rm = TRUE)
#fold gene expression
 df_lab <- df_lab %>%
   mutate(IL.13_N = 2^ - (IL.13_dct - dct_mean)) %>%
   mutate(IL.13_N = round(IL.13_N, digits = 2))
 
 ##############  IL1RN
   df_lab <- df_lab %>%
  mutate(IL1RN_dct = IL1RN - PPIB)
# mean of dct
dct_mean <- mean(df_lab$IL1RN_dct, na.rm = TRUE)
#fold gene expression
 df_lab <- df_lab %>%
   mutate(IL1RN_N = 2^ - (IL1RN_dct - dct_mean)) %>%
   mutate(IL1RN_N = round(IL1RN_N, digits = 2))
 
 ##############  CASP1
   df_lab <- df_lab %>%
  mutate(CASP1_dct = CASP1 - PPIB)
# mean of dct
dct_mean <- mean(df_lab$CASP1_dct, na.rm = TRUE)
#fold gene expression
 df_lab <- df_lab %>%
   mutate(CASP1_N = 2^ - (CASP1_dct - dct_mean)) %>%
   mutate(CASP1_N = round(CASP1_N, digits = 2))
 
 
 ##############  CXCL9
   df_lab <- df_lab %>%
  mutate(CXCL9_dct = CXCL9 - PPIB)
# mean of dct
dct_mean <- mean(df_lab$CXCL9_dct, na.rm = TRUE)
#fold gene expression
 df_lab <- df_lab %>%
   mutate(CXCL9_N = 2^ - (CXCL9_dct - dct_mean)) %>%
   mutate(CXCL9_N = round(CXCL9_N, digits = 2))
 
 
  ##############  IDO1
   df_lab <- df_lab %>%
  mutate(IDO1_dct = IDO1 - PPIB)
# mean of dct
dct_mean <- mean(df_lab$IDO1_dct, na.rm = TRUE)
#fold gene expression
 df_lab <- df_lab %>%
   mutate(IDO1_N = 2^ - (IDO1_dct - dct_mean)) %>%
   mutate(IDO1_N = round(IDO1_N, digits = 2))
 
  ##############  IRGM1
   df_lab <- df_lab %>%
  mutate(IRGM1_dct = IRGM1 - PPIB)
# mean of dct
dct_mean <- mean(df_lab$IRGM1_dct, na.rm = TRUE)
#fold gene expression
 df_lab <- df_lab %>%
   mutate(IRGM1_N = 2^ - (IRGM1_dct - dct_mean)) %>%
   mutate(IRGM1_N = round(IRGM1_N, digits = 2))
 
  ##############  MPO
   df_lab <- df_lab %>%
  mutate(MPO_dct = MPO - PPIB)
# mean of dct
dct_mean <- mean(df_lab$MPO_dct, na.rm = TRUE)
#fold gene expression
 df_lab <- df_lab %>%
   mutate(MPO_N = 2^ - (MPO_dct - dct_mean)) %>%
   mutate(MPO_N = round(MPO_N, digits = 2))
 
  ##############  MUC2
   df_lab <- df_lab %>%
  mutate(MUC2_dct = MUC2 - PPIB)
# mean of dct
dct_mean <- mean(df_lab$MUC2_dct, na.rm = TRUE)
#fold gene expression
 df_lab <- df_lab %>%
   mutate(MUC2_N = 2^ - (MUC2_dct - dct_mean)) %>%
   mutate(MUC2_N = round(MUC2_N, digits = 2))
 
  ##############  MUC5AC
   df_lab <- df_lab %>%
  mutate(MUC5AC_dct = MUC5AC - PPIB)
# mean of dct
dct_mean <- mean(df_lab$MUC5AC_dct, na.rm = TRUE)
#fold gene expression
 df_lab <- df_lab %>%
   mutate(MUC5AC_N = 2^ - (MUC5AC_dct - dct_mean)) %>%
   mutate(MUC5AC_N = round(MUC5AC_N, digits = 2))
 
 
  ##############  MYD88
   df_lab <- df_lab %>%
  mutate(MYD88_dct = MYD88 - PPIB)
# mean of dct
dct_mean <- mean(df_lab$MYD88_dct, na.rm = TRUE)
#fold gene expression
 df_lab <- df_lab %>%
   mutate(MYD88_N = 2^ - (MYD88_dct - dct_mean)) %>%
   mutate(MYD88_N = round(MYD88_N, digits = 2))
 
 
  ##############  NCR1
   df_lab <- df_lab %>%
  mutate(NCR1_dct = NCR1 - PPIB)
# mean of dct
dct_mean <- mean(df_lab$NCR1_dct, na.rm = TRUE)
#fold gene expression
 df_lab <- df_lab %>%
   mutate(NCR1_N = 2^ - (NCR1_dct - dct_mean)) %>%
   mutate(NCR1_N = round(NCR1_N, digits = 2))
 
  ##############  PRF1
   df_lab <- df_lab %>%
  mutate(PRF1_dct = PRF1 - PPIB)
# mean of dct
dct_mean <- mean(df_lab$PRF1_dct, na.rm = TRUE)
#fold gene expression
 df_lab <- df_lab %>%
   mutate(PRF1_N = 2^ - (PRF1_dct - dct_mean)) %>%
   mutate(PRF1_N = round(PRF1_N, digits = 2))
 
  ##############  RETNLB
   df_lab <- df_lab %>%
  mutate(RETNLB_dct = RETNLB - PPIB)
# mean of dct
dct_mean <- mean(df_lab$RETNLB_dct, na.rm = TRUE)
#fold gene expression
 df_lab <- df_lab %>%
   mutate(RETNLB_N = 2^ - (RETNLB_dct - dct_mean)) %>%
   mutate(RETNLB_N = round(RETNLB_N, digits = 2))
 
  ##############  SOCS1
   df_lab <- df_lab %>%
  mutate(SOCS1_dct = SOCS1 - PPIB)
# mean of dct
dct_mean <- mean(df_lab$SOCS1_dct, na.rm = TRUE)
#fold gene expression
 df_lab <- df_lab %>%
   mutate(SOCS1_N = 2^ - (SOCS1_dct - dct_mean)) %>%
   mutate(SOCS1_N = round(SOCS1_N, digits = 2))
 
  ##############  TICAM1
   df_lab <- df_lab %>%
  mutate(TICAM1_dct = TICAM1 - PPIB)
# mean of dct
dct_mean <- mean(df_lab$TICAM1_dct, na.rm = TRUE)
#fold gene expression
 df_lab <- df_lab %>%
   mutate(TICAM1_N = 2^ - (TICAM1_dct - dct_mean)) %>%
   mutate(TICAM1_N = round(TICAM1_N, digits = 2))
 
 ##############  TNF
   df_lab <- df_lab %>%
  mutate(TNF_dct = TNF - PPIB)
# mean of dct
dct_mean <- mean(df_lab$TNF_dct, na.rm = TRUE)
#fold gene expression
 df_lab <- df_lab %>%
   mutate(TNF_N = 2^ - (TNF_dct - dct_mean)) %>%
   mutate(TNF_N = round(TNF_N, digits = 2))


## --------------------------------------------------------------------------------
df_lab <- df_lab %>% 
  dplyr::select(-c(all_of(Gene_lab), PPIB, contains("_dct")))

# remove ending _N
df_lab <- df_lab %>%
  rename_with(~str_remove(.x, "_N"))

df_field <- df_field %>% 
  dplyr::select(-c(all_of(Genes_wild), GAPDH, contains("_dct")))

# remove ending _N
df_field <- df_field %>%
  rename_with(~str_remove(.x, "_N"))

# add the new genes to the complete data sets 
lab <- lab %>%
 dplyr::select(-all_of(Gene_lab)) %>%
  left_join(df_lab, by = "Mouse_ID")

field <- field %>%
 dplyr::select(-all_of(Genes_wild)) %>%
  left_join(df_field, by = "Mouse_ID")

hm_norm <- rbind(lab, field)


## --------------------------------------------------------------------------------
 ##save the imputed data 
write.csv(hm_norm, "output_data/1.1.norm_MICE_data_set.csv", row.names = FALSE)

