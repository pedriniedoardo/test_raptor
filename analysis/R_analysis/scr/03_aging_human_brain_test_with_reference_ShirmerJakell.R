# libraries ---------------------------------------------------------------
library(tidyverse)
library(RAPToR)

# functions ---------------------------------------------------------------
# stats function
sts <- function(x, y, sel=T){
  x <- x[sel]
  y <- y[sel]
  
  lmi <- lm(y~x)
  res <- residuals(lmi)
  return(list(lm=lmi, 
              r2 = summary(lmi)$adj.r.squared, 
              rho = cor(x,y, method="spearman"),
              rmse = sqrt(sum(res**2)/length(res))))
}

# -------------------------------------------------------------------------
# source("scr/R/00_tutorial_load_libraries.R")

# read in the data --------------------------------------------------------
# load the sample reference
# load("../../data/test/dschen2016.RData")

# fnpx <- "aging_hom_brain_n500_"
# ftype <- "pdf"

# read in the data fro building the reference
ref01 <- readRDS("../../out/object/pseudobulk_jakell_wholeSample_TPM_filter_ensembl.rds")
ref02 <- readRDS("../../out/object/pseudobulk_shirmer_wholeSample_TPM_filter_ensembl.rds")

# read in the test dataset
df_test <- readRDS("../../out/object/pseudobulk_all20_wholeSample_TPM_filter_ensembl.rds")

# read in the metadata
LUT_ref01 <- read_tsv("../../out/table/meta_full_jakell.tsv") %>% 
  group_by(pseudobulk2,age) %>% 
  summarise() %>% 
  ungroup() %>%
  # slice_sample(n = 20) %>% 
  # order the meta as in the column of the df_test
  dplyr::slice(match(colnames(ref01),.$pseudobulk2))

LUT_Absinta <- read_tsv("../../out/table/meta_full.tsv") %>% 
  group_by(pseudobulk2,age.x,pathology_fix,disease,sample_number_fix) %>% 
  summarise() %>% 
  ungroup() %>%
  # slice_sample(n = 20) %>% 
  # order the meta as in the column of the df_test
  dplyr::slice(match(colnames(df_test),.$pseudobulk2))

# wrangling ---------------------------------------------------------------
# define the ref using the controls samples

id_ref <- df_test %>% 
  colnames() %>% 
  str_detect(pattern = "CTRL")

df_ref <- df_test[,id_ref]

meta_ref <- LUT_sample[id_ref,] %>% 
  dplyr::rename(age = age.x) %>% 
  mutate(rowname = pseudobulk2) %>% 
  column_to_rownames()

# quantile normalization of the matrices
df_ref_norm <- log1p(limma::normalizeBetweenArrays(df_ref, method = "quantile"))
df_test_norm <- log1p(limma::normalizeBetweenArrays(df_test, method = "quantile"))

# check the distributions
boxplot(df_ref)
boxplot(df_ref_norm)
boxplot(df_test_norm)

# compute monotony
mon <- abs(apply(df_ref_norm, 1, cor, y=meta_ref$age, method="spearman"))
sel_mon <- mon>sqrt(.25)
# if it is NA make it FALSE
sel_mon2 <- case_when(sel_mon==T~T,
                      sel_mon==F~F,
                      is.na(sel_mon)~F)
pca_df_ref <- summary(prcomp(t(df_ref_norm[sel_mon2,]), center = T, scale. = F, rank. = 20))

# build the reference -----------------------------------------------------
# build references with all the samples
# selsamp <- c(T,F)

# build the model for a subset of g47 samples
m_df_ref <- ge_im(X = df_ref_norm[sel_mon2,],method = "glm",
                  p = meta_ref,
                  formula = "X ~ age",
                  dim_red = "pca", nc = 1)

# m_df_ref <- ge_im(X = df_ref_norm[sel_mon2,],method = "limma",
#                   p = meta_ref,
#                   formula = "X ~ age",
#                   dim_red = "pca", nc = 1)

# m_df_ref <- ge_im(X = df_ref_norm[sel_mon2,],method = "gam",
#                   p = meta_ref,
#                   formula = "X ~ s(age, bs = 'cr')",
#                   dim_red = "pca", nc = 1)

# build an even rank of 500 steps in ages from the min to the max
n.inter <- 500
ndat <- data.frame(age = seq(min(meta_ref$age),
                             max(meta_ref$age), l=n.inter))

# list containint both predicted expression values and relative ages
r_df_ref <- list(interpGE=predict(m_df_ref, ndat),
                 time.series=ndat$age)
boxplot(r_df_ref$interpGE[,1:10])

# stage samples on BA47 & BA11 references
# use r47 as reference
ae_df_test_norm_r_df_ref <- ae(df_test_norm, r_df_ref$interpGE, r_df_ref$time.series)

# add the estimated ages to the original dataset as metadata
LUT_sample$ae_r_df_ref <- ae_df_test_norm_r_df_ref$age.estimates[, 1]

# plotting ----------------------------------------------------------------
LUT_sample %>% 
  mutate(delta_age = ae_r_df_ref - age.x) %>% 
  mutate(pathology_fix = factor(pathology_fix,levels = c("CTRL","NAWM","CI","CA","CORE"))) %>% 
  ggplot(aes(x=pathology_fix,y=delta_age)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2),alpha = 0.5)+theme_bw()+geom_hline(yintercept = 0,linetype = "dashed",col="gray")

