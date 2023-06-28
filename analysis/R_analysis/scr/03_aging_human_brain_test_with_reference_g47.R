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
load("../../data/test/dschen2016.RData")

# fnpx <- "aging_hom_brain_n500_"
# ftype <- "pdf"

# read in the test dataset
df_test <- readRDS("../../out/object/pseudobulk_all20_wholeSample_TPM_filter_ensembl.rds")
# read in the metadata
LUT_sample <- read_tsv("../../out/table/meta_full.tsv") %>% 
  group_by(pseudobulk2,age.x,pathology_fix,disease,sample_number_fix) %>% 
  summarise() %>% 
  ungroup() %>%
  # slice_sample(n = 20) %>% 
  # order the meta as in the column of the df_test
  dplyr::slice(match(colnames(df_test),.$pseudobulk2))

# wrangling ---------------------------------------------------------------
# quantile normalization of the matrices
g47 <- log1p(limma::normalizeBetweenArrays(dschen2016$g47, method = "quantile"))
df_test_norm <- log1p(limma::normalizeBetweenArrays(df_test, method = "quantile"))

# check the distributions
boxplot(dschen2016$g47[,1:10])
boxplot(g47[,1:10])
boxplot(df_test_norm[,])

# compute monotony
mon_g47 <- abs(apply(g47, 1, cor, y=dschen2016$p47$age, method="spearman"))
sel47 <- mon_g47>sqrt(.25)
pca_47 <- summary(prcomp(t(g47[sel47,]), center = T, scale. = F, rank. = 20))

# build the reference -----------------------------------------------------
# build references with all the samples
# selsamp <- c(T,F)

# build the model for a subset of g47 samples
m47 <- ge_im(X = g47[sel47,], p = dschen2016$p47,
             formula = "X ~ s(age, bs = 'cr')",
             dim_red = "pca", nc = 1)

# build an even rank of 500 steps in ages from the min to the max
n.inter <- 500
ndat <- data.frame(age = seq(min(dschen2016$p47$age),
                             max(dschen2016$p47$age), l=n.inter))

# list containint both predicted expression values and relative ages
r47 <- list(interpGE=predict(m47, ndat),
            time.series=ndat$age)
boxplot(r47$interpGE[,1:10])

# stage samples on BA47 & BA11 references
# use r47 as reference
ae_df_test_norm_r47 <- ae(df_test_norm, r47$interpGE, r47$time.series)

# add the estimated ages to the original dataset as metadata
LUT_sample$ae47 <- ae_df_test_norm_r47$age.estimates[, 1]

# plotting ----------------------------------------------------------------
LUT_sample %>% 
  mutate(delta_age = ae47 - age.x) %>% 
  mutate(pathology_fix = factor(pathology_fix,levels = c("CTRL","NAWM","CI","CA","CORE"))) %>% 
  ggplot(aes(x=pathology_fix,y=delta_age)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2),alpha = 0.5)+theme_bw()+geom_hline(yintercept = 0,linetype = "dashed",col="gray")

