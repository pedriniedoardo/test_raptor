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
g11 <- log1p(limma::normalizeBetweenArrays(dschen2016$g11, method = "quantile"))
df_test_norm <- log1p(limma::normalizeBetweenArrays(df_test, method = "quantile"))

# check the distributions
boxplot(dschen2016$g11[,1:10])
boxplot(g11[,1:10])
boxplot(df_test_norm[,])

# compute monotony
mon_g11 <- abs(apply(g11, 1, cor, y=dschen2016$p11$age, method="spearman"))
sel11 <- mon_g11>sqrt(.25)
pca_11 <- summary(prcomp(t(g11[sel11,]), center = T, scale. = F, rank. = 20))

# build the reference -----------------------------------------------------
# build references with all the samples
# selsamp <- c(T,F)

# build the model for a subset of g47 samples
m11 <- ge_im(X = g11[sel11,], p = dschen2016$p11,
             formula = "X ~ s(age, bs = 'cr')",
             dim_red = "pca", nc = 1)

# build an even rank of 500 steps in ages from the min to the max
n.inter <- 500
ndat <- data.frame(age = seq(min(dschen2016$p11$age),
                             max(dschen2016$p11$age), l=n.inter))

# list containint both predicted expression values and relative ages
r11 <- list(interpGE=predict(m11, ndat),
            time.series=ndat$age)
boxplot(r11$interpGE[,1:10])

# stage samples on BA47 & BA11 references
# use r47 as reference
ae_df_test_norm_r11 <- ae(df_test_norm, r11$interpGE, r11$time.series)

# add the estimated ages to the original dataset as metadata
LUT_sample$ae11 <- ae_df_test_norm_r11$age.estimates[, 1]

# plotting ----------------------------------------------------------------
LUT_sample %>% 
  mutate(delta_age = ae11 - age.x) %>% 
  mutate(pathology_fix = factor(pathology_fix,levels = c("CTRL","NAWM","CI","CA","CORE"))) %>% 
  ggplot(aes(x=pathology_fix,y=delta_age)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2),alpha = 0.5)+theme_bw()+geom_hline(yintercept = 0,linetype = "dashed",col="gray")
