# AIM ---------------------------------------------------------------------
# preprocess the reference dataset

# -------------------------------------------------------------------------
# source("scr/00_tutorial_load_libraries.R")

# read in the file --------------------------------------------------------
# load the LUT for the gene to proba annotations
load("../../data/hom_hugene1.RData") # run load_hom_hugene1.R to generate

# # load chen2016 data (human brain) https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE71620
# geo_id <- "GSE71620"
# geo_obj <- getGEO(geo_id)[[1]]
# # save the object for future access
# saveRDS(geo_obj,"../../data/geo_obj_GSE71620.rds")
# load the obejct previously saved
geo_obj <- readRDS("../../data/geo_obj_GSE71620.rds")

# wrangling ---------------------------------------------------------------
# unlog the data
g <- 2^exprs(geo_obj)
# add the genes ids instead of the probes IDs
g2 <- format_ids(g, hom_hugene1, from = 2, to = 1)
# Kept 29106 out of 33297 - aggregated into 24601
# notice the diffenece in size before and after changing the ids
dim(g)
dim(g2)

# explore one specific gene for which we had more affy ids
id_test <- hom_hugene1 %>% 
  data.frame() %>% 
  dplyr::filter(hsa_ens_id %in% c("ENSG00000210049")) %>% 
  dplyr::pull(affy_id)

# the aggregation is the average of all the probes for the same gene
g[id_test,1:10]
g[id_test,1:10] %>% 
  apply(MARGIN = 2,mean)

# confim with the aggreagated dataset
g2["ENSG00000210049",1:10]

# pull the metadata from the original object and adjust
p <- pData(geo_obj)
p <- p[, c("title", "geo_accession", "age:ch1", "Sex:ch1", "pmi:ch1", "ph:ch1", "rin:ch1", "tod:ch1")]
colnames(p) <- c("title", "geo_accession", "age", "sex", "pmi", "ph", "rin", "tod")
head(p)
p$age <- as.numeric(p$age) # age in years
p$ph <- as.numeric(p$ph) # ph of sample
p$rin <- as.numeric(p$rin) # rna integrity nb
p$pmi <- as.numeric(p$pmi) # post-mortem interval (hours)
p$tod <- as.numeric(p$tod) # time of death
p$sex <- factor(p$sex, levels = c("M", "F"))
p$tissue <- factor(gsub("(BA\\d+)\\ssample\\s.*", "\\1", p$title))
p$snum <- gsub("(BA\\d+)\\ssample\\s(.*)", "\\2", p$title)

# check the distribution and the range of the age
p %>% 
  ggplot(aes(x=age))+geom_histogram()
  
# remove samples with rna integrity < 7
g3 <- g2[, p$rin > 7]
p3 <- p[p$rin > 7, ]

# remove outliers
cc <- cor(g3)
diag(cc) <- NA
qs <- apply(cc, 1, quantile, probs= 0.99, na.rm = T)
thr <- median(qs, na.rm=T) - 2*sd(qs, na.rm = T)
keep <- qs > thr
boxplot(cc, col = 1+ !keep, border= 1+!keep)
table(keep)

# select either brain region dissected
s47 <- which(p3$tissue == "BA47" & keep)
s11 <- which(p3$tissue == "BA11" & keep)

g47 <- g3[, s47]
g11 <- g3[, s11]

p47 <- p3[s47, ]
p11 <- p3[s11, ]

# build a list with expression (unlogged) and metadata
dschen2016 <- list(g47=g47, g11=g11, p47=p47, p11=p11)
save(dschen2016, file = file.path(data_folder, "dschen2016.RData"))

# # cleanup
# rm(list = ls())
# gc()
