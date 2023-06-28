library("utils")
library("parallel")
# library("readxl") # for load_tissue_specific_genesets.R and load_dshall2010.R only

library("stats")
library("ROCR")
library("ica")
library("splines")
#library("randomForest") # for qtl_dsrockman2010.R only

library("RColorBrewer")
library("viridis")
library("ggplot2")
library("ggpubr")
library("vioplot")


# Bioconductor packages (install instructions on https://www.bioconductor.org/install/)
library("limma")
library("edgeR")
library("Biobase")
library("GEOquery") 
library("biomaRt")
# library("affy")    # for load_dslehrbach2012.R only
# library("gcrma")   # for load_dslehrbach2012.R only


# RAPToR packages & associated data-packages (loaded within appropriate scripts)
library("RAPToR")
# library("wormRef")    # install from www.github.com/LBMC/wormRef 
# library("drosoRef")   #              www.github.com/LBMC/drosoRef  
# library("zebraRef")   #              www.github.com/LBMC/zebraRef
# library("mouseRef")   #              www.github.com/LBMC/mouseRef




## Utility functions / variables
data_folder <- "../../data/"

col.palette <- c("black", "firebrick", "royalblue", "forestgreen", 
                 "gold2", "magenta3", "grey50", "tan4", "cyan4", "deeppink")

source("scr/plot_functions.R")

# Convert to TPM
raw2tpm <- function(rawcounts, genelengths){
  if(nrow(rawcounts) != length(genelengths))
    stop("genelengths must match nrow(rawcounts).")
  x <- rawcounts/genelengths
  return(t( t(x) * 1e6 / colSums(x) ))
}

fpkm2tpm <- function(fpkm){
  return(exp(log(fpkm) - log(colSums(fpkm)) + log(1e6)))
}


