# AIM ---------------------------------------------------------------------
# the aim of the script is fo download and create the hom_hugene1 object. A table of conversion between ensamble gene ids and affy ids

# -------------------------------------------------------------------------
# source("scr/R/00_tutorial_load_libraries.R")

# build the query ---------------------------------------------------------
request_start <- '<?xml version="1.0" encoding="UTF-8"?>
  <!DOCTYPE Query>
  <Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
  
  <Dataset name = "hsapiens_gene_ensembl" interface = "default" >
  <Attribute name = "ensembl_gene_id" />'

probe_id <- '<Attribute name = "affy_hugene_1_0_st_v1" />'

request_end <- '</Dataset></Query>'

data_folder <- "../../data/"

# fetch the data ----------------------------------------------------------
# XML request to biomart for human affy probe-gene list
system(paste0("wget -O ", data_folder,
              'hom_affy.txt \'http://www.ensembl.org/biomart/martservice?query=',
              request_start, probe_id, request_end, "'"))

# wrangling ---------------------------------------------------------------
hom_hugene1 <- read.table(paste0(data_folder, "hom_affy.txt"), h=F, as.is = T, sep = "\t")
colnames(hom_hugene1) <- c("hsa_ens_id", "affy_id")
head(hom_hugene1)
hom_hugene1$affy_id <- as.character(hom_hugene1$affy_id)

# save the object ---------------------------------------------------------
save(hom_hugene1, file = paste0(data_folder, "hom_hugene1.RData"), compress = "xz")
