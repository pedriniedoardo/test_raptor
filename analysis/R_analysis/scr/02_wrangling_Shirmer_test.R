# libraires ---------------------------------------------------------------
library(tidyverse)
library(AnnotationDbi)
library(AnnotationHub)

# function implementation -------------------------------------------------
rpkm <- function(counts, lengths) {
  rate <- counts / lengths 
  rate / sum(counts) * 1e6
}

tpm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(rate) * 1e6
}

# # sample usage
# genes <- data.frame(
#   Gene = c("A","B","C","D","E"),
#   Length = c(100, 50, 25, 5, 1)
# )
# 
# counts <- data.frame(
#   S1 = c(80, 10,  6,  3,   1),
#   S2 = c(20, 20, 10, 50, 400)
# )
# 
# rpkms <- apply(counts, 2, function(x) rpkm(x, genes$Length))
# tpms <- apply(counts, 2, function(x) tpm(x, genes$Length))
# 
# # Sample means should be equal.
# 
# colSums(rpkms)
# colSums(tpms)
# 
# colMeans(rpkms)
# colMeans(tpms)
# 
# rate <- counts$S1/genes$Length
# rate/sum(rate)*10^6
# 
# tpms

# read in the data --------------------------------------------------------
# read in the raw total table of counts per sample
df_counts <- readRDS("../../out/object/pseudobulk_jakell_counts_wholeSample_raw.rds") %>% 
  data.frame() %>% 
  rownames_to_column("symbol")

# # read in the metadata with age informations
# LUT_sample <- read_tsv("../../out/table/meta_full.tsv") %>% 
#   group_by(pseudobulk2,age.x,pathology_fix,disease,sample_number_fix) %>% 
#   summarise()

# wrangling ---------------------------------------------------------------
# pull the gene names and collect the metadata needed to make the table as TPM
LUT_genes <- data.frame(symbol = df_counts$symbol)

# pull the annotation
ah <- AnnotationHub()
#queryt the dtabaset to identify the versino of interest
query(ah, pattern = c("Homo Sapiens", "EnsDb"))
# copy the specific code ofr the database of interest
edb_test <- ah[["AH109606"]]
# annotations available
columns(edb_test)

# using the mapids function to  pull the info of interest
LUT_genes_full <- 
  LUT_genes %>% 
  mutate(biotype = mapIds(edb_test,
                          keys = symbol,
                          column = c("GENEBIOTYPE"),
                          keytype = "SYMBOL", 
                          multiVals = "first") 
  ) %>% 
  mutate(ensembl = mapIds(edb_test,
                          keys = symbol,
                          column = c("GENEID"),
                          keytype = "SYMBOL", 
                          multiVals = "first") 
  ) %>% 
  mutate(seqEnd = mapIds(edb_test, 
                         keys = symbol, 
                         column = c("GENESEQEND"), 
                         keytype = "SYMBOL", 
                         multiVals = "first") 
  ) %>% 
  mutate(seqStart = mapIds(edb_test, 
                           keys = symbol, 
                           column = c("GENESEQSTART"), 
                           keytype = "SYMBOL", 
                           multiVals = "first") 
         
  ) %>% 
  mutate(gene_length = seqEnd - seqStart)

# save the gene annotation
LUT_genes_full %>% 
  write_tsv("../../out/table/lUT_gene_full_Jakell.tsv")

# shortlist the GOI
LUT_genes_full_filter <- LUT_genes_full %>% 
  # filter only the genes that have an ensable
  dplyr::filter(!is.na(ensembl)) %>% 
  # filter only the genes for which I have the total length
  dplyr::filter(!is.na(gene_length))

# shortlist the table of expression
df_counts_full <- LUT_genes_full_filter %>% 
  left_join(df_counts,"symbol")

saveRDS(df_counts_full,"../../out/object/pseudobulk_jakell_wholeSample_raw_filter_full.rds")

# select only the info of interest
df_counts_full_ensembl <- df_counts_full %>% 
  dplyr::select(-c(1,2,4,5,6)) %>% 
  column_to_rownames("ensembl")

saveRDS(df_counts_full_ensembl,"../../out/object/pseudobulk_jakell_wholeSample_raw_filter_ensembl.rds")

# calculate the TPM
genes <- data.frame(ensembl = df_counts_full$ensembl,
                    gene_length = df_counts_full$gene_length) %>%
  rownames_to_column("symbol") %>% 
  column_to_rownames("ensembl")

df_TPM <- apply(df_counts_full_ensembl,MARGIN = 2,function(x){
  tpm(x, genes$gene_length)
  # print(x)
})
saveRDS(df_TPM,"../../out/object/pseudobulk_jakell_wholeSample_TPM_filter_ensembl.rds")
