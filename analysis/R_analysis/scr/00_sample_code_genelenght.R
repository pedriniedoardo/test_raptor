# libraries ---------------------------------------------------------------
library(biomaRt)
library(tidyverse)

genes <- read.table(text = "          gene_id
1 ENSG00000000003
2 ENSG00000000005
3 ENSG00000000419
4 ENSG00000000457
5 ENSG00000000460
6 ENSG00000000938")

head(genes)

## Build a biomart query 
# In the example below, I use the human gene annotation from Ensembl release 82 located on "sep2015.archive.ensembl.org", More about the ensembl_build can be found on "http://www.ensembl.org/info/website/archives/index.html"

dataset <- "hsapiens_gene_ensembl"

# use a specific annotation for the conversion
# mart = biomaRt::useMart("ENSEMBL_MART_ENSEMBL", dataset = dataset, 
#                         host = paste0("sep2015", ".archive.ensembl.org"), path = "/biomart/martservice", archive = FALSE)
listMarts()
mart <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL",
                        dataset = dataset)

listFilters(mart) %>% head()

annotations <- biomaRt::getBM(mart = mart, attributes=c("ensembl_gene_id", "external_gene_name", "start_position", "end_position")) %>% 
  dplyr::mutate(gene_length = end_position - start_position)

# Filter and re-order gene.annotations to match the order in your input genes list
final.genes <- annotations %>% dplyr::filter(ensembl_gene_id %in% genes$gene_id)
final.genes


# -------------------------------------------------------------------------
# try the same using AnnotaionHub to donwnload the database
library(AnnotationDbi)
library(AnnotationHub)

ah <- AnnotationHub()
# snapshotDate(): 2022-10-26

#queryt the dtabaset to identify the versino of interest
query(ah, pattern = c("Homo Sapiens", "EnsDb"))

# copy the specific code ofr the database of interest
edb_test <- ah[["AH104864"]]

columns(edb_test)

genes

# using the selct function
df_out2 <- select(edb_test, keys=genes$gene_id,
                 columns=c("GENEBIOTYPE","ENTREZID","SYMBOL","GENESEQSTART","GENESEQEND"),
                 keytype="GENEID") %>%
  mutate(ENTREZID = as.character(ENTREZID)) %>% 
  mutate(gene_length = GENESEQEND - GENESEQSTART)

head(df_out2)

# using the mapids function
df_out3 <- 
  genes %>% 
  mutate( biotype = mapIds(edb_test, 
           keys = gene_id, 
           column = c("GENEBIOTYPE"), 
           keytype = "GENEID", 
           multiVals = "first") 
  ) %>% 
  mutate( symbol = mapIds(edb_test, 
                           keys = gene_id, 
                           column = c("SYMBOL"), 
                           keytype = "GENEID", 
                           multiVals = "first") 
  ) %>% 
  mutate(seqEnd = mapIds(edb_test, 
                          keys = gene_id, 
                          column = c("GENESEQEND"), 
                          keytype = "GENEID", 
                          multiVals = "first") 
  ) %>% 
  mutate(seqStart = mapIds(edb_test, 
                         keys = gene_id, 
                         column = c("GENESEQSTART"), 
                         keytype = "GENEID", 
                         multiVals = "first") 

           ) %>% 
  mutate(gene_length = seqEnd - seqStart)

df_out3


