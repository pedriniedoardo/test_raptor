# libraries ---------------------------------------------------------------
library(tidyverse)
library(Seurat)

# READ IN THE FILES -------------------------------------------------------
# read in the seurat object
s25 <- readRDS("../../data/new_WM_Abisnta/s25.rds")
s26 <- readRDS("../../data/new_WM_Abisnta/s26.rds")
s27 <- readRDS("../../data/new_WM_Abisnta/s27.rds")
s31 <- readRDS("../../data/new_WM_Abisnta/s31.rds")

# Merging Seurat Objects
data.combined.all <- merge(s25, y = c(s26,s27,s31), add.cell.ids = c("s25", "s26", "s27","s31"), project = "AbsintaWmNew")
data.combined.all

# make sure the default assay is the RNA
DefaultAssay(data.combined.all) <- "RNA"

# read in the LUT
LUT_sample <- read_csv("../../data/new_WM_Abisnta/LUT_sample_Absinta_new_WM.csv")

# generate pseudobulk -----------------------------------------------------

# aggregate the sameple per celltype and donor ----------------------------
# get the metadata from the other object
meta <- data.combined.all@meta.data %>%
  rownames_to_column(var = "barcodes")

meta_full <- meta %>%
  # mutate(sample_number_fix = paste0("s",sample)) %>% 
  left_join(LUT_sample,by = c("orig.ident"="sample"))

# add the new metadata
# data.combined.all$clusterCellType <- meta_full$clusterCellType
data.combined.all$pathology_fix <- meta_full$condition
data.combined.all$sample_number_fix <- meta_full$sample
data.combined.all$disease <- case_when(data.combined.all$pathology_fix=="CTRL"~"ctrl",
                                       T~"ms")


# aggregate the sameple per donor -----------------------------------------
# get the metadata from the other object

# confirm the groupign for the metadata
table(data.combined.all$sample_number_fix,data.combined.all$pathology_fix)
table(data.combined.all$disease,data.combined.all$pathology_fix)

# create the column for the aggragation
data.combined.all$pseudobulk2 <- paste0(data.combined.all$pathology_fix,
                                        ".",
                                        data.combined.all$disease,
                                        ".",
                                        data.combined.all$sample_number_fix)

# add it also to the main metadata
meta_full$pseudobulk2 <- paste0(data.combined.all$pathology_fix,
                                ".",
                                data.combined.all$disease,
                                ".",
                                data.combined.all$sample_number_fix)

# save the aggregtated table of counts raw
matrix_counts_alt <- AggregateExpression(object = data.combined.all,
                                         group.by = c("pseudobulk2"),
                                         assays = "RNA",
                                         slot = "counts",
                                         return.seurat = FALSE) %>%
  .$RNA

# save also the normalized table of counts
matrix_counts2_alt <- AggregateExpression(object = data.combined.all,
                                          group.by = c("pseudobulk2"),
                                          assays = "RNA",
                                          slot = "data",
                                          return.seurat = FALSE) %>%
  .$RNA

dim(matrix_counts_alt)
dim(matrix_counts2_alt)
matrix_counts_alt[1:10,1:4]
matrix_counts2_alt[1:10,1:4]

# save the matrix as object
matrix_counts_alt %>%
  saveRDS("../../out/object/pseudobulk_absintaNewWM_counts_wholeSample_raw.rds")
matrix_counts2_alt %>%
  saveRDS("../../out/object/pseudobulk_absintaNewWM_counts_wholeSample_norm.rds")

# save the final meta
write_tsv(meta_full,file = "../../out/table/meta_full_absintaNewWM.tsv")
