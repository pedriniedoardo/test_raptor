# libraries ---------------------------------------------------------------
library(tidyverse)
library(Seurat)

# READ IN THE FILES -------------------------------------------------------
# read in the seurat object
data.combined.all <- readRDS("../../data/jakel_clean.rds")

# confirm the identity of the dataset
DimPlot(data.combined.all,label = T)

# make sure the default assay is the RNA
DefaultAssay(data.combined.all) <- "RNA"

# read in the LUT
LUT_sample <- read_csv("../../data/LUT_sample_jakell.csv")

# generate pseudobulk -----------------------------------------------------

# aggregate the sameple per celltype and donor ----------------------------
# get the metadata from the other object
meta <- data.combined.all@meta.data %>%
  rownames_to_column(var = "barcodes")
  # add the macroclassification
  # mutate(clusterCellType = case_when(seurat_clusters %in% c(0,1,2,3,6) ~"Oligo",
  #                                    seurat_clusters %in% c(7,15) ~"Neu",
  #                                    seurat_clusters %in% c(8) ~"OPC",
  #                                    seurat_clusters %in% c(11,13) ~"Vas",
  #                                    seurat_clusters %in% c(16) ~"Lym",
  #                                    seurat_clusters %in% c(5,10,17) ~"Imm",
  #                                    seurat_clusters %in% c(4,9,12,14) ~"Ast")) %>%
  # # fix the pathology column
  # mutate(pathology_fix = case_when(pathology == "a_control"~"CTRL",
  #                                  pathology == "b_NAWM"~"NAWM",
  #                                  pathology == "c_chronic_active"~"CA",
  #                                  pathology == "d_chronic_inactive"~"CI",
  #                                  pathology == "e_core"~"CORE"))

meta_full <- meta %>%
  # mutate(sample_number_fix = paste0("s",sample)) %>% 
  left_join(LUT_sample,by = c("sample"))

# add the new metadata
# data.combined.all$clusterCellType <- meta_full$clusterCellType
data.combined.all$pathology_fix <- meta_full$condition
data.combined.all$sample_number_fix <- meta_full$sample

# confirm the groupign for the metadata
# table(data.combined.all$clusterCellType,data.combined.all$sample_number_fix,data.combined.all$pathology_fix,data.combined.all$disease)
# table(data.combined.all$sample_number_fix,data.combined.all$pathology_fix)

# # create the column for the aggragation
# data.combined.all$pseudobulk <- paste0(data.combined.all$clusterCellType,
#                                        ".",
#                                        data.combined.all$pathology_fix,
#                                        ".",
#                                        data.combined.all$disease,
#                                        ".",
#                                        data.combined.all$sample_number_fix)
# 
# # add it also to the main metadata
# meta_full$pseudobulk <- paste0(meta_full$clusterCellType,
#                                ".",
#                                meta_full$pathology_fix,
#                                ".",
#                                meta_full$disease,
#                                ".",
#                                meta_full$sample_number_fix)
# 
# # save the aggregtated table of counts raw
# matrix_counts <- AggregateExpression(object = data.combined.all,
#                                      group.by = c("pseudobulk"),
#                                      assays = "RNA",
#                                      slot = "counts",
#                                      return.seurat = FALSE) %>%
#   .$RNA
# 
# # save also the normalized table of counts
# matrix_counts2 <- AggregateExpression(object = data.combined.all,
#                                       group.by = c("pseudobulk"),
#                                       assays = "RNA",
#                                       slot = "data",
#                                       return.seurat = FALSE) %>%
#   .$RNA
# 
# dim(matrix_counts)
# dim(matrix_counts2)
# matrix_counts[1:10,1:10]
# matrix_counts2[1:10,1:10]
# 
# # save the matrix as object
# matrix_counts %>%
#   saveRDS("../../out/object/pseudobulk_all20_integrated_clean_metadata_counts_PatientCellType_raw.rds")
# matrix_counts2 %>%
#   saveRDS("../../out/object/pseudobulk_all20_integrated_clean_metadata_counts_PatientCellType_norm.rds")

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
meta_full$pseudobulk2 <- paste0(meta_full$pathology_fix,
                                ".",
                                meta_full$disease,
                                ".",
                                meta_full$sample_number_fix)

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
matrix_counts_alt[1:10,1:10]
matrix_counts2_alt[1:10,1:10]

# save the matrix as object
matrix_counts_alt %>%
  saveRDS("../../out/object/pseudobulk_jakell_counts_wholeSample_raw.rds")
matrix_counts2_alt %>%
  saveRDS("../../out/object/pseudobulk_jakell_counts_wholeSample_norm.rds")

# save the final meta
write_tsv(meta_full,file = "../../out/table/meta_full_jakell.tsv")
