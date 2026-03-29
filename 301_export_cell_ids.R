library(Seurat)

larc_obj <- readRDS("../193-SEURAT-HPC/larc_datasets/larc_merged_6k.rds")
write.csv(larc_obj@meta.data, paste0("data/cells_meta.csv"), row.names = TRUE)
