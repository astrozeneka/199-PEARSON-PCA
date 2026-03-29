library(devtools)
library(Seurat)
library(scPearsonPCA)
library(data.table)
library(ggplot2)
library(ggrepel)

larc_obj <- readRDS("../193-SEURAT-HPC/larc_datasets/larc_merged_6k.rds")

annotation_files <- list.files(
  "splitted/roi_annotated",
  pattern = "\\.csv$", full.names = TRUE
)

cell_sample_list <- lapply(annotation_files, function(f) {
  sname <- basename(f)
  sname <- sub("^LARC_[AB]_", "", sname)
  sname <- sub("_annotated\\.csv$", "", sname)
  df    <- read.csv(f, row.names = 1)
  data.frame(
    cell      = rownames(df)[df$selection_mask == 1],
    sample_id = sname,
    stringsAsFactors = FALSE
  )
})

cell_sample_df <- do.call(rbind, cell_sample_list)
larc_obj <- subset(larc_obj, cells = cell_sample_df$cell)

larc_obj <- NormalizeData(larc_obj)

# Identify highly variable features
larc_obj <- FindVariableFeatures(larc_obj, selection.method = "vst", nfeatures = 2000)
hvgs <- VariableFeatures(larc_obj)

# Calculate summary statistics from ALL genes (not just HVGs)
tc       <- Matrix::colSums(larc_obj[["RNA"]]@counts)
genefreq <- scPearsonPCA::gene_frequency(larc_obj[["RNA"]]@counts)

# Run Pearson PCA on HVGs using full-gene summary statistics
pcaobj <- scPearsonPCA::sparse_quasipoisson_pca_seurat(
    larc_obj[["RNA"]]@counts[hvgs, ]
  , totalcounts = tc
  , grate       = genefreq[hvgs]
  , scale.max   = 10    ## clip pearson residuals > 10 SDs above the mean
  , do.scale    = TRUE  ## scale each gene's pearson residuals to SD = 1
  , do.center   = TRUE  ## center each gene's pearson residuals to mean = 0
)

# ── Downstream analysis ──────────────────────────────────────────────────────

# Build UMAP and nearest-neighbor graph from PCA embedding
umapobj <- scPearsonPCA::make_umap(pcaobj)

# Store reductions and graph in Seurat object
larc_obj[["pearsonpca"]]   <- pcaobj$reduction.data
larc_obj[["pearsonumap"]]  <- umapobj$ump
larc_obj[["pearsongraph"]] <- Seurat::as.Graph(umapobj$grph)

# Unsupervised clustering on the Pearson-based neighbor graph
larc_obj <- Seurat::FindClusters(larc_obj, graph.name = "pearsongraph")
larc_obj@meta.data$pearson_clusters <- larc_obj@meta.data$seurat_clusters

# Plot UMAP coloured by cluster
umapplot <- scPearsonPCA::plot_umap(
    umapreduc  = "pearsonumap"
  , clustercol = "pearson_clusters"
  , semuse          = larc_obj
)
ggsave("plots/umapplot_S19_31776B1.png", plot = umapplot, width = 8, height = 6, dpi = 300)
print(umapplot)

# ── Export cell table ─────────────────────────────────────────────────────────

# UMAP embeddings
umap_df <- as.data.frame(Embeddings(larc_obj[["pearsonumap"]]))
colnames(umap_df) <- c("umap_1", "umap_2")

# Spatial x/y from original annotated CSVs (x_slide_mm / y_slide_mm)
coords_df <- do.call(rbind, lapply(annotation_files, function(f) {
  df <- read.csv(f, row.names = 1)
  data.frame(
    cell = rownames(df),
    x    = df$x_slide_mm,
    y    = df$y_slide_mm,
    stringsAsFactors = FALSE
  )
}))
rownames(coords_df) <- coords_df$cell

# Assemble final table
cell_table <- data.frame(
  cell    = rownames(umap_df),
  x       = coords_df[rownames(umap_df), "x"],
  y       = coords_df[rownames(umap_df), "y"],
  umap_1  = umap_df$umap_1,
  umap_2  = umap_df$umap_2,
  cluster = larc_obj@meta.data[rownames(umap_df), "pearson_clusters"],
  stringsAsFactors = FALSE
)

write.csv(cell_table, "data/cell_coordinates_S19_31776B1.csv", row.names = FALSE)
cat("Exported", nrow(cell_table), "cells to data/cell_coordinates_S19_31776B1.csv\n")
