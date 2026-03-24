# load dev library
library(devtools)
library(Seurat)
library(scPearsonPCA)

# Append to dockerfile
#remotes::install_github("Nanostring-Biostats/CosMx-Analysis-Scratch-Space",
#                         subdir = "_code/scPearsonPCA", ref = "Main")
#remotes::install_github(
#    "Nanostring-Biostats/CosMx-Analysis-Scratch-Space",
#    subdir = "_code/scPearsonPCA",
#    ref = "Main"
# )

larc_obj <- readRDS("../193-SEURAT-HPC/larc_datasets/larc_merged_6k.rds")


larc_obj <- NormalizeData(larc_obj)

# Identify highly variable features
larc_obj <- FindVariableFeatures(larc_obj, selection.method = "vst", nfeatures = 2000)
hvgs <- VariableFeatures(larc_obj)

# Calculate summary statistics from ALL genes (not just HVGs)
tc       <- Matrix::colSums(larc_obj[["RNA"]]@counts)
genefreq <- scPearsonPCA::gene_frequency(larc_obj[["RNA"]]@counts)

# Run Pearson PCA on HVGs using full-gene summary statistics
pcaobj <-
scPearsonPCA::sparse_quasipoisson_pca_seurat(
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
  , larc_objuse     = larc_obj
)
print(umapplot)

