###Single Cell Object Creation and Integration###
#This script creates a Seurat object from cellranger filtered h5 files, marks doublets, and integrates the data using Seurat's RPCA method.
#It also performs quality control, normalization, variable feature selection, scaling, PCA, clustering, and UMAP visualization.

#Code was initially created by Pankaj Chopra and adapted by Lisa Blackmer-Raynolds

#Load required packages----
library(hdf5r)
library(glmGamPoi)
library(scCustomize)
library(DoubletFinder)
library(dplyr)
library(tidyverse)
library(Seurat)
library(stringr)

#Define variables----
EXPT <- '21173FL-04'
NSAMPLES <- 15 #Number of samples (note the original object contained 15 samples, since it included samples from 2 brain regions, however, only the 8 hippocampus samples were used in subsequent analyses)
MAX_FEATURES <- 7500 # Maximum genes/cell
MIN_FEATURES <- 0 # Minimum genes/cell
MAX_MITOCHONDRIAL_PERCENTAGE <- 6 # Maximum percent of mitochondrial genes (in QC)
MIN_CELLS <- 15 # Minimum number of cells expressing a gene (in DE analysis)
MIN_PCT <- 0.1 # Minimum percent of cells expressing gene (in DE analysis)

RESOLUTION_PARAMS <- seq(0,0.8,.2)
RESOLUTION <- 0.4
PLOT.GROUPS <- c("treatment","sample","batch", "region")

SEED <- 12791
set.seed(SEED)

options(future.globals.maxSize = 8000 * 1024^2)

samples <- read_csv("SampleList.csv")

#Define functions----
## OBJECTIVE: Mark each cell as Singlet/Doublet using DoubletFinder
## INPUT: 1. seurat object
##        2. Number of PCs to be used
## OUTPUT: 1. seurat object with each cell marked as Singlet/Doublet using DoubletFinder
## Ref: pankaj.organoid.3q29.v1.39.functions.R; ssloan.admera.functions.v6.R; code from Pankaj Chopra

fmarkDoublets.v1.01 <- function(sobj.tid,npcs=30){

  # https://github.com/chris-mcginnis-ucsf/DoubletFinder
  #sobj.tid <- SCTransform(sobj.tid)
  sobj.tid <- NormalizeData(sobj.tid)
  sobj.tid <- FindVariableFeatures(sobj.tid, selection.method = "vst", nfeatures = 2000)
  sobj.tid <- ScaleData(sobj.tid)
  sobj.tid <- RunPCA(sobj.tid)
  sobj.tid <- RunUMAP(sobj.tid, dims = 1:npcs)
  sobj.tid <- FindNeighbors(sobj.tid, dims = 1:npcs, verbose = FALSE)
  sobj.tid <- FindClusters(sobj.tid, verbose = FALSE)
  
  #################
  ### DoubletFinder 
  #################
  sweep.res.list <- paramSweep(sobj.tid, PCs = 1:npcs, sct = F)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  
  annotations <- sobj.tid@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- sobj.tid@meta.data$ClusteringResults
  multiplex.rate <- 0.008 * (floor(ncol(sobj.tid)/1000)) # 10X
  nExp_poi <- round(multiplex.rate*nrow(sobj.tid@meta.data))  ## 
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  ## Default pK at this stage
  sobj.tid <- doubletFinder(sobj.tid, PCs = 1:npcs, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = F)
  pANNcol <- grep("pANN",colnames(sobj.tid@meta.data))
  pANNcol <- colnames(sobj.tid@meta.data)[pANNcol]
  sobj.tid <- doubletFinder(sobj.tid, PCs = 1:npcs, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = pANNcol, sct = F)
  
  #gidx <- grep("DF.classification",colnames(sobj.tid@meta.data))[2] ## high confidence doublets
  gidx <- grep("DF.classification",colnames(sobj.tid@meta.data))[1] ## high confidence doublets - edit 6/6/2024
  sobj.tid@meta.data$DF.classification <- sobj.tid@meta.data[,gidx]
  #sobj.tid@meta.data <- sobj.tid@meta.data[,c('orig.ident','nCount_RNA','nFeature_RNA','DF.classification')]
  
  return(sobj.tid)
  
}

#Make object----
# Initialize an empty list to store individual Seurat objects
tobj.list <- list()

for (i in 1:15) {
  # Construct the file id with leading zeros if necessary
  file_id <- paste0("21173FL-04-", sprintf("%02d", i))
  
  # Construct the full path to the file
  filename <- paste0(file_id, ".cellbender_filtered_seurat.h5")
  
  # Check if the file exists before attempting to read it
  if (file.exists(filename)) {
    # Read the file
    cell_bender_mat <- Read10X_h5(filename = filename)
    
    # Create individual Seurat object
    tobj <- CreateSeuratObject(counts = cell_bender_mat, names.field = 1, names.delim = "_", min.cells = MIN_CELLS)
    
    # Add metadata to the Seurat object from the samples data frame
    tobj@meta.data$treatment <- as.character(samples$treatment[i])
    tobj@meta.data$sample <- as.character(samples$sample[i])
    tobj@meta.data$region <- as.character(samples$region[i])
    
    # Add mitochondrial gene percentage
    tobj[["percent.mt"]] <- PercentageFeatureSet(tobj, pattern = "^mt-")
    
    # Subset the Seurat object based on quality control metrics
    tobj <- subset(tobj, subset = nFeature_RNA > MIN_FEATURES & nFeature_RNA < MAX_FEATURES & percent.mt < MAX_MITOCHONDRIAL_PERCENTAGE)
    
    # Mark and remove doublets
    tobj.marked.doublets <- fmarkDoublets.v1.01(tobj, npcs = 15)
    singlet.idx <- tobj.marked.doublets$DF.classification == "Singlet"
    tobj <- tobj[, singlet.idx]
    
    # Clean up memory
    rm(tobj.marked.doublets)
    
    # Store the Seurat object in the list
    tobj.list[[file_id]] <- tobj
    
    # Print the dimensions of the matrix
    print(paste("Dimensions of", file_id, ":", dim(cell_bender_mat)))
  } else {
    print(paste("File does not exist:", filename))
  }}

#Merge seurat objects into one----
file_ids <- paste0(EXPT,'-',str_pad(1:15,2,pad="0")) 

tobj.merged <- Merge_Seurat_List(
     list_seurat=tobj.list,
     add.cell.ids = file_ids,
     merge.data = TRUE
   )

#Run without integration----
tobj.merged <- NormalizeData(tobj.merged)
tobj.merged <- FindVariableFeatures(tobj.merged)
tobj.merged <- ScaleData(tobj.merged)
tobj.merged <- RunPCA(tobj.merged)
ElbowPlot(tobj.merged)
tobj.merged <- FindNeighbors(tobj.merged, dims = 1:20, reduction = "pca")
tobj.merged <- FindClusters(tobj.merged, resolution = RESOLUTION, cluster.name = "unintegrated_clusters")
tobj.merged <- RunUMAP(tobj.merged, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
p0 <- DimPlot(tobj.merged, reduction = "umap.unintegrated", group.by = c(PLOT.GROUPS,"unintegrated_clusters"),combine = FALSE)

#Seurat integration----
tobj.merged <- IntegrateLayers(
   object = tobj.merged, method = RPCAIntegration,
   orig.reduction = "pca", new.reduction = "integrated.rpca",
   verbose = FALSE
 )

tobj.merged.seurat <- FindNeighbors(tobj.merged, reduction = "integrated.rpca", dims = 1:30)
tobj.merged.seurat.ctree <- FindClusters(tobj.merged.seurat, resolution = RESOLUTION_PARAMS,random.seed = SEED)
tobj.merged.seurat <- FindClusters(tobj.merged.seurat, resolution = RESOLUTION, cluster.name = "rpca_clusters",random.seed = SEED)
md <- tobj.merged.seurat@meta.data
md$rpca_clusters <- factor(md$rpca_clusters,levels=sort(as.numeric(as.character(unique(md$rpca_clusters)))))
tobj.merged.seurat@meta.data <- md
tobj.merged.seurat.umap <- RunUMAP(tobj.merged.seurat, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca")

#Display integrated UMAP
DimPlot(tobj.merged.seurat, reduction = "umap.rpca", label = TRUE)

#Check QC by cell cluster----
tobj.merged.seurat$log10GenesPerUMI <- log10(tobj.merged.seurat$nFeature_RNA + 1) / log10(tobj.merged.seurat$nCount_RNA + 1)

#Feature plots for QC metrics:
FeaturePlot(tobj.merged.seurat, c("nFeature_RNA","nCount_RNA", "percent.mt", "log10GenesPerUMI"), reduction = "umap.rpca", pt.size = 0.3, order = TRUE)

#Cell cycle genes:
FeaturePlot(tobj.merged.seurat, reduction = "umap", c("Mki67", "Ccna2", "Ccnb2", "Pcna", "Ccnd1"), pt.size = 0.3, order = TRUE)

#Identify cell types----
#Re-cluster with lower resolution to get fewer cell types
tobj.merged.seurat <- FindNeighbors(tobj.merged.seurat, reduction = "integrated.rpca", dims = 1:40)
tobj.merged.seurat <- FindClusters(tobj.merged.seurat, resolution = .01)
tobj.merged.seurat <- RunUMAP(tobj.merged.seurat, dims = 1:40, reduction = "integrated.rpca")

DimPlot(tobj.merged.seurat, reduction = "umap")

#Find the markers for each cluster:
tobj.merged.seurat <- JoinLayers(tobj.merged.seurat)
markers <- FindAllMarkers(object = tobj.merged.seurat, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
write.csv(markers,'markers.csv')

#Number of cells per cluster and sample:
n_cells <- table(tobj.merged.seurat@meta.data$seurat_clusters, tobj.merged.seurat@meta.data$sample)
write.csv(n_cells, "n_cells.csv")

#Rename clusters based on markers:
  #Markers were compared to gene expression data from the Allen Brain Cell Atlas to identify cell types

tobj.merged.seurat <- RenameIdents(object = tobj.merged.seurat,
  "0" = "inhibitory_neurons",
  "1" = "excitatory_neurons",
  "2" = "m_oligos",
  "3" = "excitatory_neurons",
  "4" = "astrocytes",
  "5" = "OPCs",
  "6" = "excitatory_neurons",
  "7" = "immune",
  "8" = "vasculature",
  "9" = "excitatory_neurons",
  "10" = "vasculature",
  "11" = "vasculature"
)

tobj.merged.seurat$cell_type <- paste(Idents(tobj.merged.seurat))
DimPlot(tobj.merged.seurat, reduction = "umap", label = FALSE)

#Save seurat objects broken up by region (note, while both regions were used for clustering, the hippocampus region is the only one used in subsequent analyses due to poor dissection quality of midbrain tissue)
Idents(tobj.merged.seurat) <- "region"

hpc <- subset(tobj.merged.seurat, idents = c("hpc"))
saveRDS(hpc, file = "hpc.rds")

mid <- subset(tobj.merged.seurat, idents = c("mid"))
saveRDS(mid, file = "mid.rds")