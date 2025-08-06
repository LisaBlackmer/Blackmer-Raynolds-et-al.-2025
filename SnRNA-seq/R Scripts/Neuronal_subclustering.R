### Sub-clustering and analysis of just neurons ####
## This script performs sub-clustering and analysis on neurons subsetted from my seurat object to identify different neuronal populations and their responses to treatment.

# Code created by Lisa Blackmer-Raynolds

#Load required packages----
library(Seurat)
library(tidyverse)
library(dplyr)
library(readr)
library(writexl)
library(RColorBrewer)

#Make an object with just neurons----
Idents(hpc) <- "cell_type"
neurons_hpc <- subset(hpc, idents = c("excitatory_neurons","inhibitory_neurons"))

#Perform clustering----
neurons_hpc <- FindNeighbors(neurons_hpc, reduction = "integrated.rpca", dims = 1:40)
neurons_hpc <- FindClusters(neurons_hpc, resolution = .005)
neurons_hpc <- RunUMAP(neurons_hpc, dims = 1:40, reduction = "integrated.rpca")
DimPlot(neurons_hpc, reduction = "umap")
saveRDS(neurons_hpc, file = "neurons_hpc.rds")

#Find markers for each cluster----
neuron_marker_hpc <- FindAllMarkers(object = neurons_hpc, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
write_xlsx(neuron_marker_hpc,'neuron_subcluster_hpc12.xlsx')

n_cells <- table(neurons_hpc@meta.data$seurat_clusters, neurons_hpc@meta.data$sample)
write_xlsx(n_cells, "neuron_subcluster_hpc_n_cells12.xlsx")

#Rename clusters based on markers----
    #Markers were compared to gene expression data from the Allen Brain Cell Atlas to identify cell types
neurons_hpc <- RenameIdents(object = neurons_hpc,
                       '0' = "CA1_CA3_GLUT",
                       '1'= "DG_GLUT",
                       '2' = "GABA",
                       '3'= "NP_CT_L6b_GLUT")

DimPlot(neurons_hpc, reduction = "umap")
DimPlot(neurons_hpc, reduction = "umap", cols = RColorBrewer::brewer.pal(4, "Set2"))

neurons_hpc$cell_type_sub <- paste(Idents(neurons_hpc))

#Run DEG analysis----
neurons_hpc$cluster.treatment <- paste(Idents(neurons_hpc), neurons_hpc$treatment, sep = "_")
Idents(neurons_hpc) <- "cluster.treatment"

# Get unique cell types
cell_types_sub <- unique(neurons_hpc$cell_type_sub)

# Loop through each cell type GF vs CONV (positive lfc is up in GF)
for (cell_type in cell_types_sub) {
  # Define the groups for comparison
  group_CONV <- paste0(cell_type, "_CONV")
  group_GF <- paste0(cell_type, "_GF")
  
  # Check if both groups exist in the data set
  if (group_CONV %in% levels(Idents(neurons_hpc)) && group_GF %in% levels(Idents(neurons_hpc))) {
    # Run differential expression analysis
    DEG_df_name <- paste0("DEG_df_", cell_type)
    assign(DEG_df_name, FindMarkers(neurons_hpc,
                                        ident.1 = group_GF,
                                        ident.2 = group_CONV,
                                        test.use = "MAST",
                                        logfc.threshold = 0,
                                        only.pos = FALSE) %>%
             as.data.frame() %>%
             rownames_to_column(var = "geneID"))
    
    # Save DEG_df to a CSV file
    DEG_filename <- paste0(cell_type, "_GFvCONV_DEG.csv")
    write_csv(get(DEG_df_name), DEG_filename)
    } else {
    message(paste("Groups", group_CONV, "or", group_GF, "do not exist in the dataset for cell type", cell_type))
  }
}   

# Loop through each cell type 2wk vs GF (positive lfc is up in 2wk)
for (cell_type in cell_types_sub) {
  # Define the groups for comparison
  group_2wk <- paste0(cell_type, "_EC2wk")
  group_GF <- paste0(cell_type, "_GF")
  
  # Check if both groups exist in the data set
  if (group_2wk %in% levels(Idents(neurons_hpc)) && group_GF %in% levels(Idents(neurons_hpc))) {
    # Run differential expression analysis
    DEG_df_name_2wk <- paste0("DEG_df_2wk_", cell_type)
    assign(DEG_df_name_2wk, FindMarkers(neurons_hpc,
                                    ident.1 = group_2wk,
                                    ident.2 = group_GF,
                                    test.use = "MAST",
                                    logfc.threshold = 0,
                                    only.pos = FALSE) %>%
             as.data.frame() %>%
             rownames_to_column(var = "geneID"))
    
    # Save DEG_df to a CSV file
    DEG_filename <- paste0(cell_type, "_2wkvGF_DEG.csv")
    write_csv(get(DEG_df_name_2wk), DEG_filename)
    } else {
    message(paste("Groups", group_2wk, "or", group_GF, "do not exist in the dataset for cell type", cell_type))
  }
}

# Loop through each cell type 4wk vs GF (positive lfc is up in 4wk)
for (cell_type in cell_types_sub) {
  # Define the groups for comparison
  group_4wk <- paste0(cell_type, "_EC4wk")
  group_GF <- paste0(cell_type, "_GF")
  
  # Check if both groups exist in the data set
  if (group_4wk %in% levels(Idents(neurons_hpc)) && group_GF %in% levels(Idents(neurons_hpc))) {
    # Run differential expression analysis
    DEG_df_name_4wk <- paste0("DEG_df_4wk_", cell_type)
    assign(DEG_df_name_4wk, FindMarkers(neurons_hpc,
                                    ident.1 = group_4wk,
                                    ident.2 = group_GF,
                                    test.use = "MAST",
                                    logfc.threshold = 0,
                                    only.pos = FALSE) %>%
             as.data.frame() %>%
             rownames_to_column(var = "geneID"))
    
    # Save DEG_df to a CSV file
    DEG_filename <- paste0(cell_type, "_4wkvGF_DEG.csv")
    write_csv(get(DEG_df_name_4wk), DEG_filename)
    } else {
    message(paste("Groups", group_4wk, "or", group_GF, "do not exist in the dataset for cell type", cell_type))
  }
}