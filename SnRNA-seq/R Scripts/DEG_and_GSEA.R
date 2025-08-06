### DEG analysis and GSEA for each cell type ###

## This script performs differential expression analysis and Gene Set Enrichment Analysis (GSEA) for each cell type in the Seurat object.
## It compares GF vs CONV, 2wk vs GF, and 4wk vs GF conditions, saving results to CSV files.

## Overexpression-based pathway analysis was also performed using the online tool https://metascape.org/gp/index.html#/main/step1 once DEGs were identified. 
## Overlapping metascape pathway results were curated by hand from metascape results and graphed in R (not shown).

# Code was created by Lisa Blackmer-Raynolds

#Load R libraries
library("Seurat")
library("ggplot2")
library("cowplot")
library("dplyr")
library("clusterProfiler")
library("org.Mm.eg.db")
library("msigdbr")
library("DOSE")
library("stringr")
library("enrichplot")
library("tidyverse")
library("tibble")
library("readr")
library("biomaRt")
library("rtracklayer")

#Set up seurat object for loops:
hpc@meta.data$cell_type_treatment <- paste(hpc@meta.data$cell_type, 
                                           hpc@meta.data$treatment, sep = "_")
Idents(hpc) <- "cell_type_treatment"

# Get unique cell types
cell_types <- unique(hpc$cell_type)

#Set up mart for mapping
ensembl <- useEnsembl(
  biomart = "ensembl",
  dataset = "mmusculus_gene_ensembl",
  mirror = "useast"  # Change to "asia" or "useast" if needed
)

# Loop through each cell type GF vs CONV (positive lfc is up in GF)----
for (cell_type in cell_types) {
  # Define the groups for comparison
  group_CONV <- paste0(cell_type, "_CONV")
  group_GF <- paste0(cell_type, "_GF")
  
  # Check if both groups exist in the data set
  if (group_CONV %in% levels(Idents(hpc)) && group_GF %in% levels(Idents(hpc))) {
    # Run differential expression analysis
    DEG_df_name <- paste0("DEG_df_", cell_type)
    assign(DEG_df_name, FindMarkers(hpc,
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
    
    # Prepare gene list for GSEA
    gene_list_name <- paste0("gene_list_", cell_type)
    gene_list <- get(DEG_df_name) %>%
      dplyr::select(geneID, avg_log2FC, p_val_adj) %>%
      mutate(weighted_score = avg_log2FC * -log10(p_val_adj + 1e-300))
    
    assign(gene_list_name, gene_list)
    
    # Create gene_list_gsea
    gene_list_gsea <- gene_list$avg_log2FC
    names(gene_list_gsea) <- gene_list$geneID
    gene_list_gsea <- sort(gene_list_gsea, decreasing = TRUE)
    
    gene_list_gsea_name <- paste0("gene_list_gsea_", cell_type)
    assign(gene_list_gsea_name, gene_list_gsea)
    
    # Run GSEA
    gsea_result_name <- paste0("gsea_result_", cell_type)
    assign(gsea_result_name, gseGO(geneList = gene_list_gsea,
                                   OrgDb = org.Mm.eg.db,
                                   ont = "BP",
                                   keyType = "SYMBOL",
                                   pAdjustMethod = "BH",
                                   pvalueCutoff = 0.05,
                                   eps = 0,
                                   nPermSimple = 100000))
    
    # Save GSEA results to a CSV file
    result_df_name <- paste0("result_df_", cell_type)
    assign(result_df_name, get(gsea_result_name)@result)
    output_filename <- paste0(cell_type, "_GFvCONV_GO.csv")
    write_csv(get(result_df_name), output_filename)
    
    # Map gene symbols to Entrez IDs
    gene_mapping_ensembl <- getBM(
      attributes = c("external_gene_name", "entrezgene_id"),
      filters = "external_gene_name",
      values = names(gene_list_gsea),
      mart = ensembl)
    
    # Merge the scores from the original gene list with Entrez IDs
    gene_list_kegg <- gene_list_gsea[gene_mapping_ensembl$external_gene_name]
    names(gene_list_kegg) <- gene_mapping_ensembl$entrezgene_id
    
    # Remove NA values
    gene_list_kegg <- gene_list_kegg[!is.na(names(gene_list_kegg)) & !is.na(gene_list_kegg)]
    
    # Average the value for duplicate Entrez IDs
    gene_df <- data.frame(
      entrez_id = names(gene_list_kegg),
      score = gene_list_kegg
    ) %>%
      group_by(entrez_id) %>%
      summarize(score = mean(score), .groups = "drop")
    
    # Convert back to named numeric vector
    gene_list_kegg <- setNames(gene_df$score, gene_df$entrez_id)
    
    # Ensure it is sorted in decreasing order
    gene_list_kegg <- sort(gene_list_kegg, decreasing = TRUE)
    
    # Run KEGG GSEA
    kegg_result_name <- paste0("kegg_result_", cell_type)
    assign(kegg_result_name, gseKEGG(
      geneList = gene_list_kegg,
      organism = 'mmu',
      keyType = "ncbi-geneid",
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      eps = 0
    ))
    
    # Save KEGG results to a CSV file
    kegg_result_df_name <- paste0("kegg_result_df_", cell_type)
    assign(kegg_result_df_name, get(kegg_result_name)@result)
    kegg_output_filename <- paste0(cell_type, "_GFvCONV_KEGG.csv")
    write_csv(get(kegg_result_df_name), kegg_output_filename)
    
    # Print progress
    cat("GSEA and KEGG analysis completed for:", cell_type, "\n")
  } else {
    cat("Skipped cell type:", cell_type, "- One or both groups not found.\n")
  }
}

# Loop through each cell type 2wk vs GF (positive lfc is up in 2wk)----
for (cell_type in cell_types) {
  # Define the groups for comparison
  group_2wk <- paste0(cell_type, "_EC2wk")
  group_GF <- paste0(cell_type, "_GF")
  
  # Check if both groups exist in the data set
  if (group_2wk %in% levels(Idents(hpc)) && group_GF %in% levels(Idents(hpc))) {
    # Run differential expression analysis
    DEG_df_name <- paste0("DEG_df_", cell_type)
    assign(DEG_df_name, FindMarkers(hpc,
                                    ident.1 = group_2wk,
                                    ident.2 = group_GF,
                                    test.use = "MAST",
                                    logfc.threshold = 0,
                                    only.pos = FALSE) %>%
             as.data.frame() %>%
             rownames_to_column(var = "geneID"))
    
    # Save DEG_df to a CSV file
    deg_filename <- paste0(cell_type, "_2wkvGF_DEG.csv")
    write_csv(get(DEG_df_name), deg_filename)
    
    # Prepare gene list for GSEA
    gene_list_name <- paste0("gene_list_", cell_type)
    gene_list <- get(DEG_df_name) %>%
      dplyr::select(geneID, avg_log2FC, p_val_adj) %>%
      mutate(weighted_score = avg_log2FC * -log10(p_val_adj + 1e-300))
    
    assign(gene_list_name, gene_list)
    
    # Create gene_list_gsea
    gene_list_gsea <- gene_list$avg_log2FC
    names(gene_list_gsea) <- gene_list$geneID
    gene_list_gsea <- sort(gene_list_gsea, decreasing = TRUE)
    
    gene_list_gsea_name <- paste0("gene_list_gsea_", cell_type)
    assign(gene_list_gsea_name, gene_list_gsea)
    
    # Run GSEA
    gsea_result_name <- paste0("gsea_result_", cell_type)
    assign(gsea_result_name, gseGO(geneList = gene_list_gsea,
                                   OrgDb = org.Mm.eg.db,
                                   ont = "BP",
                                   keyType = "SYMBOL",
                                   pAdjustMethod = "BH",
                                   pvalueCutoff = 0.05,
                                   eps = 0,
                                   nPermSimple = 100000))
    
    # Save GSEA results to a CSV file
    result_df_name <- paste0("result_df_", cell_type)
    assign(result_df_name, get(gsea_result_name)@result)
    output_filename <- paste0(cell_type, "_GFv2wk_GO.csv")
    write_csv(get(result_df_name), output_filename)
    
    # Map gene symbols to Entrez IDs
    gene_mapping_ensembl <- getBM(
      attributes = c("external_gene_name", "entrezgene_id"),
      filters = "external_gene_name",
      values = names(gene_list_gsea),
      mart = ensembl)
    
    # Merge the scores from the original gene list with Entrez IDs
    gene_list_kegg <- gene_list_gsea[gene_mapping_ensembl$external_gene_name]
    names(gene_list_kegg) <- gene_mapping_ensembl$entrezgene_id
    
    # Remove NA values
    gene_list_kegg <- gene_list_kegg[!is.na(names(gene_list_kegg))]
    
    # Average the value for duplicate Entrez IDs
    gene_df <- data.frame(
      entrez_id = names(gene_list_kegg),
      score = gene_list_kegg
    ) %>%
      group_by(entrez_id) %>%
      summarize(score = mean(score), .groups = "drop")
    
    # Convert back to named numeric vector
    gene_list_kegg <- setNames(gene_df$score, gene_df$entrez_id)
    
    # Ensure it is sorted in decreasing order
    gene_list_kegg <- sort(gene_list_kegg, decreasing = TRUE)
    
    # Run KEGG GSEA
    kegg_result_name <- paste0("kegg_result_", cell_type)
    assign(kegg_result_name, gseKEGG(
      geneList = gene_list_kegg,
      organism = 'mmu',
      keyType = "ncbi-geneid",
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      eps = 0
    ))
    
    # Save KEGG results to a CSV file
    kegg_result_df_name <- paste0("kegg_result_df_", cell_type)
    assign(kegg_result_df_name, get(kegg_result_name)@result)
    kegg_output_filename <- paste0(cell_type, "_GFv2wk_KEGG.csv")
    write_csv(get(kegg_result_df_name), kegg_output_filename)
    
    # Print progress
    cat("GSEA and KEGG analysis completed for:", cell_type, "\n")
  } else {
    cat("Skipped cell type:", cell_type, "- One or both groups not found.\n")
  }
}

# Loop through each cell type 4wk vs GF (positive lfc is up in 4wk)----
for (cell_type in cell_types) {
  # Define the groups for comparison
  group_4wk <- paste0(cell_type, "_EC4wk")
  group_GF <- paste0(cell_type, "_GF")
  
  # Check if both groups exist in the data set
  if (group_4wk %in% levels(Idents(hpc)) && group_GF %in% levels(Idents(hpc))) {
    # Run differential expression analysis
    DEG_df_name <- paste0("DEG_df_", cell_type)
    assign(DEG_df_name, FindMarkers(hpc,
                                    ident.1 = group_4wk,
                                    ident.2 = group_GF,
                                    test.use = "MAST",
                                    logfc.threshold = 0,
                                    only.pos = FALSE) %>%
             as.data.frame() %>%
             rownames_to_column(var = "geneID"))
    
    # Save DEG_df to a CSV file
    deg_filename <- paste0(cell_type, "_4wkvGF_DEG.csv")
    write_csv(get(DEG_df_name), deg_filename)
    
    # Prepare gene list for GSEA
    gene_list_name <- paste0("gene_list_", cell_type)
    gene_list <- get(DEG_df_name) %>%
      dplyr::select(geneID, avg_log2FC, p_val_adj) %>%
      mutate(weighted_score = avg_log2FC * -log10(p_val_adj + 1e-300))
    
    assign(gene_list_name, gene_list)
    
    # Create gene_list_gsea
    gene_list_gsea <- gene_list$avg_log2FC
    names(gene_list_gsea) <- gene_list$geneID
    gene_list_gsea <- sort(gene_list_gsea, decreasing = TRUE)
    
    gene_list_gsea_name <- paste0("gene_list_gsea_", cell_type)
    assign(gene_list_gsea_name, gene_list_gsea)
    
    # Run GSEA
    gsea_result_name <- paste0("gsea_result_", cell_type)
    assign(gsea_result_name, gseGO(geneList = gene_list_gsea,
                                   OrgDb = org.Mm.eg.db,
                                   ont = "BP",
                                   keyType = "SYMBOL",
                                   pAdjustMethod = "BH",
                                   pvalueCutoff = 0.05,
                                   eps = 0,
                                   nPermSimple = 100000))
    
    # Save GSEA results to a CSV file
    result_df_name <- paste0("result_df_", cell_type)
    assign(result_df_name, get(gsea_result_name)@result)
    output_filename <- paste0(cell_type, "_GFv4wk_GO.csv")
    write_csv(get(result_df_name), output_filename)
    
    # Map gene symbols to Entrez IDs
    gene_mapping_ensembl <- getBM(
      attributes = c("external_gene_name", "entrezgene_id"),
      filters = "external_gene_name",
      values = names(gene_list_gsea),
      mart = ensembl)
    
    # Merge the scores from the original gene list with Entrez IDs
    gene_list_kegg <- gene_list_gsea[gene_mapping_ensembl$external_gene_name]
    names(gene_list_kegg) <- gene_mapping_ensembl$entrezgene_id
    
    # Remove NA values
    gene_list_kegg <- gene_list_kegg[!is.na(names(gene_list_kegg))]
    
    # Average the value for duplicate Entrez IDs
    gene_df <- data.frame(
      entrez_id = names(gene_list_kegg),
      score = gene_list_kegg
    ) %>%
      group_by(entrez_id) %>%
      summarize(score = mean(score), .groups = "drop")
    
    # Convert back to named numeric vector
    gene_list_kegg <- setNames(gene_df$score, gene_df$entrez_id)
    
    # Ensure it is sorted in decreasing order
    gene_list_kegg <- sort(gene_list_kegg, decreasing = TRUE)
    
    # Run KEGG GSEA
    kegg_result_name <- paste0("kegg_result_", cell_type)
    assign(kegg_result_name, gseKEGG(
      geneList = gene_list_kegg,
      organism = 'mmu',
      keyType = "ncbi-geneid",
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      eps = 0
    ))
    
    # Save KEGG results to a CSV file
    kegg_result_df_name <- paste0("kegg_result_df_", cell_type)
    assign(kegg_result_df_name, get(kegg_result_name)@result)
    kegg_output_filename <- paste0(cell_type, "_4wkvGF_KEGG.csv")
    write_csv(get(kegg_result_df_name), kegg_output_filename)
    
    # Print progress
    cat("GSEA and KEGG analysis completed for:", cell_type, "\n")
  } else {
    cat("Skipped cell type:", cell_type, "- One or both groups not found.\n")
  }
}

