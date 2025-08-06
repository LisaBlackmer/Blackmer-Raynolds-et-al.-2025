### Example Bulk RNA-seq DEG Analysis ###

#This is an example script for performing normalization, filtering, differential expression analysis, and gene set enrichment analysis on bulk RNA-seq data.
#In this example, germ-free (GF) mice were compared to E. coli (EC) mono-colonized mice, however, similar steps were used to make other pairwise comparisions in both mono-colonized and 5xFAD cohorts.
#Each pairwise comparison was normalized, filtered, and analyzed independently using these steps.

#Code is based on the following tutorial by Dr. Daniel Beiting at the University of Pennsylvania:
#https://diytranscriptomics.com/ and modified by Lisa Blackmer-Raynolds with imput and troubleshooting from Maureen Sampson


# Load packages----
library(tidyverse)
library(tximport)
library(edgeR)
library(matrixStats)
library(cowplot)
library(plotly)
library(readxl)
library(limma)
library(ggrepel)
library(fgsea)
library(dplyr)

# Importing and annotating----

# Load the sample information
targets <- read_tsv("GFvsEC_IDs.txt")

# Create a vector of sample labels
sampleLabels <- targets$Mouse_ID

# Create a path to the abundance files
path <- file.path(targets$Lib_ID, "abundance.tsv")

# Check that the files exist
all(file.exists(path)) 

# Import the transcript-to-gene mapping file
Tx.mouse <- read.delim("Gene_to_transcript.txt")

# Import Kallisto transcript counts into R using Tximport
Txi_gene <- tximport(path, #points to the abundance files (made above)
                     type = "kallisto", #type of aliment data you are reading in
                     tx2gene = Tx.mouse, #this points to my annotation file
                     txOut = FALSE, #this tells it to collapse it to gene level (rather than transcript) data
                     countsFromAbundance = "lengthScaledTPM", #how we will calculate counts from abund
                     ignoreTxVersion = TRUE) #this makes it ignore the transcript version

# Check the data and generate summary statistics----

myTPM <- Txi_gene$abundance # Save the TPM values
myCounts <- Txi_gene$counts # Save the raw counts
colSums(myTPM) # Check the total TPM values for each sample
colSums(myCounts) # Check the total counts for each sample

myCounts.df <- as_tibble(myCounts, rownames = "geneID") # Convert counts to a tibble
colnames(myCounts.df) <- c("geneID", sampleLabels) # Add sample labels to the counts dataframe

myTPM.df <- as_tibble(myTPM, rownames = "geneID") # Convert to a tibble

# Calculate summary statistics for the TPM data
myTPM.stats <- transform(myTPM, 
                         SD=rowSds(myTPM), 
                         AVG=rowMeans(myTPM),
                         MED=rowMedians(myTPM))

# Make a DGElist from counts
myDGEList <- DGEList(myCounts)
save(myDGEList, file = "GFvsEC_mixed_sex_DGEList")

# Filter and normalize the data ----

# Use the 'cpm' function from EdgeR to get log2 counts per million
cpm <- cpm(myDGEList) 
colSums(cpm)
log2.cpm <- cpm(myDGEList, log=TRUE)

# Coerce the data matrix to a dataframe
log2.cpm.df <- as_tibble(log2.cpm, rownames = "geneID")

# Add sample names to this dataframe
colnames(log2.cpm.df) <- c("geneID", sampleLabels)

# Pivot the dataframe (from wide to long)
log2.cpm.df.pivot <- pivot_longer(log2.cpm.df,
                                  cols = GFm43:ECf33, # column names to be stored as a SINGLE variable **Change this to match sample names**
                                  names_to = "samples", # name of that new variable (column)
                                  values_to = "expression") # name of new variable (column) storing all the values (data)

log2.cpm.df.pivot

# Visualize the data
ggplot(log2.cpm.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="unfiltered, non-normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

# Filter the data to remove lowly expressed genes.
# Here we are keeping genes that have at least 1 counts per million in at least 3 samples. 
# ** The number of samples should be adjusted to reflect the smallest sample size in a single treatment group.**

keepers <- rowSums(cpm>1)>=3 # This creates a logical vector indicating which genes to keep
myDGEList.filtered <- myDGEList[keepers,] # This subsets the DGEList to keep only the genes that meet the criteria
dim(myDGEList.filtered)

# Recalculate log2 CPM for the filtered data
log2.cpm.filtered <- cpm(myDGEList.filtered, log=TRUE)
log2.cpm.filtered.df <- as_tibble(log2.cpm.filtered, rownames = "geneID")
colnames(log2.cpm.filtered.df) <- c("geneID", sampleLabels)

# Pivot and graph the filtered data
log2.cpm.filtered.df.pivot <- pivot_longer(log2.cpm.filtered.df,
                                           cols = GFm43:ECf33, # column names to be stored as a SINGLE variable **Change this to match sample names**
                                           names_to = "samples",
                                           values_to = "expression")

ggplot(log2.cpm.filtered.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="filtered, non-normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

# TMM normalize the data
myDGEList.filtered.norm <- calcNormFactors(myDGEList.filtered, method = "TMM")

# Calculate log2 CPM for the normalized data
log2.cpm.filtered.norm <- cpm(myDGEList.filtered.norm, log=TRUE)
log2.cpm.filtered.norm.df <- as_tibble(log2.cpm.filtered.norm, rownames = "geneID")
colnames(log2.cpm.filtered.norm.df) <- c("geneID", sampleLabels)

# Pivot and graph the normalized data
log2.cpm.filtered.norm.df.pivot <- pivot_longer(log2.cpm.filtered.norm.df,
                                                cols = GFm43:ECf33, # column names to be stored as a SINGLE variable **Change this to match sample names**
                                                names_to = "samples",
                                                values_to = "expression")

ggplot(log2.cpm.filtered.norm.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="filtered, TMM normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

# Save the normalized data
write_tsv(log2.cpm.filtered.norm.df, "normalized log2 cpm.txt")

# PCA plots -----------
# Identify variables of interest in study design file
targets
group <- targets$Treatment
group <- factor(group)

sex <- targets$Sex
sex <- factor(sex)

# Hierarchical clustering
distance <- dist(t(log2.cpm.filtered.norm), method = "maximum")
clusters <- hclust(distance, method = "average") 
plot(clusters, labels=sampleLabels)

# Principal component analysis (PCA)
pca.res <- prcomp(t(log2.cpm.filtered.norm), scale.=F, retx=T)
ls(pca.res)
summary(pca.res) # Prints variance summary for all principal components.
pca.res$rotation #$rotation shows you how much each gene influenced each PC (called 'scores')
pca.res$x # 'x' shows you how much each sample influenced each PC (called 'loadings')
#note that these have a magnitude and a direction (this is the basis for making a PCA plot)
screeplot(pca.res) # A screeplot is a standard way to view eigenvalues for each PCA
pc.var<-pca.res$sdev^2 # sdev^2 captures these eigenvalues from the PCA result
pc.per<-round(pc.var/sum(pc.var)*100, 1) # we can then use these eigenvalues to calculate the percentage variance explained by each PC
pc.per

# Visualize your PCA result
pca.res.df <- as_tibble(pca.res$x)

PCA <-ggplot(pca.res.df) +
  aes(x=PC1, y=PC2, label=sampleLabels, color=group, shape=sex) +
  geom_point(size=4) +
  # geom_label() + #If you want to label the points with sample labels
  # stat_ellipse() + #If you want to circle clusters on the PCA
  # coord_fixed() + #If you want to apply the correct aspect ratio
  xlab(paste0("PC1 (",pc.per[1],"%",")")) + 
  ylab(paste0("PC2 (",pc.per[2],"%",")")) +
  labs(title="PCA plot",
       caption=paste0("produced on ", Sys.time())) +
  coord_fixed() +
  theme_bw()

ggplotly(PCA)

# Create a PCA 'small multiples' chart
pca.res.df <- pca.res$x[,1:4] %>%
  as_tibble() %>%
  add_column(sample = sampleLabels,
             group = group,
             sex = sex)

pca.pivot <- pivot_longer(pca.res.df,
                          cols = PC1:PC4,
                          names_to = "PC",
                          values_to = "loadings")

ggplot(pca.pivot) +
  aes(x=sample, y=loadings, fill=group) +
  geom_bar(stat="identity") +
  facet_wrap(~PC) +
  labs(title="PCA 'small multiples' plot GF vs E. coli mixed sex",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw() +
  coord_flip()

#Add or remove alpha=sex under aes to see sex and treatment or change fill=group to fill=sex to see more clearly

# Add a column with the average expression of each group
mydata.df <- log2.cpm.filtered.norm.df %>% 
  mutate(GF.AVG = (GFm43 + GFf44 + GFf45 + GFf66 + GFf67 + GFm68 + GFm69)/7, # Note these need to be changed to match your sample names
         EC.AVG = (ECm27 + ECm28 + ECm29 + ECm30 + ECf31 + ECf32 + ECf33)/7, # Note these need to be changed to match your sample names
         LogFC = (EC.AVG - GF.AVG)) %>% #with this order genes that are higher in my EC group will have a positive number
  mutate_if(is.numeric, round, 2)

# Sort the data by LogFC
mydata.sort <- mydata.df %>%
  dplyr::arrange(desc(LogFC)) %>% 
  dplyr::select(geneID, LogFC)


#DEG Analysis-----

# Set up design matrix
group <- factor(targets$Treatment)
design <- model.matrix(~0 + group) #start at 0 and made comparisons based on group (so in this case I'll have 0 and 1)
colnames(design) <- levels(group) #renames column names

# Use VOOM function from Limma package to model the mean-variance relationship
v.DEGList.filtered.norm <- voom(myDGEList.filtered.norm, design, plot = TRUE)
# Fit a linear model to the data
fit <- lmFit(v.DEGList.filtered.norm, design)

# Make a contrast matrix
contrast.matrix <- makeContrasts(colonization = EC - GF,
                                 #sex = M - F,
                                 levels=design)

# Extract the linear model fit (get bayesian stats)
fits <- contrasts.fit(fit, contrast.matrix)
ebFit <- eBayes(fits)
write.fit(ebFit, file="lmfit_results_GFvsEC_mixed.txt")

# View DEGs
myTopHits <- topTable(ebFit, adjust ="BH", coef=1, number=10000000000, sort.by="logFC")

# Convert to a tibble
myTopHits.df <- myTopHits %>%
  as_tibble(rownames = "geneID")

# Save
write_tsv(myTopHits.df,"DEGs_GFvsEC_mixed_sex.txt")

#Filter data to remove genes that aren't enriched in microglia----

#Load the microglia enrichment data (from Zhang et al. 2014 PMID: 25186741): https://brainrnaseq.org/mouse-microglia-bennett-et-al-2016/?1240214697=1539028313
Brain_RNAseq_Data <- read_tsv("Brain_RNAseq_Data.txt")

# Filter the microglia enrichment data to only include genes that are enriched in microglia.
# MG_Enrichment is calculated as (average expression in microglia cells) / (average expression all non-microglia brain cells) therefore, values >1.1 indicate genes more highly expressed in microglia
Brain_RNAseq_Data_filtered <- dplyr::filter(Brain_RNAseq_Data, MG_Enrichment > 1.1)

# To make sure all gene symbols are the same format, we convert them both to uppercase
Brain_RNAseq_Data_filtered$gene.symbol = toupper(Brain_RNAseq_Data_filtered$gene.symbol)
myTopHits.df$geneID = toupper(myTopHits.df$geneID)

# Filter the DEGs to only include genes that are enriched in microglia
myTopHits_filtered <- subset(myTopHits.df, geneID %in% Brain_RNAseq_Data_filtered$gene.symbol)

# Save the filtered DEGs
write_tsv(myTopHits_filtered,"DEGs_MG_filtered_GFvsEC_mixed_sex.txt")

# Volcano plots----
vplot <- ggplot(myTopHits_filtered) +
  aes(y=-log10(adj.P.Val), x=logFC, text = paste("Symbol:", geneID)) +
  geom_point(size=1) +
  geom_point(data = myTopHits_filtered %>% filter(adj.P.Val < 0.05), color = "deeppink1") +
  geom_point(data = myTopHits_filtered %>% filter(adj.P.Val >= 0.05), color = "darkslateblue") +
  #geom_label_repel(data = myTopHits.select.df,
  #aes(label = geneID),
  #max.overlaps = Inf,
  #size = 2, color = "black",
  #box.padding = unit(0.35, "lines"),
  #point.padding = unit(0.3, "lines"),
  #min.segment.length = unit(0, 'lines'))+
  geom_label_repel(data = subset(myTopHits_filtered, adj.P.Val<0.05),
                   aes(label = geneID),
                   max.overlaps = 40,
                   size = 2, color = "black",
                   min.segment.length = unit(0, 'lines'))+
  ylim(0,3.02) +
  #box.padding = unit(0.35, "lines"),
  #point.padding = unit(0.3, "lines"))+
  geom_hline(yintercept = -log10(0.05), linetype="longdash", colour="black", size=0.5) +
  #geom_vline(xintercept = 1, linetype="longdash", colour="#BE684D", size=1) +
  #geom_vline(xintercept = -1, linetype="longdash", colour="#2C467A", size=1) +
  #annotate("rect", xmin = 1, xmax = 12, ymin = -log10(0.01), ymax = 7.5, alpha=.2, fill="#BE684D") +
  #annotate("rect", xmin = -1, xmax = -12, ymin = -log10(0.01), ymax = 7.5, alpha=.2, fill="#2C467A") +
  #geom_point(data = myTopHits.df %>% filter(geneID == "LDHA"), size=1, color = "purple2") +
  #labs(title="Volcano plot",
  #    subtitle = "Descriptive Title",
  #   caption=paste0("produced on ", Sys.time())) +
    theme_bw()
  
  ggplotly(vplot)

#Run GSEA analysis----

#Set up data for GSEA
mydata.df.sub <- dplyr::select(myTopHits.df, geneID, logFC, adj.P.Val)
mydata.df.sub <- mydata.df.sub %>%
  mutate(weighted_score = logFC * -log10(adj.P.Val))
mydata.gsea <- mydata.df.sub$weighted_score
names(mydata.gsea) <- as.character(mydata.df.sub$geneID)
mydata.gsea <- sort(mydata.gsea, decreasing = TRUE)

#Load gene list
LPS_GO.df <- read_excel("LPS GO.xlsx")

# Convert each column to a gene vector, removing NAs, and normalizing capitalization
LPS_gene_sets <- lapply(LPS_GO.df, function(x) na.omit(x))

LPS_gene_sets <- lapply(LPS_gene_sets, toupper)

#Run GSEA using fgsea
fgsea_results <- fgsea(
  pathways = LPS_gene_sets,
  stats = mydata.gsea)

print(fgsea_results)

#Make an enrichment tabel
plotGseaTable(
  pathways = LPS_gene_sets,
  stats = mydata.gsea,
  fgseaRes = fgsea_results,
  gseaParam = 0.5)



  
  
  
