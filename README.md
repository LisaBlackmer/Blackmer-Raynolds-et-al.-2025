# Blackmer-Raynolds et al., 2025 Code Repository

This GitHub repository contains the code for the paper titled: "Indigenous gut microbes modulate neural cell state and neurodegenerative disease
susceptibility." For a detailed description of the methods and results, please refer to the [paper](#). *Note: The paper is currently under review and the links will be updated once it is published.*

While this code is meant to offer insights into our study methods and approach, it is not expected to function directly out-of-the-box for reproduction of our results. Should you have any questions or want further clarification, feel free to reach out.

## Data Availability:
The data used in this study is available through the Gene Expression Omnibus (GEO) database:
- **SnRNA-seq data:** [GSE289589](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289589)
- **Mono-colonized bulk RNA-seq data:** [GSE289591](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289591)
- **5xFAD bulk RNA-seq data:** [GSE289590](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289590)

---

## SnRNA-seq Data Processing:

1. **Alignment:**  
    Alignment to the mouse reference genome (GRCm39) was performed using the Cell Ranger pipeline. Ambient RNA was removed using Cellbender with the shell script `Cellranger_cellbender.sh`. The dependencies for this step are listed in `Cellranger_cellbender_dependencies.txt`.

2. **R Packages:**  
    R packages can be found and installed using the renv package manager with the files located in the `renv` folder located within the associated `R Scripts` folder.

3. **Data Integration & Clustering:**  
    The data was imported into R, integrated, and clustered with the script `Make_seurat_object.R`.

4. **Differential Expression & Enrichment Analysis:**  
    Differential expression analysis and gene set enrichment analysis were performed with the script `DEG_and_GSEA.R`.

5. **Subclustering and Analysis:**  
    - Neurons: `Neuronal_subclustering.R`
    - Immune cells: `Immune_subclustering.R`

---

## Bulk RNA-seq Data Processing (Mono-Colonized & 5xFAD):

1. **Alignment:**  
    Alignment to the mouse reference genome (GRCm39) was performed using Kallisto with the script `kallisto.sh`. The dependencies for this step are listed in `kallisto_dependencies.txt`.

2. **R Packages:**  
    R packages can be found and installed using the renv package manager with the files located in the `renv` folder.

3. **Data Analysis:**  
    The data was imported into R, normalized, and analyzed using the steps outlined in `Example_DEG_analysis.R`. *Note: This script is an example; however all pairwise comparisons were performed using the same method.*

---

## Acknowledgments:

The contributors and resources used to make each script are detailed in the header comments of the respective script files. Packages/programs are listed in the `dependencies.txt` files (for alignment scripts) and in the `renv.lock` file (for R scripts) and cited in the paper.

If you use this code, please cite the paper as follows:
> Citation


