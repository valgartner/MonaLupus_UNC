# Overview
This repository contains code for analyzing RNA-Sequencing data for the **M**ulti-**O**mic i**N**tegrated **A**nalysis in Lupus (**MONA Lupus**) project at University of North Carolina Chapel Hill.

# Data
**Note on batch in this dataset**
These samples were sequenced across 2 batches, and unfortunately most of the Lupus samples were in batch 1 and all of the controls were in batch 2. A handful of Lupus samples were also sequenced in batch 2, so these samples can at least be compared using PCA. In the plot (**Fig 1**) below, samples in "U" (unknown) batch set were split across batch 1 and 2. These "U" samples cluster by disease vs control, not by batch, which gives us some confidence that batch effects are not the primary driver of variation.

[Fig 1.](results/MonaLupus_BatchEffectsPCA.png)
**Fig 1. Batch effect PCA.** *This plot was generated in `20253010_MonaLupus_DifferentialExpression.Rmd`*

----
# Directory structure
**counts/**
- RNA-Sequencing counts file (output from RNA-Snakemake pipeline): `Fureylab_Lupus_20221220.counts.txt`

**metadata/**
- coldata file: `mona-lupus_colfile.txt`
- original metadata excel sheet from Shruti: `MONA Lupus Complete De-Identified Dataset.xlsx`
	- each tab of this table was saved as a text file (and manually edited for better readability in R)
	- separate files were joined into `mona-lupus_colfile.txt` using the script `Organize_metadata.Rmd` also located in this folder

**analysis/**
*See "Analyses" section below for a description of each script*
- `20253010_MonaLupus_DifferentialExpression.Rmd`
- `webgestaltR.Rmd`
- `upset-plots.Rmd`
- `query_DEG_results.Rmd`
- `DEG_Count.Rmd`

**results/**
- Plots and files not specific to differential expression output (e.g. data quality, PCA plots, etc.)
- Output for each DE contrast saved in its own directory:
	- `stable_disease-All_ESRD/`
	- `control-All_ESRD/`
	- `control-stable_disease/`
	- `control-disease/`
- output from `webgestaltR.Rmd` saved in this directory:
	- `gsea/`

----
# Analyses
## 1. Differential Expression Analysis
**Script:** `analyses/20253010_MonaLupus_DifferentialExpression.Rmd`
- This script takes two inputs:
	- metadata: `mona-lupus_colfile.txt`
	- RNA-Sequencing counts data: `Fureylab_Lupus_20221220.counts.txt`
- Analyses performed:
	- **Differential Expression**
		- **Part 1**  - read in data, filter out samples that are low quality or unresolved outliers, make sure coldata & counts dataframes have the same samples included, find x&y chromosome genes for removal from the negative control gene set
		- **Part 2** - perform the initial differential expression analysis with only sex + group in the design (no correction for unwanted variation)
		- **Part 3** - correct for unwanted variation using RUVg
			- Find negative control genes - not differential between groups. In this dataset, it was difficult to find a "true" negative control dataset. We would hope to see negative control genes from all samples well mixed & no separation between conditions (black and red correspond to disease and control) . 
			- Next, test out multiple values of k to determine how many RUV factors to remove.
				- We plot metadata variables correlated with RUV factors to understand possible sources of variation besides disease. In this dataset, we see that both TIN score (quality score) and disease status correlated with RUV factor one. 
				- While we want to remove variation due to TIN score, we ideally wouldn't remove variation related to disease status! Because of the way the batches were run, we unfortunately can't separate these two factors. K values > 1 don't correspond with any metadata variables. 
				- **We decided to use k=1 to remove only the first factor of variation using RUVg to remove the variation due to TIN score.**
				- This code chunk also plots how many differentially expressed genes (DEGs) there are when DESeq is run after removing *k* factors of variation
		- **Part 4** - Finalize correction using RUVg and visualize on PCA
			- with and without control group
		- **Part 5** - Differential expression for different disease subtype contrasts. Note that before running `DESeq()` using the Wald test, I filter each dataframe to include only the samples for that particular contrast (e.g. see lines 524-540)
			- DE analysis A: by "group" variable (control vs disease)
			- DE analysis B: by "condition" variable (control vs stable disease)
			- DE analysis C: by "ESRD" variable (control vs All_ESRD, and stable disease vs All_ESRD)
			- Code chunks in part 5 will output files into contrast subdirectories:
				- FULL DESeq results
				- SIGNIFICANT DESeq results
				- volcano plots illustrating DEG distribution
				- files formatted for use in enrichr and gsea pathway analysis tools

	- **Other plots generated in this script**
		- Correlation between sclerotic glomeruli and TIN score (no correlation seen!)
		- TIN median by condition
		- sample sizes of each contrast (for slides)
		- expression heatmap
		- PCA by batch (**Fig 1 in this document**)
		- Interactive PCA plot to make it easier to identify outliers
		- Plot gene counts for each disease group
			- *Useful for discussions & Shehzad's favorite plots!*


## 2. Pathway Analysis  and upset plots with Web Gestalt
**Script:** `analyses/webgestaltR.Rmd`
- This script takes as input all the `*.rnk` files generated by script 1 (found in `/results/analysis_date/CONTRAST/gsea/`) and performs pathway analysis with the `WebGestaltR` package.
- For each contrast, run `WebGestaltR` function on each gene list and outputs the full pathway analysis results. The output files to care about are:
	- `enriched_geneset_wsc_topsets_*_gsea.txt` which contains the top 20 pathways using weighted set cover to reduce redundancy.
	- `enrichment_results__gsea.txt` which contains the full output (including many similar pathways)

-----
## 3. Visualize shared genes and pathways using upset plots
**Script:** `analyses/upset-plots.Rmd`
- This script takes several inputs: 
	- significant DEG results for each contrast e.g. `results/analysis_date/CONTRAST/analysis_date_Fureylab_Mona-Lupus_DESeq_SIGIFICANT_results_CONTRASTpadj0.05_kval1.txt`
	- the table containing the top 20 pathways for each contrast generated by `webgestaltR.Rmd` saved here:
		- `/results/analysis_date/gsea/all_gsea_weighted-set-cover_df.txt`
- This script uses the package `ComplexUpset` to make upset plots. I highly recommend this package over `UpsetR`, as `ComplexUpset` makes it very easy to pull out members of each intersection, and also I found it easier to get the data into the correct shape with this package compared to `UpsetR`. 
- My original code using `UpsetR` is included but commented out at the end of this script. 

**Pulling out intersection members**
https://github.com/krassowski/complex-upset/issues/125
exclusive intersections:
```r
upset_data(movies, genres)$with_sizes %>%
    filter(intersection == "Comedy-Drama" & in_exclusive_intersection)
```

inclusive intersections:
```r
upset_data(movies, genres)$with_sizes %>%
    filter(intersection == "Comedy-Drama" & in_inclusive_intersection)
```

*Understanding inclusive vs exclusive intersections:*
https://krassowski.github.io/complex-upset/articles/Examples_R.html#region-selection-modes


------
## Query DE results
**Script:** `analyses/DEG_Count.Rmd`
- This script pulls out specific lists of genes from each contrast so it's easy to see if your gene of interest is significantly differentially expressed or not.

------
## DEG Count
**Script:** `analyses/query_DEG_results.Rmd`
- Short script to summarize the number of DEGs in each contrast and whether they're protein coding or not (for use in presentations).
*Example output:*
```r
[1] "ctrl_stable_dis_df 
total positive DEG count: 6201 
positive coding DEG count: 4286 
positive non-coding DEG count: 1915 percent 
positive genes that are protein coding: 69.1178842122238 
```


---


## Session Info
```r
R version 4.3.1 (2023-06-16)
Platform: aarch64-apple-darwin20 (64-bit)
Running under: macOS 15.3.2

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: America/New_York
tzcode source: internal

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] assertthat_0.2.1            ggpmisc_0.5.5               ggpp_0.5.6                  EnhancedVolcano_1.20.0      ggrepel_0.9.5               svglite_2.1.3               plotly_4.10.4              
 [8] pheatmap_1.0.12             corrplot_0.92               dendextend_1.17.1           vsn_3.70.0                  RUVSeq_1.36.0               edgeR_4.0.16                limma_3.58.1               
[15] EDASeq_2.36.0               ShortRead_1.60.0            GenomicAlignments_1.38.2    Rsamtools_2.18.0            Biostrings_2.70.2           XVector_0.42.0              BiocParallel_1.36.0        
[22] DESeq2_1.42.0               SummarizedExperiment_1.32.0 MatrixGenerics_1.14.0       matrixStats_1.2.0           GenomicRanges_1.54.1        GenomeInfoDb_1.38.6         org.Hs.eg.db_3.18.0        
[29] AnnotationDbi_1.64.1        IRanges_2.36.0              S4Vectors_0.40.2            Biobase_2.62.0              BiocGenerics_0.48.1         ComplexUpset_1.3.3          ggplotify_0.1.2            
[36] cowplot_1.1.3               gdata_3.0.0                 data.table_1.15.0           WebGestaltR_1.0.0           lubridate_1.9.3             forcats_1.0.0               stringr_1.5.1              
[43] dplyr_1.1.4                 purrr_1.0.2                 readr_2.1.5                 tidyr_1.3.1                 tibble_3.2.1                ggplot2_3.5.0               tidyverse_2.0.0            

loaded via a namespace (and not attached):
  [1] splines_4.3.1           BiocIO_1.12.0           bitops_1.0-7            filelock_1.0.3          R.oo_1.26.0             preprocessCore_1.64.0   XML_3.99-0.16.1         lifecycle_1.0.4        
  [9] doParallel_1.0.17       lattice_0.22-5          MASS_7.3-60.0.1         magrittr_2.0.3          rmarkdown_2.25          yaml_2.3.8              doRNG_1.8.6             DBI_1.2.2              
 [17] RColorBrewer_1.1-3      abind_1.4-5             poolr_1.1-1             zlibbioc_1.48.0         R.utils_2.12.3          RCurl_1.98-1.14         yulab.utils_0.1.4       rappdirs_0.3.3         
 [25] GenomeInfoDbData_1.2.11 MatrixModels_0.5-3      codetools_0.2-19        DelayedArray_0.28.0     xml2_1.3.6              tidyselect_1.2.0        farver_2.1.1            viridis_0.6.5          
 [33] BiocFileCache_2.10.1    mathjaxr_1.6-0          jsonlite_1.8.8          survival_3.5-8          iterators_1.0.14        systemfonts_1.0.5       foreach_1.5.2           tools_4.3.1            
 [41] progress_1.2.3          Rcpp_1.0.12             glue_1.7.0              gridExtra_2.3           SparseArray_1.2.4       xfun_0.42               withr_3.0.0             BiocManager_1.30.22    
 [49] fastmap_1.1.1           latticeExtra_0.6-30     fansi_1.0.6             SparseM_1.81            digest_0.6.35           timechange_0.3.0        R6_2.5.1                gridGraphics_0.5-1     
 [57] colorspace_2.1-0        gtools_3.9.5            jpeg_0.1-10             biomaRt_2.58.2          RSQLite_2.3.5           R.methodsS3_1.8.2       utf8_1.2.4              generics_0.1.3         
 [65] rtracklayer_1.62.0      prettyunits_1.2.0       httr_1.4.7              htmlwidgets_1.6.4       S4Arrays_1.2.0          whisker_0.4.1           pkgconfig_2.0.3         gtable_0.3.4           
 [73] blob_1.2.4              hwriter_1.3.2.1         htmltools_0.5.7         scales_1.3.0            png_0.1-8               knitr_1.45              rstudioapi_0.15.0       tzdb_0.4.0             
 [81] rjson_0.2.21            curl_5.2.0              cachem_1.0.8            parallel_4.3.1          restfulr_0.0.15         pillar_1.9.0            grid_4.3.1              vctrs_0.6.5            
 [89] dbplyr_2.4.0            cluster_2.1.6           evaluate_0.23           GenomicFeatures_1.54.4  cli_3.6.2               locfit_1.5-9.8          compiler_4.3.1          rlang_1.1.3            
 [97] crayon_1.5.2            rngtools_1.5.2          interp_1.1-6            aroma.light_3.32.0      affy_1.80.0             fs_1.6.3                stringi_1.8.3           viridisLite_0.4.2      
[105] deldir_2.0-2            munsell_0.5.0           lazyeval_0.2.2          quantreg_5.97           Matrix_1.6-5            hms_1.1.3               patchwork_1.3.0.9000    bit64_4.0.5            
[113] KEGGREST_1.42.0         statmod_1.5.0           apcluster_1.4.13        memoise_2.0.1           affyio_1.72.0           bit_4.0.5               polynom_1.4-1      
```