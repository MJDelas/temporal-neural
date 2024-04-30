ATAC_2_time_vs_domain
================

# ATAC analysis

Differential accessibility in groups of samples

Generate the comparison and export the data.

Start the analysis in the next script by importing these data so it is
more efficient for exploring data.

``` r
rm(list=ls())

library(DESeq2)
library(RColorBrewer)
library(tidyverse)
library(ComplexHeatmap)
```

### Load settings

Colors, main directory

``` r
source('./r_inputs/TemporalSpatialNeuralTube_settings.R')
```

### Set dirs

``` r
outdir="outputs_glialatac_2_time_space/"
subworkinput="inputs_glialatac_1_eda_pca/"
ifelse(!dir.exists(file.path(workingdir,outdir)), dir.create(file.path(workingdir,outdir)), "Directory exists")
```

    ## [1] "Directory exists"

``` r
suboutdir1="output_Domain_Specific/"
ifelse(!dir.exists(file.path(workingdir,outdir,suboutdir1)), dir.create(file.path(workingdir,outdir,suboutdir1)), "Directory exists")
```

    ## [1] "Directory exists"

``` r
suboutdir2="output_Time_Specific/"
ifelse(!dir.exists(file.path(workingdir,outdir,suboutdir2)), dir.create(file.path(workingdir,outdir,suboutdir2)), "Directory exists")
```

    ## [1] "Directory exists"

## Load data

``` r
#counts table
count_table <- read.table(file=paste0(workingdir,subworkinput,"consensus_peaks.mLb.clN.featureCounts.txt"),header=TRUE, stringsAsFactors = FALSE)

# clean colnames
colnames(count_table) <- gsub(".mLb.clN.bam","",colnames(count_table))

# we do not need coordinates
count_table <- count_table %>%
  select("Geneid", starts_with(c("MUT","WT")))


## Annotation table
ann_table <- read.table(file=paste0(workingdir,subworkinput,"consensus_peaks.mLb.clN.annotatePeaks.txt"), header=TRUE, stringsAsFactors = FALSE, sep = "\t")
colnames(ann_table)[1] <- "Peakid"


ann_table_clean <- ann_table %>% 
  select(c("Peakid","Chr","Start","End","Strand","Annotation","Distance.to.TSS","Nearest.PromoterID")) %>%
  separate(Annotation, into = "Annotation_brief", sep = " ", remove = FALSE)
```

    ## Warning: Expected 1 pieces. Additional pieces discarded in 87054 rows [2, 3, 4, 8, 11,
    ## 12, 15, 16, 17, 20, 24, 28, 29, 32, 33, 35, 36, 37, 40, 42, ...].

## Differential analysis between domains for each timepoint

Targeted diff analysis in subsets of samples:

Wild type only: D5: pairwise for p1, p2, pMN D7: pairwise for p1, p2,
pMN D9: pairwise for p1, p2, pMN D11: pairwise for p1, p2, pMN

``` r
#subset 
timepoint=c("_D5_","_D7_","_D9_","_D11_")
#matrix
allgates=matrix(c("_pM_","_p2_","_pM_","_p1_","_p1_","_p2_"),
                 nrow=2,
                 ncol=3)
geno="WT_"

comparisons <- allgates

count_matrix <- count_table %>%
  column_to_rownames("Geneid")


PairWiseDEseq <- lapply(c(1:length(timepoint)),function (x) {
  lapply(c(1:ncol(allgates)), function (y) {
      timepoints <- timepoint[x]
      celltypes <- allgates[,y]
      sub_counts <- count_matrix %>%
        dplyr::select(starts_with(geno)  & contains(celltypes) & contains(timepoints))
      
      ## Make metadata file for DESeq
      genecolData_sub <- data.frame(Sample_ID = colnames(sub_counts))
      genecolData_sub <- genecolData_sub %>% 
        separate(Sample_ID,into=c("Genotype","Day","Gate","NFIAgate","Rep"), sep="_", remove=FALSE) %>%
        mutate(Condition=paste(Genotype,Day,Gate,NFIAgate, sep="_"),
         DayNFIA=paste(Day,NFIAgate,Genotype,sep = "_"),
         DayGate=paste(Day,Gate,sep="_"),
         Experiment=paste(Genotype,Rep,sep="_"),
         NFIAstatus=paste(NFIAgate,Genotype,sep="_"))
      genecolData_sub <- as.data.frame(unclass(genecolData_sub))
      
      dds_sub <- DESeqDataSetFromMatrix(countData = sub_counts,
                                    colData = genecolData_sub,
                                    design = ~ Gate)
      
      dds_sub <- DESeq(dds_sub)
      
      vsd_sub <- varianceStabilizingTransformation(dds_sub,blind = FALSE)
      
      # Export normalized tables for plotting elsewhere
      dds_sub_counts <- counts(dds_sub, normalized = TRUE)
      vsd_sub_data <- assay(vsd_sub)
      
      results_sub <- results(dds_sub)

      # plotMA(results_sub,ylim=c(-8,8))

      ## Export files
      
      write.table(dds_sub_counts,
      file = paste0(workingdir,outdir,suboutdir1,"CountsNormalized_",timepoints, resultsNames(dds_sub)[2],".txt"),
          quote = FALSE, row.names = TRUE)
      write.csv(vsd_sub_data,
          paste0(workingdir,outdir,suboutdir1,"VSData_",timepoints, resultsNames(dds_sub)[2],".csv"),
          quote = FALSE)
      write.table(results_sub,
          file = paste0(workingdir,outdir,suboutdir1,"Results_DESeq_",timepoints, resultsNames(dds_sub)[2],".txt"),
          quote = FALSE, row.names = TRUE)

      results_return <- results_sub %>% as.data.frame() %>% rownames_to_column("Geneid")
      results_return$Comparison <- paste0("Comp_",timepoints, resultsNames(dds_sub)[2])
      results_return

  })
}) 
```

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

## Differential analysis between timepoints for each domain

Targeted diff analysis in subsets of samples:

Wild type only: p1: pairwise D5-D7, D7-D9, D9-D11, D5-D9, D5-D11, D7-D11
p2: pairwise D5-D7, D7-D9, D9-D11, D5-D9, D5-D11, D7-D11 pM: pairwise
D5-D7, D7-D9, D9-D11, D5-D9, D5-D11, D7-D11

``` r
#subset 
allgates=c("_p2_","_pM_","_p1_")
#matrix
timepoints=matrix(c("_D5_","_D7_","_D7_","_D9_","_D9_","_D11_","_D5_","_D9_","_D5_","_D11_","_D7_","_D11_"),
                 nrow=2,
                 ncol=6)
geno="WT_"

comparisons <- timepoints

count_matrix <- count_table %>%
  column_to_rownames("Geneid")


PairWiseDEseq <- lapply(c(1:length(allgates)),function (x) {
  lapply(c(1:ncol(timepoints)), function (y) {
      celltypes <- allgates[x]
      days <- timepoints[,y]
      sub_counts <- count_matrix %>%
        dplyr::select(starts_with(geno)  & contains(celltypes) & contains(days))
      
      ## Make metadata file for DESeq
      genecolData_sub <- data.frame(Sample_ID = colnames(sub_counts))
      genecolData_sub <- genecolData_sub %>% 
        separate(Sample_ID,into=c("Genotype","Day","Gate","NFIAgate","Rep"), sep="_", remove=FALSE) %>%
        mutate(Condition=paste(Genotype,Day,Gate,NFIAgate, sep="_"),
         DayNFIA=paste(Day,NFIAgate,Genotype,sep = "_"),
         DayGate=paste(Day,Gate,sep="_"),
         Experiment=paste(Genotype,Rep,sep="_"),
         NFIAstatus=paste(NFIAgate,Genotype,sep="_"))
      genecolData_sub <- as.data.frame(unclass(genecolData_sub))
      
      dds_sub <- DESeqDataSetFromMatrix(countData = sub_counts,
                                    colData = genecolData_sub,
                                    design = ~ Day)
      
      dds_sub <- DESeq(dds_sub)
      
      vsd_sub <- varianceStabilizingTransformation(dds_sub,blind = FALSE)
      
      # Export normalized tables for plotting elsewhere
      dds_sub_counts <- counts(dds_sub, normalized = TRUE)
      vsd_sub_data <- assay(vsd_sub)
      
      results_sub <- results(dds_sub)

      #plotMA(results_sub,ylim=c(-8,8))

      ## Export files
      
      write.table(dds_sub_counts,
      file = paste0(workingdir,outdir,suboutdir2,"CountsNormalized_",celltypes, resultsNames(dds_sub)[2],".txt"),
          quote = FALSE, row.names = TRUE)
      write.csv(vsd_sub_data,
          paste0(workingdir,outdir,suboutdir2,"VSData_",celltypes, resultsNames(dds_sub)[2],".csv"),
          quote = FALSE)
      write.table(results_sub,
          file = paste0(workingdir,outdir,suboutdir2,"Results_DESeq_",celltypes, resultsNames(dds_sub)[2],".txt"),
          quote = FALSE, row.names = TRUE)

      results_return <- results_sub %>% as.data.frame() %>% rownames_to_column("Geneid")
      results_return$Comparison <- paste0("Comp_",celltypes, resultsNames(dds_sub)[2])
      results_return

  })
}) 
```

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

``` r
sessionInfo()
```

    ## R version 4.4.0 (2024-04-24)
    ## Platform: aarch64-apple-darwin20
    ## Running under: macOS Sonoma 14.4.1
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRblas.0.dylib 
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## time zone: Europe/London
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] grid      stats4    stats     graphics  grDevices utils     datasets 
    ## [8] methods   base     
    ## 
    ## other attached packages:
    ##  [1] ComplexHeatmap_2.19.0       lubridate_1.9.3            
    ##  [3] forcats_1.0.0               stringr_1.5.1              
    ##  [5] dplyr_1.1.4                 purrr_1.0.2                
    ##  [7] readr_2.1.5                 tidyr_1.3.1                
    ##  [9] tibble_3.2.1                ggplot2_3.5.1              
    ## [11] tidyverse_2.0.0             RColorBrewer_1.1-3         
    ## [13] DESeq2_1.43.5               SummarizedExperiment_1.33.3
    ## [15] Biobase_2.63.1              MatrixGenerics_1.15.1      
    ## [17] matrixStats_1.3.0           GenomicRanges_1.55.4       
    ## [19] GenomeInfoDb_1.39.14        IRanges_2.37.1             
    ## [21] S4Vectors_0.41.7            BiocGenerics_0.49.1        
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] shape_1.4.6.1           circlize_0.4.16         gtable_0.3.5           
    ##  [4] rjson_0.2.21            xfun_0.43               GlobalOptions_0.1.2    
    ##  [7] lattice_0.22-6          tzdb_0.4.0              vctrs_0.6.5            
    ## [10] tools_4.4.0             generics_0.1.3          parallel_4.4.0         
    ## [13] fansi_1.0.6             cluster_2.1.6           pkgconfig_2.0.3        
    ## [16] Matrix_1.7-0            lifecycle_1.0.4         GenomeInfoDbData_1.2.12
    ## [19] compiler_4.4.0          munsell_0.5.1           codetools_0.2-20       
    ## [22] clue_0.3-65             htmltools_0.5.8.1       yaml_2.3.8             
    ## [25] pillar_1.9.0            crayon_1.5.2            BiocParallel_1.37.1    
    ## [28] DelayedArray_0.29.9     iterators_1.0.14        foreach_1.5.2          
    ## [31] abind_1.4-5             tidyselect_1.2.1        locfit_1.5-9.9         
    ## [34] digest_0.6.35           stringi_1.8.3           fastmap_1.1.1          
    ## [37] colorspace_2.1-0        cli_3.6.2               SparseArray_1.3.7      
    ## [40] magrittr_2.0.3          S4Arrays_1.3.7          utf8_1.2.4             
    ## [43] withr_3.0.0             scales_1.3.0            UCSC.utils_0.99.7      
    ## [46] timechange_0.3.0        rmarkdown_2.26          XVector_0.43.1         
    ## [49] httr_1.4.7              GetoptLong_1.0.5        png_0.1-8              
    ## [52] hms_1.1.3               evaluate_0.23           knitr_1.46             
    ## [55] doParallel_1.0.17       rlang_1.1.3             Rcpp_1.0.12            
    ## [58] glue_1.7.0              rstudioapi_0.16.0       jsonlite_1.8.8         
    ## [61] R6_2.5.1                zlibbioc_1.49.3
