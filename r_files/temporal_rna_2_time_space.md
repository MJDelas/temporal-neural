RNA_2_time_vs_domain
================

# RNA analysis

Differential accessibility in groups of samples

Generate the comparison and export the data.

Start the analysis in the next script by importing these data so it is
more efficient.

``` r
rm(list=ls())

library(DESeq2)
library(RColorBrewer)
library(ComplexHeatmap)
library(tximport)
library(tidyverse)
```

### Load settings

Colors, main directory

``` r
source('./r_inputs/TemporalSpatialNeuralTube_settings.R')
```

### Set dirs

``` r
outdir="outputs_glialrna_2_time_space/"
subworkinput="inputs_glialRNA_1_qc/"
ifelse(!dir.exists(file.path(workingdir,outdir)), dir.create(file.path(workingdir,outdir)), "Directory exists")
```

    ## [1] TRUE

``` r
suboutdir1="output_Domain_Specific/"
ifelse(!dir.exists(file.path(workingdir,outdir,suboutdir1)), dir.create(file.path(workingdir,outdir,suboutdir1)), "Directory exists")
```

    ## [1] TRUE

``` r
suboutdir2="output_Time_Specific/"
ifelse(!dir.exists(file.path(workingdir,outdir,suboutdir2)), dir.create(file.path(workingdir,outdir,suboutdir2)), "Directory exists")
```

    ## [1] TRUE

## Load data

For RNA analysis, we are using the output of star_salmon from
`nf-core rnaseq`. Pipeline outpout directory: `results/star_salmon`

``` r
#salmon counts from pipeline, import
path_files =  list.files(paste0(workingdir,subworkinput)) 
samples = data.frame(run=path_files, stringsAsFactors = FALSE) %>%
  filter(str_detect(run, "MUT|WT"))

files <- file.path(paste0(workingdir,subworkinput), samples$run, "quant.sf")
names(files) <- samples$run
all(file.exists(files))
```

    ## [1] TRUE

``` r
#from pipeline
tx2gene = read_tsv(paste0(paste0(workingdir,subworkinput),"/salmon_tx2gene.tsv"))
```

    ## New names:
    ## Rows: 35210 Columns: 3
    ## ── Column specification
    ## ──────────────────────────────────────────────────────── Delimiter: "\t" chr
    ## (3): NM_001011874, Xkr4...2, Xkr4...3
    ## ℹ Use `spec()` to retrieve the full column specification for this data. ℹ
    ## Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## • `Xkr4` -> `Xkr4...2`
    ## • `Xkr4` -> `Xkr4...3`

``` r
txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene)
```

    ## reading in files with read_tsv
    ## 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 
    ## transcripts missing from tx2gene: 1
    ## summarizing abundance
    ## summarizing counts
    ## summarizing length

``` r
#head(txi.salmon$counts)
```

## Differential analysis between domains for each timepoint

Targeted diff analysis in subsets of samples:

Wild type only: D5: pairwise for p1, p2, pMN D7: pairwise for p1, p2,
pMN D9: pairwise for p1, p2, pMN D11: pairwise for p1, p2, pMN

This DESeq2 analysis is done by subsetting samples. I subset the
`txi.salmon$counts` table and then use `DESeqDataSetFromMatrix` with
`round(sub_counts)` to perform the differential analysis.

``` r
count_table <- txi.salmon$counts

count_matrix <- count_table %>%
  as.data.frame() %>%
  rownames_to_column("GeneID") %>%
  filter(!str_starts(GeneID, "ERCC")) %>%
  column_to_rownames("GeneID")


#subset 
timepoint=c("_D5_","_D7_","_D9_","_D11_")
#matrix
allgates=matrix(c("_pM_","_p2_","_pM_","_p1_","_p1_","_p2_"),
                 nrow=2,
                 ncol=3)
geno="WT_"

comparisons <- allgates

# count_matrix <- count_table %>%
#   column_to_rownames("Geneid")


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
         NFIAstatus=paste(NFIAgate,Genotype,sep="_"),
         Gate=factor(Gate, levels = gates))
      genecolData_sub <- as.data.frame(unclass(genecolData_sub))
      
      
      dds_sub <- DESeqDataSetFromMatrix(countData =  round(sub_counts),
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
      file = paste0(workingdir,outdir,suboutdir1,"CountsNormalized_",timepoints,resultsNames(dds_sub)[2],".txt"),
          quote = FALSE, row.names = TRUE)
      write.csv(vsd_sub_data,
          paste0(workingdir,outdir,suboutdir1,"VSData_",timepoints,resultsNames(dds_sub)[2],".csv"),
          quote = FALSE)
      write.table(results_sub,
          file = paste0(workingdir,outdir,suboutdir1,"Results_DESeq_",timepoints,resultsNames(dds_sub)[2],".txt"),
          quote = FALSE, row.names = TRUE)

      results_return <- results_sub %>% as.data.frame() %>% rownames_to_column("Geneid")
      results_return$Comparison <- paste0("Comp_",timepoints,resultsNames(dds_sub)[2])
      results_return

  })
}) 
```

    ## converting counts to integer mode

    ## factor levels were dropped which had no samples

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## converting counts to integer mode

    ## factor levels were dropped which had no samples

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## converting counts to integer mode

    ## factor levels were dropped which had no samples

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## converting counts to integer mode

    ## factor levels were dropped which had no samples

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## converting counts to integer mode

    ## factor levels were dropped which had no samples

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## converting counts to integer mode

    ## factor levels were dropped which had no samples

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## converting counts to integer mode

    ## factor levels were dropped which had no samples

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## converting counts to integer mode

    ## factor levels were dropped which had no samples

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## converting counts to integer mode

    ## factor levels were dropped which had no samples

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## converting counts to integer mode

    ## factor levels were dropped which had no samples

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## converting counts to integer mode

    ## factor levels were dropped which had no samples

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## converting counts to integer mode

    ## factor levels were dropped which had no samples

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
      
      dds_sub <- DESeqDataSetFromMatrix(countData = round(sub_counts),
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
      file = paste0(workingdir,outdir,suboutdir2,"CountsNormalized_",celltypes,resultsNames(dds_sub)[2],".txt"),
          quote = FALSE, row.names = TRUE)
      write.csv(vsd_sub_data,
          paste0(workingdir,outdir,suboutdir2,"VSData_",celltypes,resultsNames(dds_sub)[2],".csv"),
          quote = FALSE)
      write.table(results_sub,
          file = paste0(workingdir,outdir,suboutdir2,"Results_DESeq_",celltypes,resultsNames(dds_sub)[2],".txt"),
          quote = FALSE, row.names = TRUE)

      results_return <- results_sub %>% as.data.frame() %>% rownames_to_column("Geneid")
      results_return$Comparison <- paste0("Comp_",celltypes,resultsNames(dds_sub)[2])
      results_return

  })
}) 
```

    ## converting counts to integer mode

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## converting counts to integer mode

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## converting counts to integer mode

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## converting counts to integer mode

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## converting counts to integer mode

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## converting counts to integer mode

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## converting counts to integer mode

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## converting counts to integer mode

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## converting counts to integer mode

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## converting counts to integer mode

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## converting counts to integer mode

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## converting counts to integer mode

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## converting counts to integer mode

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## converting counts to integer mode

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## converting counts to integer mode

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## converting counts to integer mode

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## converting counts to integer mode

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## converting counts to integer mode

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
    ##  [1] lubridate_1.9.3             forcats_1.0.0              
    ##  [3] stringr_1.5.1               dplyr_1.1.4                
    ##  [5] purrr_1.0.2                 readr_2.1.5                
    ##  [7] tidyr_1.3.1                 tibble_3.2.1               
    ##  [9] ggplot2_3.5.1               tidyverse_2.0.0            
    ## [11] tximport_1.31.1             ComplexHeatmap_2.19.0      
    ## [13] RColorBrewer_1.1-3          DESeq2_1.43.5              
    ## [15] SummarizedExperiment_1.33.3 Biobase_2.63.1             
    ## [17] MatrixGenerics_1.15.1       matrixStats_1.3.0          
    ## [19] GenomicRanges_1.55.4        GenomeInfoDb_1.39.14       
    ## [21] IRanges_2.37.1              S4Vectors_0.41.7           
    ## [23] BiocGenerics_0.49.1        
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] tidyselect_1.2.1        fastmap_1.1.1           digest_0.6.35          
    ##  [4] timechange_0.3.0        lifecycle_1.0.4         cluster_2.1.6          
    ##  [7] magrittr_2.0.3          compiler_4.4.0          rlang_1.1.3            
    ## [10] tools_4.4.0             utf8_1.2.4              yaml_2.3.8             
    ## [13] knitr_1.46              S4Arrays_1.3.7          bit_4.0.5              
    ## [16] DelayedArray_0.29.9     abind_1.4-5             BiocParallel_1.37.1    
    ## [19] withr_3.0.0             fansi_1.0.6             colorspace_2.1-0       
    ## [22] scales_1.3.0            iterators_1.0.14        cli_3.6.2              
    ## [25] rmarkdown_2.26          crayon_1.5.2            generics_0.1.3         
    ## [28] rstudioapi_0.16.0       httr_1.4.7              tzdb_0.4.0             
    ## [31] rjson_0.2.21            zlibbioc_1.49.3         parallel_4.4.0         
    ## [34] XVector_0.43.1          vctrs_0.6.5             Matrix_1.7-0           
    ## [37] jsonlite_1.8.8          GetoptLong_1.0.5        hms_1.1.3              
    ## [40] bit64_4.0.5             clue_0.3-65             locfit_1.5-9.9         
    ## [43] foreach_1.5.2           glue_1.7.0              codetools_0.2-20       
    ## [46] stringi_1.8.3           gtable_0.3.5            shape_1.4.6.1          
    ## [49] UCSC.utils_0.99.7       munsell_0.5.1           pillar_1.9.0           
    ## [52] htmltools_0.5.8.1       GenomeInfoDbData_1.2.12 circlize_0.4.16        
    ## [55] R6_2.5.1                doParallel_1.0.17       vroom_1.6.5            
    ## [58] evaluate_0.23           lattice_0.22-6          png_0.1-8              
    ## [61] Rcpp_1.0.12             SparseArray_1.3.7       xfun_0.43              
    ## [64] pkgconfig_2.0.3         GlobalOptions_0.1.2
