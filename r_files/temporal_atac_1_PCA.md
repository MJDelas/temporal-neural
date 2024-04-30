ATAC_1 PCA plots
================

# ATAC analysis

``` r
rm(list=ls())

library(DESeq2)
library(RColorBrewer)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
```

Before running, change your working directory in
`/r_inputs/TemporalSpatialNeuralTube_settings.R`

The ATAC analysis starts from the output of nf-core atacseq. File:
`consensus_peaks.mLb.clN.featureCounts.txt` in dir:
`results/bwa/mergedLibrary/macs/broadPeak/consensus/`

### Load settings

Colors, main directory

``` r
source('./r_inputs/TemporalSpatialNeuralTube_settings.R')
```

### Set dirs

``` r
subworkinput="inputs_glialatac_1_eda_pca/"
outdir="outputs_glialatac_1/"
ifelse(!dir.exists(file.path(workingdir,outdir)), dir.create(file.path(workingdir,outdir)), "Directory exists")
```

    ## [1] TRUE

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

## Differential analysis - ALL SAMPLES

``` r
count_matrix <- count_table %>%
  column_to_rownames("Geneid")

## Make metadata file for DESeq

genecolData_first <- data.frame(Sample_ID = colnames(count_matrix))
genecolData_first <- genecolData_first %>% 
  separate(Sample_ID,into=c("Genotype","Day","Gate","NFIAgate","Rep"), sep="_", remove=FALSE) %>%
  mutate(Condition=paste(Genotype,Day,Gate,NFIAgate, sep="_"),
         DayNFIA=paste(Day,NFIAgate,Genotype,sep = "_"),
         DayGate=paste(Day,Gate,sep="_"),
         Experiment=paste(Genotype,Rep,sep="_"),
         NFIAstatus=paste(NFIAgate,Genotype,sep="_"))
genecolData_first <- as.data.frame(unclass(genecolData_first))


dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = genecolData_first,
                              design = ~ Gate)
```

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

``` r
dds <- DESeq(dds)
```

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## -- replacing outliers and refitting for 1182 genes
    ## -- DESeq argument 'minReplicatesForReplace' = 7 
    ## -- original counts are preserved in counts(dds)

    ## estimating dispersions

    ## fitting model and testing

``` r
vsd <- varianceStabilizingTransformation(dds,blind = FALSE)
```

## Export files

Useful for ploting heatmaps elsewhere

``` r
# Export normalized tables for plotting elsewhere
dds_counts <- counts(dds, normalized = TRUE)
vsd_data <- assay(vsd)

write.table(dds_counts, file = paste0(workingdir,outdir,"consensus_peaks.mLb.clN.normCounts.txt"), quote = FALSE, row.names = TRUE)
write.csv(vsd_data, file = paste0(workingdir,outdir,"/consensus_peaks.mLb.vsd.csv"), quote = FALSE)
```

## Differential analysis WT samples

``` r
count_matrix <- count_table %>%
  column_to_rownames("Geneid") %>% 
  select(starts_with("WT"))

## Make metadata file for DESeq

genecolData_first <- data.frame(Sample_ID = colnames(count_matrix))
genecolData_first <- genecolData_first %>% 
  separate(Sample_ID,into=c("Genotype","Day","Gate","NFIAgate","Rep"), sep="_", remove=FALSE) %>%
  mutate(Condition=paste(Genotype,Day,Gate,NFIAgate, sep="_"),
         DayNFIA=paste(Day,NFIAgate,Genotype,sep = "_"),
         DayGate=paste(Day,Gate,sep="_"),
         Experiment=paste(Genotype,Rep,sep="_"),
         NFIAstatus=paste(NFIAgate,Genotype,sep="_"))
genecolData_first <- as.data.frame(unclass(genecolData_first))


dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = genecolData_first,
                              design = ~ Gate)
```

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

``` r
dds <- DESeq(dds)
```

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## -- replacing outliers and refitting for 465 genes
    ## -- DESeq argument 'minReplicatesForReplace' = 7 
    ## -- original counts are preserved in counts(dds)

    ## estimating dispersions

    ## fitting model and testing

``` r
vsd <- varianceStabilizingTransformation(dds,blind = FALSE)
```

## Export files

Useful for ploting heatmaps elsewhere - not for now.

``` r
# Export normalized tables for plotting elsewhere
dds_counts <- counts(dds, normalized = TRUE)
vsd_data <- assay(vsd)
# 
# write.table(dds_counts, file = paste0(workingdir,outdir,"consensus_peaks.mLb.clN.normCounts.txt"), quote = FALSE, row.names = TRUE)
# write.csv(vsd_data, file = paste0(workingdir,outdir,"/consensus_peaks.mLb.vsd.csv"), quote = FALSE)
```

## Plot PCAs

``` r
## to get more than 2 components

## to go it in the top ntop variable genes 
# calculate the variance for each gene
rv <- rowVars(vsd_data)
# select the ntop genes by variance
ntop=30000
select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]

#pca <- prcomp(t(assay(object)[select,]))
t_vsd <- t(vsd_data[select,])
vsd_pca <- prcomp(t_vsd, retx=TRUE, scale. = FALSE)

var_explained <-vsd_pca$sdev^2/sum(vsd_pca$sdev^2)
```

``` r
vsd_pca_plot <- vsd_pca$x %>% 
  as.data.frame %>%
  rownames_to_column("Sample") %>%
  separate(Sample,into=c("Genotype","Day","Gate","NFIAgate","Rep"), sep="_", remove=FALSE) %>%
  mutate(Condition=paste(Genotype,Day,Gate,NFIAgate, sep="_"),
         DayNFIA=factor(paste(Day,NFIAgate,Genotype,sep = "_"), levels=sorted.dayNfia),
         NFIAstatus=paste(NFIAgate,Genotype,sep="_"),
         DayGate=factor(paste(Day,Gate,sep="_"), levels = sorted.DayGate),
         Experiment=paste(Genotype,Rep,sep="_"),
         Day=factor(Day, levels = sorted.day))
  


ggplot(vsd_pca_plot, aes(x=-PC1,y=-PC2,fill=Day,shape=Gate)) +
#  scale_fill_manual(values = colorIZ) +
  geom_point(size=4, alpha=0.9) +
  guides(fill = guide_legend(override.aes=list(shape=21))) +
  scale_shape_manual(values = shapes4_fill_manual) +
  scale_fill_manual(values = color_days) +
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) +
  theme_bw(base_size=16)
```

![](temporal_atac_1_PCA_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

``` r
ggplot(vsd_pca_plot, aes(x=-PC1,y=-PC2,fill=Gate,shape=Day)) +
  scale_fill_manual(values = colorIZ[4:7]) +
  geom_point(size=4, alpha=0.9) +
  guides(fill = guide_legend(override.aes=list(shape=21))) +
  scale_shape_manual(values = shapes4_fill_manual) +
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) +
  theme_bw(base_size=16)
```

![](temporal_atac_1_PCA_files/figure-gfm/unnamed-chunk-10-2.png)<!-- -->

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
    ##  [1] circlize_0.4.16             ComplexHeatmap_2.19.0      
    ##  [3] lubridate_1.9.3             forcats_1.0.0              
    ##  [5] stringr_1.5.1               dplyr_1.1.4                
    ##  [7] purrr_1.0.2                 readr_2.1.5                
    ##  [9] tidyr_1.3.1                 tibble_3.2.1               
    ## [11] ggplot2_3.5.1               tidyverse_2.0.0            
    ## [13] RColorBrewer_1.1-3          DESeq2_1.43.5              
    ## [15] SummarizedExperiment_1.33.3 Biobase_2.63.1             
    ## [17] MatrixGenerics_1.15.1       matrixStats_1.3.0          
    ## [19] GenomicRanges_1.55.4        GenomeInfoDb_1.39.14       
    ## [21] IRanges_2.37.1              S4Vectors_0.41.7           
    ## [23] BiocGenerics_0.49.1        
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] tidyselect_1.2.1        farver_2.1.1            fastmap_1.1.1          
    ##  [4] digest_0.6.35           timechange_0.3.0        lifecycle_1.0.4        
    ##  [7] cluster_2.1.6           magrittr_2.0.3          compiler_4.4.0         
    ## [10] rlang_1.1.3             tools_4.4.0             utf8_1.2.4             
    ## [13] yaml_2.3.8              knitr_1.46              S4Arrays_1.3.7         
    ## [16] labeling_0.4.3          DelayedArray_0.29.9     abind_1.4-5            
    ## [19] BiocParallel_1.37.1     withr_3.0.0             fansi_1.0.6            
    ## [22] colorspace_2.1-0        scales_1.3.0            iterators_1.0.14       
    ## [25] cli_3.6.2               rmarkdown_2.26          crayon_1.5.2           
    ## [28] generics_0.1.3          rstudioapi_0.16.0       httr_1.4.7             
    ## [31] tzdb_0.4.0              rjson_0.2.21            zlibbioc_1.49.3        
    ## [34] parallel_4.4.0          XVector_0.43.1          vctrs_0.6.5            
    ## [37] Matrix_1.7-0            jsonlite_1.8.8          hms_1.1.3              
    ## [40] GetoptLong_1.0.5        clue_0.3-65             locfit_1.5-9.9         
    ## [43] foreach_1.5.2           glue_1.7.0              codetools_0.2-20       
    ## [46] stringi_1.8.3           gtable_0.3.5            shape_1.4.6.1          
    ## [49] UCSC.utils_0.99.7       munsell_0.5.1           pillar_1.9.0           
    ## [52] htmltools_0.5.8.1       GenomeInfoDbData_1.2.12 R6_2.5.1               
    ## [55] doParallel_1.0.17       evaluate_0.23           lattice_0.22-6         
    ## [58] highr_0.10              png_0.1-8               Rcpp_1.0.12            
    ## [61] SparseArray_1.3.7       xfun_0.43               pkgconfig_2.0.3        
    ## [64] GlobalOptions_0.1.2
