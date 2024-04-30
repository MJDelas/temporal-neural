RNA Deseq
================

# RNA analysis

Re-analysis after nf-core/rnaseq 3.5

``` r
rm(list=ls())

library(DESeq2)
library(RColorBrewer)
library(tidyverse)
library(ComplexHeatmap)
library(tximport)
```

### Load settings

Colors, main directory

``` r
source('./r_inputs/TemporalSpatialNeuralTube_settings.R')
```

### Set dirs

``` r
subworkinput="inputs_glialRNA_1_qc/"

outdir="outputs_glialRNA_1/"
ifelse(!dir.exists(file.path(workingdir,outdir)), dir.create(file.path(workingdir,outdir)), "Directory exists")
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

genes named Geneid

## qc samples

``` r
count.table <- txi.salmon$counts
```

# Diff expression

The `DESeqDataSetFromTximport` was just re-loading the full matrix (with
ERCC). I am trying rounding the data after subsetting and importing as
usual with `DESeqDataSetFromMatrix`.

``` r
count.matrix <- count.table %>%
  as.data.frame() %>%
  rownames_to_column("GeneID") %>%
  filter(!str_starts(GeneID, "ERCC")) %>%
  column_to_rownames("GeneID")

## Make metadata file for DESeq

genecolData_first <- data.frame(Sample_ID = colnames(count.matrix))
genecolData_first <- genecolData_first %>% 
  separate(Sample_ID,into=c("Genotype","Day","Gate","NFIAgate","Rep"), sep="_", remove=FALSE) %>%
  mutate(Condition=paste(Genotype,Day,Gate,NFIAgate, sep="_"),
         DayNFIA=paste(Day,NFIAgate,Genotype,sep = "_"),
         DayGate=paste(Day,Gate,sep="_"),
         Experiment=paste(Genotype,Rep,sep="_"),
         NFIAstatus=paste(NFIAgate,Genotype,sep="_"))
genecolData_first <- as.data.frame(unclass(genecolData_first))


# dds <- DESeqDataSetFromTximport(txi = txi.salmon,
#                               colData = genecolData_first,
#                               design = ~ Gate)

dds <- DESeqDataSetFromMatrix(countData = round(count.matrix),
                              colData = genecolData_first,
                              design = ~ Gate)
```

    ## converting counts to integer mode

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

    ## -- replacing outliers and refitting for 606 genes
    ## -- DESeq argument 'minReplicatesForReplace' = 7 
    ## -- original counts are preserved in counts(dds)

    ## estimating dispersions

    ## fitting model and testing

``` r
vsd <- varianceStabilizingTransformation(dds,blind = FALSE)
```

## Export files

``` r
# Export normalized tables for plotting elsewhere
dds_counts <- counts(dds, normalized = TRUE)
vsd_data <- assay(vsd)

write.table(dds_counts, file = paste0(workingdir,outdir,"featurecounts.normCounts.txt"), quote = FALSE, row.names = TRUE)
write.csv(vsd_data, file = paste0(workingdir,outdir,"featurecounts.vsd.csv"), quote = FALSE)
```

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
    ##  [1] tximport_1.31.1             ComplexHeatmap_2.19.0      
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
    ## [37] jsonlite_1.8.8          hms_1.1.3               GetoptLong_1.0.5       
    ## [40] bit64_4.0.5             clue_0.3-65             locfit_1.5-9.9         
    ## [43] foreach_1.5.2           glue_1.7.0              codetools_0.2-20       
    ## [46] stringi_1.8.3           gtable_0.3.5            shape_1.4.6.1          
    ## [49] UCSC.utils_0.99.7       munsell_0.5.1           pillar_1.9.0           
    ## [52] htmltools_0.5.8.1       GenomeInfoDbData_1.2.12 circlize_0.4.16        
    ## [55] R6_2.5.1                doParallel_1.0.17       vroom_1.6.5            
    ## [58] evaluate_0.23           lattice_0.22-6          png_0.1-8              
    ## [61] Rcpp_1.0.12             SparseArray_1.3.7       xfun_0.43              
    ## [64] pkgconfig_2.0.3         GlobalOptions_0.1.2
