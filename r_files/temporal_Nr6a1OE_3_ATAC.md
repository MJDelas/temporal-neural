Nr6a1 down elements in ATAC
================

# Exploring elements that were mis-regulated by Nr6a1 overexpression

``` r
rm(list=ls())
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
subworkinput="inputs_glialnr6a1_3/"
outdir="outputs_glialnr6a1_3/"
ifelse(!dir.exists(file.path(workingdir,outdir)), dir.create(file.path(workingdir,outdir)), "Directory exists")
```

    ## [1] "Directory exists"

## Load data

BED tables from previous script intersected with the CaTS ATAC consensus
to find the consensus elements and plot their accessibility.

Counts/vsd from full GLIAL ATAC

``` r
#counts table
count_atac <- read.table(file = paste0(workingdir,"outputs_glialatac_1/","consensus_peaks.mLb.clN.normCounts.txt"),header=TRUE, stringsAsFactors = FALSE)
vsd_atac <- read.csv(file=paste0(workingdir,"outputs_glialatac_1/","consensus_peaks.mLb.vsd.csv"),header=TRUE, stringsAsFactors = FALSE)


## Annotation table
ann_table <- read.table(file=paste0(workingdir,"inputs_glialatac_1_eda_pca/","consensus_peaks.mLb.clN.annotatePeaks.txt"), header=TRUE, stringsAsFactors = FALSE, sep = "\t")
colnames(ann_table)[1] <- "Peakid"


ann_table_clean <- ann_table %>% 
  select(c("Peakid","Chr","Start","End","Strand","Annotation","Distance.to.TSS","Nearest.PromoterID")) %>%
  separate(Annotation, into = "Annotation_brief", sep = " ", remove = FALSE)
```

    ## Warning: Expected 1 pieces. Additional pieces discarded in 87054 rows [2, 3, 4, 8, 11,
    ## 12, 15, 16, 17, 20, 24, 28, 29, 32, 33, 35, 36, 37, 40, 42, ...].

``` r
nr6a1down_elements <- read.table(file = paste0(workingdir,subworkinput,"Nr6a1_down__within__consensus.bed"),header=FALSE, stringsAsFactors = FALSE)

nr6a1down_elements_ann <- nr6a1down_elements %>% 
  left_join(ann_table_clean, by=c("V4"="Peakid"))
```

## Make ave per condition

For the subset of elements

``` r
## make averages by condition

vsd_sub <- vsd_atac %>%
  filter(X %in% nr6a1down_elements$V4) %>%
  pivot_longer(values_to = "vsd", names_to = "sampleid", MUT_D7_pM_NFIAn_R1:WT_D9_p2_NFIAn_R1) %>%
  separate(sampleid, into = c("geno","day","gate","nfia","rep")) %>%
  mutate(condition=paste(geno,day,gate, sep="_")) %>%
  group_by(X,condition) %>%
  summarise(vsd_ave=mean(vsd)) %>%
  pivot_wider(values_from = vsd_ave, names_from = condition)
```

    ## `summarise()` has grouped output by 'X'. You can override using the `.groups`
    ## argument.

## Plot original ATAC for this subset of elements

``` r
# make a heatmap of CaTS-RNA: this heatmap filters rows of genes with sd zero
# filter genes
vsd_hm <- vsd_sub %>%
  column_to_rownames("X") %>%
  select(starts_with("WT"))

dim(vsd_hm)
```

    ## [1] 40 12

``` r
vsd_hm_z <- t(scale(t(vsd_hm))) 

# order the columns
vsd_hm_z <- vsd_hm_z[,sorted_WTDayGate]

# metadata for the heatmap
genecolData_first <- data.frame(Sample_ID = colnames(vsd_hm_z))
genecolData_first <- genecolData_first %>% 
  separate(Sample_ID,into=c("Genotype","Day","Gate"), sep="_", remove=FALSE) %>%
  mutate(DayGate=paste(Day,Gate,sep="_"))
genecolData_first <- as.data.frame(unclass(genecolData_first))

phen_data <- genecolData_first %>%
  select(c("Sample_ID","DayGate","Day")) %>%
  remove_rownames() %>%
  column_to_rownames("Sample_ID")

ann_color_JD <- list(
  DayGate = c(D5_p1="#abdff4",D5_p2="#f1df9a", D5_pM="#f19aac",
              D7_p1="#55bee8",D7_p2="#e6c444",D7_pM="#e64466",
              D9_p1="#1a91c1",D9_p2="#c19e1a",D9_pM="#c11a3d",
              D11_p1="#0e506b",D11_p2="#6b570e",D11_pM="#7c1127"),
  Day = c(D5="#fadede",D7="#f3aaaa",D9="#e96666",D11="#cf1e1e"))

# Build the annotation for the complex heatmap
colAnn <- HeatmapAnnotation(
    df = phen_data,
    which = 'col', # 'col' (samples) or 'row' (gene) annotation?
    na_col = 'white', # default colour for any NA values in the annotation data-frame, 'ann'
    col = ann_color_JD,
    annotation_height = 0.6,
    annotation_width = unit(1, 'cm'),
    gap = unit(1, 'mm'))

# Annotated heatmap with selected colors
hm_colors = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100)


hmap_catsATAC <- Heatmap(vsd_hm_z,
    name = 'Z-score',

    col = hm_colors,

    # row (gene) parameters
      cluster_rows = TRUE,
      show_row_dend = TRUE,
      #row_title = 'Statistically significant genes',
      row_title_side = 'left',
      row_title_gp = gpar(fontsize = 12,  fontface = 'bold'),
      row_title_rot = 90,
      show_row_names = FALSE,
      row_names_gp = gpar(fontsize = 12, fontface = 'bold'),
      row_names_side = 'left',
      row_dend_width = unit(25,'mm'),

    # column (sample) parameters
      cluster_columns = FALSE,
      show_column_dend = TRUE,
      column_title = '',
      column_title_side = 'bottom',
      column_title_gp = gpar(fontsize = 12, fontface = 'bold'),
      column_title_rot = 0,
      show_column_names = FALSE,
      column_names_gp = gpar(fontsize = 8),
      column_names_max_height = unit(10, 'cm'),
      column_dend_height = unit(25,'mm'),

    # cluster methods for rows and columns
      clustering_distance_columns = function(x) as.dist(1 - cor(t(x))),
      clustering_method_columns = 'ward.D2',
      clustering_distance_rows = function(x) as.dist(1 - cor(t(x))),
      clustering_method_rows = 'ward.D2',

    # specify top and bottom annotations
      top_annotation = colAnn)


 
draw(hmap_catsATAC,
    heatmap_legend_side = 'left',
    annotation_legend_side = 'left',
    row_sub_title_side = 'left')
```

![](temporal_Nr6a1OE_3_ATAC_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
# print heatmap
pdf(paste0(workingdir,outdir,"Heatmap_Nr6a1down_elements_in_CaTSATAC.pdf"), width = 10, height = 2) 

draw(hmap_catsATAC,
    heatmap_legend_side = 'left',
    annotation_legend_side = 'left',
    row_sub_title_side = 'left')

dev.off()
```

    ## quartz_off_screen 
    ##                 2

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
    ## [1] grid      stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] ComplexHeatmap_2.19.0 lubridate_1.9.3       forcats_1.0.0        
    ##  [4] stringr_1.5.1         dplyr_1.1.4           purrr_1.0.2          
    ##  [7] readr_2.1.5           tidyr_1.3.1           tibble_3.2.1         
    ## [10] ggplot2_3.5.1         tidyverse_2.0.0       RColorBrewer_1.1-3   
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] utf8_1.2.4          generics_0.1.3      shape_1.4.6.1      
    ##  [4] stringi_1.8.3       hms_1.1.3           digest_0.6.35      
    ##  [7] magrittr_2.0.3      evaluate_0.23       timechange_0.3.0   
    ## [10] iterators_1.0.14    circlize_0.4.16     fastmap_1.1.1      
    ## [13] foreach_1.5.2       doParallel_1.0.17   GlobalOptions_0.1.2
    ## [16] fansi_1.0.6         scales_1.3.0        codetools_0.2-20   
    ## [19] cli_3.6.2           rlang_1.1.3         crayon_1.5.2       
    ## [22] munsell_0.5.1       withr_3.0.0         yaml_2.3.8         
    ## [25] tools_4.4.0         parallel_4.4.0      tzdb_0.4.0         
    ## [28] colorspace_2.1-0    BiocGenerics_0.49.1 GetoptLong_1.0.5   
    ## [31] vctrs_0.6.5         R6_2.5.1            png_0.1-8          
    ## [34] stats4_4.4.0        matrixStats_1.3.0   lifecycle_1.0.4    
    ## [37] S4Vectors_0.41.7    IRanges_2.37.1      clue_0.3-65        
    ## [40] cluster_2.1.6       pkgconfig_2.0.3     pillar_1.9.0       
    ## [43] gtable_0.3.5        glue_1.7.0          highr_0.10         
    ## [46] xfun_0.43           tidyselect_1.2.1    rstudioapi_0.16.0  
    ## [49] knitr_1.46          rjson_0.2.21        htmltools_0.5.8.1  
    ## [52] rmarkdown_2.26      Cairo_1.6-2         compiler_4.4.0
