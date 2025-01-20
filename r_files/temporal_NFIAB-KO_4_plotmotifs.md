ATAC_NFIAB_KO_motifs
================

# ATAC analysis

Differential accessibility defects in NIFA/B KO.

``` r
rm(list=ls())


library(RColorBrewer)
library(tidyverse)
```

### Load settings

Colors, main directory

``` r
source('./r_inputs/TemporalSpatialNeuralTube_settings.R')
```

### Set dirs

``` r
outdir="outputs_NFIAB_ATAC/"
subworkinput="outputs_NFIAB_motifs/"
ifelse(!dir.exists(file.path(workingdir,outdir)), dir.create(file.path(workingdir,outdir)), "Directory exists")
```

    ## [1] "Directory exists"

## Load data

``` r
ame_results_arch_topmotif <- read.table(file=paste0(workingdir,subworkinput,"Motifs_top_archtype_NFIA_depen_and_independent.txt"),header=TRUE, stringsAsFactors = FALSE)
```

``` r
top_motifs <- ame_results_arch_topmotif %>% 
  select(rank,motif_id,pvalue,adj.pvalue,cluster,Cluster_ID, Database, Name, Seed_motif, Archetype) %>% 
  group_by(cluster) %>%
  slice_min(order_by=adj.pvalue, n=5)

top_motifs
```

    ## # A tibble: 11 × 10
    ## # Groups:   cluster [2]
    ##     rank motif_id            pvalue adj.pvalue cluster Cluster_ID Database Name 
    ##    <int> <chr>                <dbl>      <dbl> <chr>        <int> <chr>    <chr>
    ##  1     1 SOX4_MOUSE.H11M… 3.53e-144  1.62e-140 NFIA I…         89 HOCOMOC… SOX/1
    ##  2     2 SOX4_HUMAN.H11M… 3.53e-144  1.62e-140 NFIA I…         89 HOCOMOC… SOX/1
    ##  3     7 FOXP1_HUMAN.H11… 4.94e-105  2.3 e-101 NFIA I…         79 HOCOMOC… FOX/4
    ##  4     8 RFX3_MOUSE.H11M… 2.61e- 92  1.85e- 89 NFIA I…        247 HOCOMOC… RFX/1
    ##  5     9 EGR2_HUMAN.H11M… 3.21e- 89  9.19e- 86 NFIA I…        109 HOCOMOC… KLF/…
    ##  6     1 NFIX_NFI_1       1.07e-144  9.65e-142 NFIA d…        189 Taipale… NFI/1
    ##  7     2 NFIC_HUMAN.H11M… 1.82e-125  2.94e-122 NFIA d…        128 HOCOMOC… NFI/3
    ##  8    10 SOX2_HUMAN.H11M… 3.52e- 40  5.05e- 37 NFIA d…         89 HOCOMOC… SOX/1
    ##  9    14 FOXP1_HUMAN.H11… 5.25e- 36  7.94e- 33 NFIA d…         79 HOCOMOC… FOX/4
    ## 10    15 ETS2_MOUSE.H11M… 5.42e- 35  1.04e- 31 NFIA d…         96 HOCOMOC… ETS/2
    ## 11    16 ETS2_HUMAN.H11M… 5.42e- 35  1.04e- 31 NFIA d…         96 HOCOMOC… ETS/2
    ## # ℹ 2 more variables: Seed_motif <chr>, Archetype <chr>

``` r
sessionInfo()
```

    ## R version 4.4.0 (2024-04-24)
    ## Platform: aarch64-apple-darwin20
    ## Running under: macOS 15.2
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
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] lubridate_1.9.3    forcats_1.0.0      stringr_1.5.1      dplyr_1.1.4       
    ##  [5] purrr_1.0.2        readr_2.1.5        tidyr_1.3.1        tibble_3.2.1      
    ##  [9] ggplot2_3.5.1      tidyverse_2.0.0    RColorBrewer_1.1-3
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] gtable_0.3.5      compiler_4.4.0    tidyselect_1.2.1  scales_1.3.0     
    ##  [5] yaml_2.3.8        fastmap_1.2.0     R6_2.5.1          generics_0.1.3   
    ##  [9] knitr_1.47        munsell_0.5.1     pillar_1.9.0      tzdb_0.4.0       
    ## [13] rlang_1.1.4       utf8_1.2.4        stringi_1.8.4     xfun_0.44        
    ## [17] timechange_0.3.0  cli_3.6.2         withr_3.0.0       magrittr_2.0.3   
    ## [21] digest_0.6.35     grid_4.4.0        rstudioapi_0.16.0 hms_1.1.3        
    ## [25] lifecycle_1.0.4   vctrs_0.6.5       evaluate_0.23     glue_1.7.0       
    ## [29] fansi_1.0.6       colorspace_2.1-0  rmarkdown_2.27    tools_4.4.0      
    ## [33] pkgconfig_2.0.3   htmltools_0.5.8.1
