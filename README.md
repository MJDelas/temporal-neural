# The cis-regulatory logic integrating spatial and temporal patterning in the vertebrate neural tube

Code associated with the analysis of the manuscript Zhang et al (2024).

## Data availability

- [GEO superseries]:(https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE264172)

## Data processing

1. ATAC

2. RNA

## Code in this repo

- [ATAC 1](r_files/temporal_atac_1_PCA.md): 
    - Fig 1 ATAC PCA for all cell types across time.
- [ATAC 2](r_files/temporal_atac_2_time_space.md): 
    - Calculate ATAC diff accessibility for all pairwise comparisons
- [ATAC 3](r_files/temporal_atac_2_time_space_plot.md): 
    - Fig 1 ATAC space vs time quantifications
    - Fig 1 temporal heatmap

- [RNA 1](r_files/temporal_rna_1_export_tables.md)
    - Export normalized tables for plotting heatmaps and genes
- [RNA 2](r_files/temporal_rna_2_time_space.md)
    - Calculate RNA diff expression for all pairwise comparisons
- [RNA 3](r_files/temporal_rna_2_time_space_import_plot.md)
    - Fig S1 RNAseq GO enrichment
    - Tables Supp, temporal genes for all domains, temporally dynamic and cell type specific genes.

- [Fig 2 footprints and TF expression]
- [Fig 2 screen results]
- [Fig 2 Eed, Ezh2, Brd8 mutant effects]
- [Fig 3 Nr6a1 mutant effects]
- [Fig 3 Nr6a1 overexpression ATACseq]
