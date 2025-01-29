# The cis-regulatory logic integrating spatial and temporal patterning in the vertebrate neural tube

Code associated with the manuscript [Zhang et al (2024)](https://www.biorxiv.org/content/10.1101/2024.04.17.589864v1.full).

## Data availability

- [GEO superseries](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE264172)

## Data processing

### 1. ATAC-seq processing 
Using [nf-core/atacseq](https://nf-co.re/atacseq)
```
nextflow run nf-core/atacseq \
	--input design.csv \
        --genome mm10 \
    	--skip_diff_analysis \
    	--min_reps_consensus 2 \
        -r 1.2.1
```

### 2. RNA-seq processing
Using [nf-core/rnaseq](https://nf-co.re/rnaseq)
```
nextflow run nf-core/rnaseq \
        --input design.csv \
        --genome mm10 \
        --aligner star_salmon \
        --star_index false \
        --skip_biotype_qc \
    	--deseq2_vst \
        -r 3.5
```

### 3. Footprinting
Using [TOBIAS nextflow pipeline](https://github.com/luslab/briscoe-nf-tobias)
```
nextflow run luslab/briscoe-nf-tobias \
  -r master \
  --skip_bam_index true \
  --design design.csv \
  --genome genome.fa \
  --regions results/bwa/mergedLibrary/macs/broadPeak/consensus/consensus_peaks.mLb.clN.bed \
  --peaks results/bwa/mergedLibrary/macs/broadPeak/consensus/consensus_peaks.mLb.clN.bed \
  --blacklist ~/.nextflow/assets/nf-core/atacseq/assets/blacklists/mm10-blacklist.bed \
  --motifs motifs_archetypes.meme \
  --bd_hmem true \
  --bd_cpus 64 \
```

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
- [RNA 2a](r_files/temporal_rna_2_time_space.md)
    - Calculate RNA diff expression for all pairwise comparisons
- [RNA 2b](r_files/temporal_rna_2_time_space_import_plot.md)
    - Fig S1 RNAseq GO enrichment
    - Tables Supp, temporal genes for all domains, temporally dynamic and cell type specific genes
- [RNA 3](r_files/temporal_rna_3_plotgenes_import.md)
    - Gene expression plot for various figures
- [Footprinting 1](r_files/temporal_footprint_1_WT.md)
    - Fig 2 and Fig S1 temporal footprints and correlated TFs
- [Screen](r_files/temporal_screen_1_hitselection.md)
    - Fig 2 and S2: screen results and hit highlight
- [Screen hits flow analysis](r_files/temporal_screen_3_Flow.md)
    - Fig S3: proportion of NFIA progenitors after PRC or Brd8 perturbation
    - Fig 5: proportion of NFIA progenitors after Nr6a1 perturbation
    - Fig 6: proportion of NFIA neurons after Nr6a1 perturbation
- [Mutant RNAseq](r_files/temporal_screen_2_hitsRNAseq.md)
    - Fig 3 and S3 Polycomb and Brd8 mutants: temporal and spatial patterning defects
    - Fig 4 and S4 Nr6a1 mutant: temporal defects
- [Nr6a1 OE ATAC 1](r_files/temporal_Nr6a1OE_2_ATAC.md)
    - Fig 4 Accessibility changes upon Nr6a1 overexpression
- [Nr6a1 OE ATAC 2](r_files/temporal_Nr6a1OE_3_ATAC.md)
    - Fig 4 Accessibility in wild type CaTS-ATAC for elements changing in Nr6a1 overexpression
- [Time and space integration](r_files/temporal_multi_2_time_and_space_genes.md)
    - Fig 5 and S5 cell type specific and temporally dynamic genes and associated elements
- [NFIA/B KO ATAC](r_files/temporal_NFIAB-KO_1_ATAC.md)
    - Fig 5 heatmap
- [NIFA/B KO ATAC motifs](r_files/temporal_NFIAB-KO_4_plotmotifs.md)
    - Fig 5 motif enrichment plotting 
- [NFIA/B KO RNA](r_files/temporal_NFIAB-KO_2_RNA.md)
    - Fig 5 bar plots 
