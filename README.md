# Vitamin-D Murine Microbiome and Single-cell Analysis Study

## Overview
This respository contains code and analysis to replicate the single-cell analysis findings in Ryan, et al.
Here, the effets of a Vitamin D-restricted diet on the murine oral microbiome and gingival/buccal transcriptomics were investigated.

## Contents
* The raw_data folder contains the count matrix files generated via BeeNet software package
* vitd_samples.R - contains code for loading matrix files, making seurat objects, and annotating the objects
* buccal_analysis - contains code for analysis of buccal tissue (Figure 4)
* gingival_analysis - contains code for analysis of gingival tissue (Figure 5)
* analysis_buccal_subset - contains code for analysis of major cell types from the buccal tissue (Figure 4, Supplement 4)
* analysis_gingival_subset - contains code for analysis of major cell types from the gingival tissue (Figure 5, Supplement 5)
