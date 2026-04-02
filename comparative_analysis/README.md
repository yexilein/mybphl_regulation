# GRN inference from multi-tissue single-cell multi-ome datasets

This subdirectory contains the scripts, notebooks, and helper scripts used to identify key Mybphl regulators.

Scripts 01-04 (run in order):
 - 01_make_pseudobulks.R: harmonize ATAC peaks over the three datasets of interest (EOM, limb, larynx), find motifs, compute metacells.
 - 02_run_chromvar.R: compute chromVAR scores for the transcription factor motif found in JASPAR2024.
 - 03_initiate_grn.R: initiate GRN using linear methods (accessibility-corrected Spearman correlation).
 - 04_run_pando.R/sh: prune GRN based on motif presence and Pando’s XGB algorithm.

Notebooks 11-13 (order does not matter):
 - 11_mybp_regulation.Rmd: exploration of regulatory structure around Mybp family genes (based on unpruned GRN).
 - 12_regulon_analysis.Rmd: visualization of key TFs (based on pruned GRN).
 - 13_deg_analysis.Rmd: cross-tissue comparison of expression of key genes.

Helper scripts:
 - aurocs.R: quick computation of AUC values (used for differential expression).
 - colors.R: colors used in preliminary visualizations.
 - dataset.R: helper functions to load and manipulate the datasets.
 - differential_analysis: helper functions for differential expression analysis.
 - grn.R: helper functions for GRN initiation.