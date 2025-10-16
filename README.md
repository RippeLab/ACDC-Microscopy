# ACDC-Microscopy
Collection of scripts used for Preprocessing and Analysis of padFISH and Immunofluorescence (IF) data sets in the study "Two different chromatin modules regulate proinflammatory gene expression" from Seufert et al. (BioRxiv) (https://doi.org/10.1101/2024.08.03.606159)


# padFISH - padlock probe-based smFISH
padFISH is a multiplexed single-molecule fluorescence in situ hybridization (smFISH) method that combines elements from hybridization-based in situ sequencing (HybISS) [1] and single-cell resolution in situ hybridization on tissues (SCRINSHOT) [2]. In this study, we employed padFISH to trace nascent RNAs using intronic padlock probes (PLPs) targeting their cDNA and to detect their rolling circle amplification (RCA) product through subsequent Bridge Probe (BP) and Detection Oligo (DO) hybridization. padFISH data was used for (1) the co-expression analysis at the CXCL gene cluster and (2) the bursting kinetics-related analysis at selected genes.

# Immunofluorescence staining (IF)
Immunostaining was performed to detect the endogenous protein NFkB/p65 before and after TNF-alpha induction. 

# Preprocessing
1.	For the preprocessing of raw images (Imaris format) and metadata, image stacks were first transformed into maximum projected TIF files by running the RunIJ_imstoTif_GPU_v3.R  script and imstoTif_headless_v2.ijm FIJI macro.
2.	The 02_Preprocessing_v1.6.ijm FIJI macro performed flatfield correction, chromatic aberration correction and stitching (via Grid/Collection Stitching FIJI plugin). Stitched images were used as input in all subsequent analysis.
3.	Nuclei segmentation on DAPI images was performed with Cellpose 2 [3] using the pretrained cyto model with diameter 150 and 200 for 60x and 100x objectives, respectively. To run Cellpose in batch the RunCellpose_batch_GPU_v3.R was used.
4.	Cell nuclei at the borders of the image or that displayed overexposure in individual channels were filtered out in R prior to further analysis. See IF_Integrated_intensity_v5.R and padFISH_CXCL_Coexpression_analysis_v2.R scripts.

# Analysis
1.	Individual channel images and Cellpose nuclear as well as calculated cytoplasmic masks were used directly in FIJI or using the quantNuclei function in R to quantify image features (e. g. integrated or mean intensity, area). See see IF_Integrated_intensity_v5.R and padFISH_CXCL_Coexpression_analysis_v2.R for the usage of custom function quantNuclei_v01.R [4]. 
3.	For padFISH-based co-expression analysis at the CXCL cluster: the padFISH_CXCL_Coexpression_analysis_v2.R script was first used to process padFISH data, then a combined analysis with sequencing data was performed (see Zenodo repository).
4.  For Immunostaining of NFkB/p65, we quantified the density distribution of the signal with IF_Integrated_intensity_v5.R script and the ratios of signal with the IF_CytoNucRatio_v2.R script.
5.  For padFISH-based bursting kinetics analysis, the padFISH_BurstingKinetics_v1.R script was used. 

# References 
1.	Gyllborg, D. et al. Hybridization-based in situ sequencing (HybISS) for spatially resolved transcriptomics in human and mouse brain tissue. Nucleic Acids Res 48, e112 (2020)
2.	Sountoulidis, A. et al. SCRINSHOT enables spatial mapping of cell states in tissue sections with single-cell resolution. PLoS Biol 18, e3000675 (2020)
3.	Pachitariu, M. & Stringer, C. Cellpose 2.0: how to train your own model. Nat Methods 19, 1634-1641 (2022)
4.	Frank, L. Perturbing and imaging nuclear compartments to reveal mechanisms of transcription regulation and telomere maintenance. PhD thesis, Heidelberg University (2023)
