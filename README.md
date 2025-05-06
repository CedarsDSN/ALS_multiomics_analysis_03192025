**ALS Multiomics Project**
==========================

# Project Description
This project aims to identify and characterize ALS molecular subtypes by integrating transcriptomics, proteomics, epigenomics, and clinical data from the Answer ALS dataset. Patients will be grouped based on ALSFRS-R progression rates, and multi-omics differential analysis will be used to extract subtype-defining features. Clustering, including consensus clustering, will be used to identify stable subtypes, which will then be evaluated for associations with clinical outcomes. If subtypes lack clear clinical separation, re-clustering with clinical variables will be performed. Polygenic risk scores (PRS) will be computed using GWAS data and tested for association with key genes and proteins. eQTL and pQTL analyses will assess whether these molecular features are genetically regulated. Pathway enrichment will provide a biological interpretation of each subtype’s molecular signature.

# Data
AnswerALS contains clinical and omics data (genomics, transcriptomics, epigenomics, proteomics) and is available upon providing information using a research proposal and signing a data user agreement. Data was requested from AnswerALS on March 27th, 2025. Data user agreement was received on 9th April, 2025. Access to data was granted on 28th April, 2025. 

# Analysis
04/29/2025 - The first step in analysis is categorizing the patients into fast, slow, and intermediate progressors, using the clinical data made available. More information about the steps, files, and corresponding columns used for calculating the progression rate is available in the readme file under analysis folder. The file containing the categorization of patients is under analysis, named "ALS_Progression_Rates_Full.csv". The script used to generate this is under scripts folder "calculate_progression.py".

05/06/2025 - This update compares two approaches for identifying differentially expressed genes (DEGs) across ALS progression subgroups: Fast, Intermediate, and Slow.

Approach 1: Standard DESeq2 Pipeline
Used DESeq2’s default model (DESeq()), which includes built-in filtering and dispersion estimation.

Designed for quick, exploratory analyses.

Detected a modest number of DEGs, with the largest signal observed in the Fast vs Intermediate comparison.

Approach 2: Manual Dispersion-Aware Workflow
Applied stricter gene filtering: kept only genes expressed (≥10 normalized counts) in at least 4 samples.

Performed dispersion estimation and Wald testing manually with increased iteration limits.

Included only genes with stable model convergence (betaConv = TRUE).

Produced fewer DEGs but with higher statistical confidence.

There are very few significant DE genes at both levels. I will now explore dividing the patients into fast and slow progressors only and check if more differentially expressed genes show up. For a more detailed summary, refer to the README section in analysis folder.

# Results

