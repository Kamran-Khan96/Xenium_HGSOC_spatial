# Xenium_HGSOC_spatial

This repository contains representative Python and R scripts used in the analyses accompanying our study of **122 treatment-naïve high-grade serous ovarian cancer (HGSOC) samples** profiled using the **10x Genomics Xenium in situ spatial transcriptomics platform** with a targeted 480-gene panel.

The study investigated associations between **tumour-intrinsic gene expression, the spatial composition of tumour-neighbouring cell populations, and overall survival (OS)**.

The code provided here illustrates the key analytical workflows used throughout the study, including spatial data processing, transcript reassignment, cell-type annotation, tumour pseudobulk generation, spatial neighbourhood analysis, survival modelling, and downstream statistical analyses. These scripts are provided as representative examples to support reproducibility and facilitate adaptation to similar spatial transcriptomic datasets.

**Publication:** [Link to the published paper]

## Repository contents

1. **Python scripts** – Spatial data processing, transcript reassignment, quality control, cell-type annotation, spatial neighbourhood analyses, and cell-level analyses.
2. **R scripts** – Statistical analyses, survival modelling, risk-score generation, visualisation, and downstream analyses.

## Folder and script descriptions

### 1. Xenium data processing and transcript reassignment with MisTIC

**`1.MisTIC_on_xenium_TMA6`**
Adjusts potentially misassigned transcripts using MisTIC, starting from the raw Xenium outputs. Code for one representative Xenium TMA dataset (TMA6) is provided.

### 2. Preprocessing, QC and supervised annotation

**`2.1.TMA6 tma meta info incorporation`**
Prepares and incorporates metadata describing individual TMA cores into the Xenium spatial data.

**`2.2.QC and Supervised annotation`**
Performs quality control of the Xenium count matrix and supervised cell-type annotation.

**`2.3.Resolvi for high-dimensional tSNE`**
Generates a high-dimensional t-SNE representation of merged TMA5 and TMA6 data using the latent representation obtained from scVI/resolVI.

### 3. Tumour-cell pseudobulk generation and cell-type abundance calculation

**`3.1. getting relative cell type abundances and tumour cell pseudobulk expression`**
Calculates patient-level relative cell-type abundances and generates tumour-cell pseudobulk expression profiles.

**`3.2. batch adjustments of tumour pseudobulk data with ComBat`**
Adjusts TMA-specific technical batch effects in the tumour-cell pseudobulk expression data using ComBat.

### 4. Tumour-centric spatial neighbourhood calculation

**`4. Estimating the abundances of tumour neighbouring non-tumour cell-types`**
Uses Squidpy to quantify the relative abundance of non-tumour cell populations neighbouring tumour cells across patients and at incremental spatial distances.

### 5. Survival analysis

**`5.1. univariate Cox regression tumour intrinsic genes`**
Performs univariable Cox regression to identify tumour-intrinsic genes associated with overall survival.

**`5.2. prioritisation of survival associated genes with elastic net regularisation`**
Uses elastic-net regularisation to prioritise survival-associated tumour-intrinsic genes.

**`5.3. ssGSEA gene set scoring with elastic net prioritised genes`**
Calculates ssGSEA scores for gene sets comprising elastic-net-prioritised tumour-intrinsic genes associated with shorter or longer survival, and evaluates their associations with overall survival.

**`5.4. unified risk scoring with EN prioritised genes`**
Calculates a coefficient-weighted unified risk score using the elastic-net-prioritised genes.

**`5.5. Screening for survival association of tumour neighbouring cell type distribution`**
Evaluates associations between tumour-neighbouring cell-type distributions and overall survival across incremental spatial distances.

**`5.6. joint survival association of survival-associated genesets and cell-type distribution`**
Evaluates joint associations between survival-associated tumour-intrinsic gene sets and tumour-neighbouring cell populations.

## Notes

1. This repository contains **representative analysis scripts** rather than a complete end-to-end analysis pipeline.
2. Patient-level clinical and molecular data are not included because of ethical, privacy, and data-access restrictions.
3. The scripts may require modification to accommodate different Xenium datasets, cell-type annotation schemes, file structures, or computing environments.
4. File paths and computational environment settings may need to be adapted before running the scripts.
5. The analyses shown in the repository correspond to the workflows described in the accompanying publication.

## Citation

If you use these scripts or adapt the workflows for your research, please cite the accompanying publication:

[biorxiv]
