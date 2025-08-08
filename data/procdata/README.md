# Processed Data Directory

## Purpose

This directory contains **intermediate processed data objects** generated from raw datasets in `data/rawdata/`.  
The files here are the result of **preprocessing, quality control, harmonization, subtype annotation, and batch correction** steps. They are designed for downstream analyses in the STS biomarker discovery pipeline.

Processed objects are stored as `.qs` or `.rda` files, typically containing:
- **Normalized gene expression matrices** (TPM or normalized microarray intensities)
- **Drug sensitivity data** (AAC, drug meta-data, cell information)
- **Clinical metadata** (histology, subtype, treatment, outcome)
- **Subtype mappings** (harmonized labels across datasets)
- Packaged as `SummarizedExperiment` or `MultiAssayExperiment` objects

---

## How to Generate Processed Data

To prepare the  `.qs` or `.rda` files from raw inputs, follow these steps:

### 1. Download Raw Data

Download the raw dataset from ORCESTRA or GEO as explained in `data/rawdata/` and follow scripts TBD.

### 2. Expression Processing

- Retained **protein-coding** genes only and removed duplicated gene names
- Harmonized **gene identifiers** to HGNC symbols

### 3. Clinical Metadata Curation

- Removed **normal tissue** samples from GEO
- Harmonized histological subtypes:
  - *Leiomyosarcoma*, *Liposarcoma*, *UPS*, *MFH*
- Added `Primary.Metastasis` and `type` fields for integrative modeling
- Retained original subtype labels in `subtype_original`

### 4. Preclinical Metadata Curation

- Added cellosaurus annotations for subtype mapping
- TBD

### 5. Batch Correction
- Created combined STS matrices for clinical and preclinical data
- Separated **Soft Tissue** from **Bone** lineages for downstream analysis
- Prepared STS-only datasets:
  - `PGx_gse_rna_sts.qs` — matched genes, harmonized metadata

---
