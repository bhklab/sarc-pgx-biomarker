# Processed Data Directory

## Purpose

This directory contains **intermediate processed data objects** generated from raw datasets in `data/rawdata/`.  
These files are the result of **preprocessing, quality control, harmonization, subtype annotation, and batch correction** steps.  
They are designed for downstream analyses in the STS biomarker discovery pipeline.

Processed objects are stored as `.qs` (fast R serialization) or `.rda` (R binary data) files, typically containing:
- **Normalized gene expression matrices** (TPM, log2-CPM, or normalized microarray intensities)
- **Drug sensitivity data** (AAC, drug metadata, cell annotations)
- **Clinical metadata** (histology, subtype, treatment, outcome)
- **Subtype mappings** (harmonized labels across datasets)
- Packaged as `SummarizedExperiment` or `MultiAssayExperiment` objects

---

## How to Generate Processed Data

### 1. Download Raw Data

Follow the instructions in [`data/rawdata/README.md`](data/rawdata/README.md)  
and use the acquisition scripts in [`workflow/scripts/README.md`](workflow/scripts/README.md) 


### 2. Expression Processing
- Retained **protein-coding** genes only and removed duplicates
- Harmonized **gene identifiers** to HGNC symbols
- Removed control probes (e.g., AFFX for microarray)
- Ensured feature alignment across GEO and PSets

### 3. Clinical Metadata Curation
- Removed **normal tissue** samples from GEO
- Harmonized histological subtypes:
  - *Leiomyosarcoma*, *Liposarcoma*, *UPS*, *MFH*
- Added `Primary.Metastasis` and `type` fields for integrative modeling
- Preserved original subtype labels in `subtype_original`

### 4. Preclinical Metadata Curation
- Added Cellosaurus-derived subtype annotations
- Mapped tissue lineage (uterus → Soft Tissue)
- Removed “Other” lineage category

### 5. Harmonization & Batch Correction
- Created combined STS matrices for clinical + preclinical data
- Removed non-STS lineages (Bone, Other)
- Retained intersecting feature space (11,049 protein-coding genes)
- Generated STS-only harmonized datasets:
- `PGx_gse_rna_sts.qs` — matched genes, harmonized metadata

---
