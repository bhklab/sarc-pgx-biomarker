# Processed Data Directory

## Purpose

This directory contains **intermediate processed data objects** generated from raw inputs located in `data/rawdata/`.  
These objects represent the outputs of standardized **preprocessing, quality control, harmonization, subtype curation**, and **batch correction** steps.

These processed files are formatted for downstream use in the **soft-tissue sarcoma (STS) biomarker discovery pipeline**.

---

## Contents

Data are saved as `.qs` (fast R serialization) or `.rda` files and typically packaged as `SummarizedExperiment` or `MultiAssayExperiment` objects. They include:

- **Normalized gene expression matrices**  
  - RMA-normalized (microarrays) or log2 TPM (RNA-seq)
  - Harmonized to HGNC gene symbols via GENCODE v33
  - Cross-platform alignment using quantile normalization

- **Drug sensitivity profiles**  
  - AAC values from 3-parameter Hill curve modeling (`PharmacoGx`)
  - Matched drug and cell metadata
  - Quality-controlled: drugs with >60% missing values excluded

- **Clinical and preclinical metadata**  
  - Subtype, histology, treatment status, outcome data
  - Cellosaurus-based lineage annotations
  - Harmonized across datasets and platforms

- **Subtype harmonization maps**  
  - Mapping between original and curated subtype labels
  - Used for cross-cohort integration and stratification

---

## How to Generate Processed Data

### 1. Download Raw Data

Follow instructions in [`data/rawdata/README.md`](/rawdata/README.md)  
Use acquisition scripts from [`workflow/scripts/`](/workflow/scripts)

---

### 2. Expression Processing

- Retained **protein-coding genes** only (n = 11,049)
- Gene-level reannotation via **GENCODE v33**
- Control probes (e.g., AFFX) and duplicates removed
- Harmonized identifiers to **HGNC symbols**
- Platforms:
  - Affymetrix U133A / Plus 2.0 → `affy` + BrainArray CDF v25.0.0
  - Exon 1.0 ST → AROMA workflow
  - RNA-seq (TCGA, NCT-MASTER) → log2(TPM + 1)

---

### 3. Clinical Metadata Curation

- Removed **normal tissue** controls (e.g., from GSE21122)
- Harmonized subtype labels across datasets:
  - e.g., *Leiomyosarcoma*, *Liposarcoma*, *UPS*, *MFH*
- Added fields:
  - `Primary.Metastasis` — binary metastatic status
  - `type` — clinical vs. preclinical
  - `subtype_original` — preserved raw annotations

---

### 4. Preclinical Metadata Curation

- Cell line annotations from **CCLE**, **GDSC**, **NCI-Sarcoma**
- Lineage derived from **Cellosaurus**, with manual curation
- Removed ambiguous/“Other” tissue origins
- Standardized drug identifiers and targets using **DrugBank** and **ChEMBL**

---

### 5. Cross-Dataset Harmonization

- Expression matrices merged across platforms and sources
- Applied **quantile normalization** to align distributions
- Batch effects minimized for integrative modeling
- Filtered to retain:
  - Only STS-relevant samples/cell lines
  - Intersecting gene feature set (11,049 protein-coding genes)

**Key processed object:**
- `PGx_gse_rna_sts.qs` — harmonized clinical + preclinical expression data with metadata

---

## Additional Notes

- All processed files in this directory are **derived artifacts** and not meant to be edited manually.
- For transparency, each processing step is logged in the pipeline outputs.
- Data versions and transformations are documented in [`docs/data_sources.md`](/docs/data_sources.md)

---