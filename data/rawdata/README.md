# Raw Data Directory

## Purpose

This directory contains **immutable raw data files** used as original inputs to the biomarker discovery pipeline for soft-tissue sarcoma. These files are **excluded from version control** and must be obtained manually from public repositories to ensure full reproducibility.


---

## Data Access Instructions

**Note:** No raw data files are stored in this repository.  
To reproduce results, download the source datasets as outlined below or use available retrieval scripts (TBD).

---

### Clinical Transcriptomic Data (GEO)

Transcriptomic datasets from soft-tissue sarcoma patients were retrieved from the NCBI Gene Expression Omnibus (GEO):

- [**GSE21122**](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE21122) – 149 STS samples + 9 normal tissues (Affymetrix U133A)
- [**GSE21050**](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE21050) – 310 STS samples with metastasis data (Affymetrix U133 Plus 2.0)
- [**GSE30929**](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE30929) – 140 liposarcomas with DRFS data (Affymetrix U133A)

**Preprocessing:**
- Normalized using Robust Multi-array Average (RMA)
- Gene-level re-annotation using GENCODE v33 for consistent comparison across cohorts

---

### Preclinical Transcriptomic Data — CCLE, GDSC, NCI-Sarcoma

Microarray data from large-scale pharmacogenomic efforts:

- **CCLE** – HG-U133 Plus 2.0 array (Broad Institute)
- **GDSC** – HG-U219 array (Wellcome Sanger Institute)
- **NCI-Sarcoma Panel** – Human Exon 1.0 ST array (NCI)

**Processing Workflow:**
- Affymetrix CEL files processed using RMA (via `affy` package v1.87.0)
- BrainArray custom CDF annotations (v25.0.0) applied where appropriate
- NCI-Sarcoma data normalized using the `AROMA` pipeline
- Expression harmonized via **quantile normalization**
- Gene symbols mapped from Ensembl IDs using GENCODE v33
- Data represented on log2 RMA scale

---

### Drug Sensitivity Profiles 

Drug response data were obtained from:

- **CCLE**, **GDSC**, **CTRP**, and **NCI-Sarcoma** repositories

**Key Processing:**
- Dose–response curves modeled using a 3-parameter Hill function (`PharmacoGx` v3.13.2)
- Summary metric: **Area Above the Curve (AAC)**
- Quality control applied to exclude artifacts and poor-quality replicates
- Minimum of 10 STS cell lines required per drug for inclusion
- Drug target annotations obtained from **DrugBank v5.1.X** and **ChEMBL** (manual curation for conflicts)

```

## Signature Sets: Curated Gene Expression Signatures

For downstream enrichment analyses (ORA, GSEA), we used curated gene signatures relevant to drug response and sarcoma biology. These are available from:

- [**Hallmark & GO Terms** — MSigDB v2025.1](https://www.gsea-msigdb.org/gsea/msigdb/)
- [**Precompiled RData File** — Zenodo DOI: TBD](TBD)

---

## Inclusion & Exclusion Criteria

Raw datasets and cell line profiles were selected based on:

- Availability of expression data 
- Relevance to soft-tissue sarcoma (based on histology or tissue label)
- TBD ---> based on validation cohorts: TCGA and MASTER

Refer to the **Materials and Methods** section of the manuscript for full details.

---

## Additional Notes

- GEO datasets were curated to remove samples with missing annotations or technical artifacts.
- Only soft-tissue sarcoma–relevant cell lines were retained from ORCESTRA.
- This directory is **read-only** during pipeline execution. All transformations and analyses are performed downstream in `data/procdata/`.
- For full documentation of data versions and access dates, refer to `docs/data_sources.md`.

