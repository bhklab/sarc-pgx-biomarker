# Raw Data Directory

## Purpose

This directory is reserved for **immutable raw data files** that serve as the original input to the biomarker discovery pipeline for soft-tissue sarcoma. These files are not tracked by Git and must be obtained separately to ensure full reproducibility.

---

## Data Access Instructions

**No raw data files are included in this repository.**  
To reproduce the results, you must manually download the original data using the links below or follow the shared scripts.

### Clinical Datasets: Gene Expression Omnibus (GEO)

Curated gene expression microarray datasets for clinical soft-tissue sarcoma samples were obtained from the GEO database:

- [**GSE21122**](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE21122)
- [**GSE21050**](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE21050)
- [**GSE30929**](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE30929)

Each dataset includes:
- Microarray-based gene expression data
- Sample metadata (e.g., histology, subtype, time-to-evetn outcomes)
- Downloaded using `GEOquery` and curated using custom scripts TBD

---

### Cell Line Pharmacogenomics Data: ORCESTRA Platform

Pharmacogenomic profiles for soft-tissue sarcoma cell lines were retrieved from [**ORCESTRA**](https://www.orcestra.ca/pset):

This dataset includes:
- Gene expression microarray data for sarcoma cell lines
- Drug sensitivity profiles (e.g., AAC, IC50, AUC)
- Curated metadata for lineage annotation and quality control

---

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

