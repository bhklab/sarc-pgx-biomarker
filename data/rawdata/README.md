# Raw Data Directory

## Purpose

This directory contains **immutable raw data files** used as original inputs to the biomarker discovery pipeline for soft-tissue sarcoma. These files are **excluded from version control** and must be obtained manually from public repositories to ensure full reproducibility.


---

## Data Access Instructions

**Note:** No raw data files are stored in this repository.  
To reproduce results, download the source datasets as outlined below or use available retrieval scripts (TBD).

---

## Clinical Transcriptomic Data (GEO)

Transcriptomic datasets from soft-tissue sarcoma patients were retrieved from the NCBI Gene Expression Omnibus (GEO):

- [**GSE21122**](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE21122) – 149 STS samples + 9 normal tissues (Affymetrix U133A)
- [**GSE21050**](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE21050) – 310 STS samples with metastasis data (Affymetrix U133 Plus 2.0)
- [**GSE30929**](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE30929) – 140 liposarcomas with DRFS data (Affymetrix U133A)

**Preprocessing:**
- Normalized using Robust Multi-array Average (RMA)
- Gene-level re-annotation using GENCODE v33 for consistent comparison across cohorts

---

## Preclinical Transcriptomic Data 

Microarray data from large-scale pharmacogenomic efforts:

- **CCLE** – HG-U133 Plus 2.0 array (Broad Institute)
- **GDSC** – HG-U219 array (Wellcome Sanger Institute)
- **NCI-Sarcoma** – Human Exon 1.0 ST array (NCI)

**Processing Workflow:**
- Affymetrix CEL files processed using RMA (via `affy` package v1.87.0)
- BrainArray custom CDF annotations (v25.0.0) applied where appropriate
- NCI-Sarcoma data normalized using the `AROMA` pipeline
- Expression harmonized via **quantile normalization**
- Gene symbols mapped from Ensembl IDs using GENCODE v33
- Data represented on log2 RMA scale

---

## Drug Sensitivity Profiles 

Drug response data were obtained from:

- **CCLE**, **GDSC**, **CTRP**, and **NCI-Sarcoma** repositories

**Key Processing:**
- Dose–response curves modeled using a 3-parameter Hill function (`PharmacoGx` v3.13.2)
- Summary metric: **Area Above the Curve (AAC)**
- Quality control applied to exclude artifacts and poor-quality replicates
- Minimum of 10 STS cell lines required per drug for inclusion
- Drug target annotations obtained from **DrugBank v5.1.X** and **ChEMBL** (manual curation for conflicts)

---

## Gene Signature Sets

Curated gene sets used for enrichment analysis:

- [**MSigDB v2025.1** – Hallmark & GO Terms](https://www.gsea-msigdb.org/gsea/msigdb/)
- [**Custom RData Signature File (Zenodo DOI: TBD)**](TBD)

---

## Validation Datasets: TCGA-SARC & NCT-MASTER

Used for external validation of candidate biomarkers:

- **TCGA-SARC (RNA-seq)** – [The Cancer Genome Atlas Sarcoma cohort](https://portal.gdc.cancer.gov/projects/TCGA-SARC)
- **NCT-MASTER (RNA-seq)** – Precision oncology cohort from the DKFZ/NCT/DKTK [ClinicalTrials.gov: NCT05852522](https://clinicaltrials.gov/study/NCT05852522?locStr=Heidelberg,%20Germany&country=Germany&state=Baden-W%C3%BCrttemberg&city=Heidelberg&cond=Cancer&aggFilters=status:rec&rank=2)

---

## Inclusion Criteria

Datasets were included based on:

- Availability of transcriptomic and/or pharmacologic data
- Relevance to soft-tissue sarcoma (histological or lineage annotation)
- Sufficient sample/cell line representation (n ≥ 10 for modeling)
- Compatibility with preprocessing workflows

See the **Materials and Methods** section of the manuscript for full criteria.

---

## Directory Policy

- This directory is **read-only** during pipeline execution
- All downstream outputs are written to `data/procdata/`
- Data versioning and access dates are documented in [`docs/data_sources.md`](/docs/data_sources.md)

---

