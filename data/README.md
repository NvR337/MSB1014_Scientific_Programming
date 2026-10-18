# Data

All input data come from NASA Open Science Data Repository (OSDR) study **OSD-245**
(Rodent Research-6, mouse liver, bulk RNA-seq):
https://osdr.nasa.gov/bio/repo/data/studies/OSD-245

Download the following files from the study's *Files* tab and place them in this folder
**with these exact names** (the notebook reads them from `data/`):

| File | Section on OSDR | Used for |
|---|---|---|
| `GLDS-245_rna_seq_RSEM_Unnormalized_Counts_GLbulkRNAseq.csv` | GeneLab processed RNA-seq | Raw (unnormalized) gene counts, genes × samples |
| `s_OSD-245.txt` | Study metadata (ISA-Tab) | Sample metadata / experimental factors |
| `a_OSD-245_transcription-profiling_rna-sequencing-(rna-seq)_Illumina NovaSeq 6000.txt` | Study metadata (ISA-Tab) | Assay-level metadata |
| `GLDS-245_rna_seq_qc_metrics_GLbulkRNAseq.csv` | GeneLab processed RNA-seq | Per-sample QC metrics (read depth, RIN, mapping rates, …) |

Data accessed on: **18 October 2026**.

<!-- If you committed the files themselves, replace the instructions above with:
"The files are included in this folder as downloaded from OSDR on <date>." -->
