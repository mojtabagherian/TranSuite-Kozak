# TranSuite-Kozak

Scripts to evaluate **Kozak consensus enrichment (RNNATGGV)** in start codons re-annotated by TranSuite.  
This repository accompanies our manuscript and provides reproducible code to check whether revised ORFs are supported by stronger Kozak sequences than the original Araport11 reference.

---

## 📖 Overview
- Extracts the −3..+4 sequence context around annotated AUG codons.  
- Scores strict Kozak matches using the plant-friendly definition **RNNATGGV** (R = A/G, V = A/C/G).  
- Compares Reference vs Revised annotations for transcripts where ORFs change.  
- Outputs per-transcript results as a CSV file.

---

## 🛠 Requirements
- Python 3.7+  
- Standard libraries only (`argparse`, `csv`, `re`, `collections`)  

No additional packages are required.

---

## ▶️ Usage

### 1. Inputs
- **Transcriptome FASTA** (`Athaliana_transcripts.fa`)  
- **Reference GTF** (`Athaliana_447_Araport11_Reference_ORF.gtf`)  
- **Revised GTF** (`Transfix.gtf`)  
- **dest-compare.csv** with at least 3 columns:
  - `Transcript_ID`
  - `is_PTC50nt_ref`
  - `is_PTC50nt_revise`

Example `dest-compare.csv`:
```csv
Transcript_ID,is_PTC50nt_ref,is_PTC50nt_revise
AT1G01020.3,False,True
AT2G03740.2,True,False

Example Run code:

  python3 kozak.py \
    --compare_csv dest-compare.csv \
    --ref_gtf Athaliana_447_Araport11_Reference_ORF.gtf \
    --rev_gtf Transfix.gtf \
    --tx_fasta Athaliana_transcripts.fa \
    --tx_col Transcript_ID \
    --ref_col is_PTC50nt_ref \
    --revise_col is_PTC50nt_revise \
    --out_prefix kozak_changed



## 📊 Output

The script writes a CSV (e.g. `kozak_changed_results.csv`) containing:

| transcript_id | ref_start | ref_8mer  | ref_kozak | rev_start | rev_8mer  | rev_kozak |
|---------------|-----------|-----------|-----------|-----------|-----------|-----------|
| AT1G01020.3   | 583       | CAGATGAG  | False     | 465       | GCAATGGC  | True      |


MIT 



