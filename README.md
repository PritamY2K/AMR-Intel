\# AMR-Intel



AMR-Intel is a Streamlit-based antimicrobial resistance prediction app.



\## Main workflow



Genome FASTA upload → FASTA cleaning → Prodigal protein prediction → CARD-RGI protein-mode AMR screening → AMR summary dashboard.



\## Features



\- Upload bacterial genome FASTA files

\- Clean genome FASTA headers automatically

\- Predict proteins using Prodigal

\- Run CARD-RGI in protein mode

\- Summarize AMR genes, drug classes, and resistance mechanisms

\- Generate downloadable Excel reports



\## Local setup



This app requires:



\- Python

\- Streamlit

\- WSL Ubuntu

\- Conda environment named `rgi\_env`

\- CARD-RGI

\- Prodigal

\- CARD database loaded locally



\## Run locally



```bash

streamlit run app.py

