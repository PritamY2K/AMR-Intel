import streamlit as st
import pandas as pd
import plotly.express as px
from pathlib import Path
import subprocess
import re
from textwrap import wrap


# =========================
# App configuration
# =========================

st.set_page_config(page_title="AMR-Intel", layout="wide")

BASE_DIR = Path("D:/AMR_App")
UPLOAD_DIR = BASE_DIR / "uploads"
RESULT_DIR = BASE_DIR / "results"
CARD_DB_DIR = BASE_DIR / "card_database"

UPLOAD_DIR.mkdir(parents=True, exist_ok=True)
RESULT_DIR.mkdir(parents=True, exist_ok=True)

st.title("AMR-Intel: Genome-Based Antimicrobial Resistance Prediction App")

st.markdown(
    """
    This app accepts bacterial genome FASTA files and predicts antimicrobial resistance
    using a genome-to-protein workflow:

    **Genome FASTA → FASTA cleaning → Prodigal protein prediction → CARD-RGI protein-mode AMR screening**
    """
)


# =========================
# Helper functions
# =========================

def safe_sample_name(name):
    """
    Convert sample name into a safe filename.
    """
    name = name.strip()
    name = re.sub(r"[^A-Za-z0-9_]+", "_", name)
    return name if name else "sample"


def win_to_wsl(path):
    """
    Convert Windows path D:/AMR_App/... to WSL path /mnt/d/AMR_App/...
    """
    path = str(path).replace("\\", "/")
    path = path.replace("D:", "/mnt/d")
    return path


def clean_genome_fasta(input_path, output_path):
    """
    Clean genome FASTA for downstream Prodigal/RGI use.

    Actions:
    - Removes problematic NCBI headers
    - Writes simple headers: >c1, >c2, etc.
    - Keeps only A, T, G, C, N
    """
    records = []
    current_seq = []

    with open(input_path, "r", errors="ignore") as f:
        for line in f:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if current_seq:
                    records.append("".join(current_seq))
                    current_seq = []
            else:
                cleaned = "".join([c for c in line.upper() if c in "ACGTN"])
                if cleaned:
                    current_seq.append(cleaned)

    if current_seq:
        records.append("".join(current_seq))

    with open(output_path, "w") as f:
        for i, seq in enumerate(records, start=1):
            f.write(f">c{i}\n")
            for part in wrap(seq, 80):
                f.write(part + "\n")

    total_bases = sum(len(seq) for seq in records)

    return len(records), total_bases


def read_rgi_table(file_path):
    """
    Read CARD-RGI tab-separated output.
    """
    return pd.read_csv(file_path, sep="\t")


def summarize_rgi(df):
    """
    Display RGI result summary, charts, and Excel download.
    """
    st.success("AMR result loaded successfully.")

    st.subheader("Raw RGI Result Table")
    st.dataframe(df, width="stretch")

    possible_gene_cols = [
        "Best_Hit_ARO",
        "ARO",
        "AMR Gene Family",
        "AMR_gene_family",
        "Model_name"
    ]

    possible_drug_cols = [
        "Drug Class",
        "Drug_class",
        "drug_class"
    ]

    possible_mech_cols = [
        "Resistance Mechanism",
        "Resistance_mechanism",
        "resistance_mechanism"
    ]

    gene_col = next((c for c in possible_gene_cols if c in df.columns), None)
    drug_col = next((c for c in possible_drug_cols if c in df.columns), None)
    mech_col = next((c for c in possible_mech_cols if c in df.columns), None)

    st.subheader("AMR Summary")

    col1, col2, col3 = st.columns(3)

    with col1:
        if gene_col:
            st.metric("Unique AMR entries", df[gene_col].nunique())
        else:
            st.metric("AMR hits", len(df))

    with col2:
        if drug_col:
            st.metric("Drug classes", df[drug_col].nunique())
        else:
            st.metric("Drug classes", "NA")

    with col3:
        if mech_col:
            st.metric("Resistance mechanisms", df[mech_col].nunique())
        else:
            st.metric("Resistance mechanisms", "NA")

    if drug_col:
        st.subheader("Drug Class Distribution")
        drug_counts = df[drug_col].fillna("Unknown").value_counts().reset_index()
        drug_counts.columns = ["Drug class", "Count"]

        fig = px.bar(
            drug_counts,
            x="Drug class",
            y="Count",
            title="Distribution of predicted AMR drug classes"
        )

        st.plotly_chart(fig, width="stretch")

    if mech_col:
        st.subheader("Resistance Mechanism Distribution")
        mech_counts = df[mech_col].fillna("Unknown").value_counts().reset_index()
        mech_counts.columns = ["Resistance mechanism", "Count"]

        fig2 = px.pie(
            mech_counts,
            names="Resistance mechanism",
            values="Count",
            title="Distribution of predicted resistance mechanisms"
        )

        st.plotly_chart(fig2, width="stretch")

    st.subheader("Download Excel Report")

    report_path = RESULT_DIR / "AMR_Intel_Report.xlsx"
    df.to_excel(report_path, index=False)

    with open(report_path, "rb") as f:
        st.download_button(
            label="Download Excel Report",
            data=f,
            file_name="AMR_Intel_Report.xlsx",
            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
        )


# =========================
# Sidebar/status information
# =========================

st.sidebar.header("AMR-Intel Setup")
st.sidebar.write("Required backend tools:")
st.sidebar.write("1. WSL Ubuntu")
st.sidebar.write("2. Conda environment: rgi_env")
st.sidebar.write("3. Prodigal")
st.sidebar.write("4. CARD-RGI 6.0.5")
st.sidebar.write("5. CARD localDB")


# =========================
# Analysis mode
# =========================

analysis_mode = st.radio(
    "Select analysis mode",
    [
        "Upload genome FASTA and run AMR prediction",
        "Upload existing RGI result file"
    ]
)


# =========================
# Mode 1: Upload existing RGI result
# =========================

if analysis_mode == "Upload existing RGI result file":

    uploaded_result = st.file_uploader(
        "Upload CARD-RGI tab-separated result file",
        type=["txt", "tsv", "csv"]
    )

    if uploaded_result is not None:

        result_path = RESULT_DIR / uploaded_result.name

        with open(result_path, "wb") as f:
            f.write(uploaded_result.getbuffer())

        try:
            if uploaded_result.name.endswith(".csv"):
                df = pd.read_csv(result_path)
            else:
                df = pd.read_csv(result_path, sep="\t")

            summarize_rgi(df)

        except Exception as e:
            st.error(f"Could not read uploaded RGI file: {e}")


# =========================
# Mode 2: Genome FASTA analysis
# =========================

if analysis_mode == "Upload genome FASTA and run AMR prediction":

    st.info(
        """
        This mode runs Prodigal and CARD-RGI through WSL Ubuntu.

        Please ensure that:
        - `rgi_env` exists in WSL
        - RGI version is 6.0.5
        - CARD database is loaded in `D:/AMR_App/card_database/localDB`
        - Prodigal is available inside `rgi_env`
        """
    )

    uploaded_fasta = st.file_uploader(
        "Upload bacterial genome FASTA file",
        type=["fasta", "fa", "fna"]
    )

    sample_name_input = st.text_input("Sample name", value="sample_01")

    run_button = st.button("Run Genome-Based AMR Prediction")

    if uploaded_fasta is not None and run_button:

        sample_name = safe_sample_name(sample_name_input)

        original_fasta = UPLOAD_DIR / f"{sample_name}_original.fasta"
        cleaned_fasta = UPLOAD_DIR / f"{sample_name}_cleaned.fasta"

        protein_faa = RESULT_DIR / f"{sample_name}_proteins.faa"
        protein_clean_faa = RESULT_DIR / f"{sample_name}_proteins_clean.faa"

        rgi_output_prefix = RESULT_DIR / f"{sample_name}_rgi"
        rgi_txt = RESULT_DIR / f"{sample_name}_rgi.txt"

        # Save uploaded genome
        with open(original_fasta, "wb") as f:
            f.write(uploaded_fasta.getbuffer())

        st.success("Genome uploaded successfully.")

        # Clean FASTA
        try:
            record_count, total_bases = clean_genome_fasta(
                original_fasta,
                cleaned_fasta
            )

            st.subheader("Cleaned FASTA Summary")
            st.write(f"Number of sequence records: {record_count}")
            st.write(f"Total nucleotide bases: {total_bases}")

            if record_count == 0 or total_bases == 0:
                st.error("The uploaded file does not contain valid nucleotide FASTA sequence.")
                st.stop()

        except Exception as e:
            st.error(f"FASTA cleaning failed: {e}")
            st.stop()

        # Convert Windows paths to WSL paths
        wsl_cleaned_fasta = win_to_wsl(cleaned_fasta)
        wsl_protein_faa = win_to_wsl(protein_faa)
        wsl_protein_clean_faa = win_to_wsl(protein_clean_faa)
        wsl_rgi_output_prefix = win_to_wsl(rgi_output_prefix)
        wsl_card_db_dir = win_to_wsl(CARD_DB_DIR)

        # Remove previous output files for same sample
        for old_file in RESULT_DIR.glob(f"{sample_name}_rgi*"):
            try:
                old_file.unlink()
            except Exception:
                pass

        # WSL bash command
        bash_command = (
            f"source ~/miniconda3/etc/profile.d/conda.sh && "
            f"conda activate rgi_env && "
            f"prodigal "
            f"-i {wsl_cleaned_fasta} "
            f"-a {wsl_protein_faa} "
            f"-o /dev/null "
            f"-p single && "
            f"tr -d '*' < {wsl_protein_faa} > {wsl_protein_clean_faa} && "
            f"cd {wsl_card_db_dir} && "
            f"rgi main "
            f"-i {wsl_protein_clean_faa} "
            f"-o {wsl_rgi_output_prefix} "
            f"-t protein "
            f"--local "
            f"--clean"
        )

        wsl_command = ["wsl", "bash", "-lc", bash_command]

        st.subheader("Running backend command")
        st.code(" ".join(wsl_command))

        with st.spinner("Running Prodigal and CARD-RGI. Please wait..."):
            result = subprocess.run(
                wsl_command,
                shell=False,
                capture_output=True,
                text=True
            )

        if result.returncode != 0:
            st.error("AMR prediction failed.")

            st.subheader("Error message")
            st.code(result.stderr)

            st.subheader("Command output")
            st.code(result.stdout)

            st.stop()

        st.success("AMR prediction completed successfully.")

        if rgi_txt.exists():
            df = read_rgi_table(rgi_txt)
            summarize_rgi(df)

        else:
            st.error("RGI finished, but the expected RGI .txt file was not found.")

            st.write("Expected file:")
            st.code(str(rgi_txt))

            st.subheader("Command output")
            st.code(result.stdout)

            st.subheader("Error message")
            st.code(result.stderr)