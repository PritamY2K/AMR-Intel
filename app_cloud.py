import streamlit as st
import pandas as pd
import plotly.express as px
from pathlib import Path
import subprocess
import re
from textwrap import wrap
import shutil
import zipfile


# =========================
# App configuration
# =========================

st.set_page_config(page_title="AMR-Intel Cloud", layout="wide")

BASE_DIR = Path(__file__).resolve().parent
UPLOAD_DIR = BASE_DIR / "uploads"
RESULT_DIR = BASE_DIR / "results"
CARD_DB_DIR = BASE_DIR / "card_database"

UPLOAD_DIR.mkdir(parents=True, exist_ok=True)
RESULT_DIR.mkdir(parents=True, exist_ok=True)

st.title("AMR-Intel Cloud: Integrated Genome-Based AMR Prediction App")

st.markdown(
    """
    **AMR-Intel** predicts antimicrobial resistance from bacterial genome data using
    multiple AMR detection engines.

    **Workflow:**  
    Genome FASTA / NCBI accession / WGS reads → FASTA preparation → Prodigal protein prediction →
    CARD-RGI → AMRFinderPlus → ABRicate-ResFinder → Comparative AMR dashboard.

    This cloud version runs directly on Linux/Docker and does not require WSL.
    """
)


# =========================
# Basic helper functions
# =========================

def safe_sample_name(name):
    """Convert a user-provided sample name into a safe file-name stem."""
    name = str(name).strip()
    name = re.sub(r"[^A-Za-z0-9_]+", "_", name)
    return name if name else "sample"


def win_to_wsl(path):
    """Convert Windows D:/... paths to WSL /mnt/d/... paths."""
    path = str(path).replace("\\", "/")
    path = path.replace("D:", "/mnt/d")
    return path


def run_wsl_command(bash_command):
    """Run one bash command inside WSL and return result + command list."""
    wsl_command = ["wsl", "bash", "-lc", bash_command]
    result = subprocess.run(
        wsl_command,
        shell=False,
        capture_output=True,
        text=True
    )
    return result, wsl_command


def clean_genome_fasta(input_path, output_path):
    """
    Clean genome FASTA for downstream Prodigal/RGI/ABRicate use.

    Actions:
    - Removes problematic original FASTA headers
    - Writes simple headers: >c1, >c2, etc.
    - Keeps only A, T, G, C, N characters
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


def normalize_colname(col):
    """Normalize a column name for flexible matching."""
    return re.sub(r"[^a-z0-9]+", "", str(col).lower())


def find_first_column(df, possible_columns):
    """
    Flexible column finder.
    First tries exact matching, then normalized matching.
    """
    for col in possible_columns:
        if col in df.columns:
            return col

    normalized_map = {normalize_colname(c): c for c in df.columns}

    for col in possible_columns:
        key = normalize_colname(col)
        if key in normalized_map:
            return normalized_map[key]

    return None


def clean_gene_name(value):
    """Standardize gene names for comparison."""
    if pd.isna(value):
        return "Unknown"

    value = str(value).strip()

    if value == "":
        return "Unknown"

    value = value.replace("allele", "").strip()
    value = re.sub(r"\s+", " ", value)

    return value


# =========================
# Readers
# =========================

def read_rgi_table(file_path):
    return pd.read_csv(file_path, sep="\t")


def read_amrfinder_table(file_path):
    return pd.read_csv(file_path, sep="\t")


def read_abricate_table(file_path):
    return pd.read_csv(file_path, sep="\t")


# =========================
# Standardization functions
# =========================

def standardize_rgi(df):
    gene_col = find_first_column(
        df,
        [
            "Best_Hit_ARO",
            "Best Hit ARO",
            "ARO",
            "ARO Name",
            "ARO_name",
            "Model_name",
            "Model Name",
            "AMR Gene Family",
            "AMR_gene_family"
        ]
    )

    drug_col = find_first_column(
        df,
        [
            "Drug Class",
            "Drug_class",
            "drug_class"
        ]
    )

    mech_col = find_first_column(
        df,
        [
            "Resistance Mechanism",
            "Resistance_mechanism",
            "resistance_mechanism"
        ]
    )

    family_col = find_first_column(
        df,
        [
            "AMR Gene Family",
            "AMR_gene_family",
            "Model_type",
            "Model Type"
        ]
    )

    out = pd.DataFrame(index=df.index)

    if gene_col:
        out["gene"] = df[gene_col].apply(clean_gene_name)
    else:
        out["gene"] = "Unknown"

    if drug_col:
        out["drug_class"] = df[drug_col].fillna("Unknown").astype(str)
    else:
        out["drug_class"] = "Unknown"

    if mech_col:
        out["mechanism"] = df[mech_col].fillna("Unknown").astype(str)
    else:
        out["mechanism"] = "Unknown"

    if family_col:
        out["gene_family"] = df[family_col].fillna("Unknown").astype(str)
    else:
        out["gene_family"] = "Unknown"

    out["tool"] = "CARD-RGI"
    out["present"] = 1

    return out[["tool", "gene", "drug_class", "mechanism", "gene_family", "present"]]


def standardize_amrfinder(df):
    gene_col = find_first_column(
        df,
        [
            "Gene symbol",
            "Gene Symbol",
            "gene_symbol",
            "Element symbol",
            "Element Symbol",
            "Protein name",
            "Protein Name",
            "Name"
        ]
    )

    class_col = find_first_column(
        df,
        [
            "Class",
            "class",
            "Drug Class",
            "drug_class"
        ]
    )

    subclass_col = find_first_column(
        df,
        [
            "Subclass",
            "subclass"
        ]
    )

    element_type_col = find_first_column(
        df,
        [
            "Element type",
            "Element Type",
            "element_type",
            "Scope",
            "scope"
        ]
    )

    method_col = find_first_column(
        df,
        [
            "Method",
            "method"
        ]
    )

    out = pd.DataFrame(index=df.index)

    if gene_col:
        out["gene"] = df[gene_col].apply(clean_gene_name)
    else:
        out["gene"] = "Unknown"

    if class_col:
        out["drug_class"] = df[class_col].fillna("Unknown").astype(str)
    else:
        out["drug_class"] = "Unknown"

    if element_type_col:
        out["mechanism"] = df[element_type_col].fillna("Unknown").astype(str)
    else:
        out["mechanism"] = "AMR determinant"

    if subclass_col:
        out["gene_family"] = df[subclass_col].fillna("Unknown").astype(str)
    elif method_col:
        out["gene_family"] = df[method_col].fillna("Unknown").astype(str)
    else:
        out["gene_family"] = "Unknown"

    out["tool"] = "AMRFinderPlus"
    out["present"] = 1

    return out[["tool", "gene", "drug_class", "mechanism", "gene_family", "present"]]


def standardize_abricate(df):
    gene_col = find_first_column(
        df,
        [
            "GENE",
            "#GENE",
            "Gene",
            "gene"
        ]
    )

    product_col = find_first_column(
        df,
        [
            "PRODUCT",
            "Product",
            "product"
        ]
    )

    resistance_col = find_first_column(
        df,
        [
            "RESISTANCE",
            "Resistance",
            "resistance"
        ]
    )

    db_col = find_first_column(
        df,
        [
            "DATABASE",
            "Database",
            "db"
        ]
    )

    out = pd.DataFrame(index=df.index)

    if gene_col:
        out["gene"] = df[gene_col].apply(clean_gene_name)
    else:
        out["gene"] = "Unknown"

    if resistance_col:
        out["drug_class"] = df[resistance_col].fillna("Unknown").astype(str)
    else:
        out["drug_class"] = "Unknown"

    if product_col:
        out["mechanism"] = df[product_col].fillna("AMR determinant").astype(str)
    else:
        out["mechanism"] = "AMR determinant"

    if db_col:
        out["gene_family"] = df[db_col].fillna("ResFinder").astype(str)
    else:
        out["gene_family"] = "ResFinder"

    out["tool"] = "ABRicate-ResFinder"
    out["present"] = 1

    return out[["tool", "gene", "drug_class", "mechanism", "gene_family", "present"]]


# =========================
# Comparative analysis
# =========================

def create_combined_df(rgi_df=None, amrfinder_df=None, abricate_df=None):
    frames = []

    if rgi_df is not None and not rgi_df.empty:
        frames.append(standardize_rgi(rgi_df))

    if amrfinder_df is not None and not amrfinder_df.empty:
        frames.append(standardize_amrfinder(amrfinder_df))

    if abricate_df is not None and not abricate_df.empty:
        frames.append(standardize_abricate(abricate_df))

    if not frames:
        return pd.DataFrame(columns=[
            "tool", "gene", "drug_class", "mechanism", "gene_family", "present"
        ])

    combined_df = pd.concat(frames, ignore_index=True)

    combined_df["gene"] = combined_df["gene"].fillna("Unknown").astype(str)
    combined_df["tool"] = combined_df["tool"].fillna("Unknown").astype(str)
    combined_df["drug_class"] = combined_df["drug_class"].fillna("Unknown").astype(str)
    combined_df["mechanism"] = combined_df["mechanism"].fillna("Unknown").astype(str)
    combined_df["gene_family"] = combined_df["gene_family"].fillna("Unknown").astype(str)

    combined_df = combined_df[combined_df["gene"].str.strip() != ""]
    combined_df = combined_df[combined_df["gene"] != "Unknown"]

    return combined_df


def create_consensus_table(combined_df):
    if combined_df.empty:
        return pd.DataFrame(columns=["gene", "Detected tools", "Consensus confidence"])

    matrix = combined_df.pivot_table(
        index="gene",
        columns="tool",
        values="present",
        aggfunc="max",
        fill_value=0
    ).reset_index()

    tool_columns = [c for c in matrix.columns if c != "gene"]

    for col in tool_columns:
        matrix[col] = matrix[col].astype(int)

    def assign_confidence(row):
        detected_count = sum(row[col] > 0 for col in tool_columns)

        if detected_count >= 3:
            return "High"
        elif detected_count == 2:
            return "Moderate"
        elif detected_count == 1:
            return "Tool-specific / Low"
        else:
            return "Not detected"

    matrix["Detected tools"] = matrix[tool_columns].sum(axis=1)
    matrix["Consensus confidence"] = matrix.apply(assign_confidence, axis=1)

    for col in tool_columns:
        matrix[col] = matrix[col].map({1: "Yes", 0: "No"})

    return matrix


def create_tool_count_table(combined_df):
    if combined_df.empty:
        return pd.DataFrame(columns=["Tool", "Detected AMR entries"])

    return (
        combined_df.groupby("tool")["gene"]
        .nunique()
        .reset_index()
        .rename(columns={"tool": "Tool", "gene": "Detected AMR entries"})
    )


def create_drug_class_table(combined_df):
    if combined_df.empty:
        return pd.DataFrame(columns=["Tool", "Drug class", "Detected genes"])

    return (
        combined_df.groupby(["tool", "drug_class"])["gene"]
        .nunique()
        .reset_index()
        .rename(columns={
            "tool": "Tool",
            "drug_class": "Drug class",
            "gene": "Detected genes"
        })
    )


def create_mechanism_table(combined_df):
    if combined_df.empty:
        return pd.DataFrame(columns=["Tool", "Mechanism", "Detected genes"])

    return (
        combined_df.groupby(["tool", "mechanism"])["gene"]
        .nunique()
        .reset_index()
        .rename(columns={
            "tool": "Tool",
            "mechanism": "Mechanism",
            "gene": "Detected genes"
        })
    )


# =========================
# Display functions
# =========================

def show_individual_summary(raw_df, tool_name):
    st.subheader(f"{tool_name} Raw Result Table")
    st.dataframe(raw_df, width="stretch")

    if tool_name == "CARD-RGI":
        std_df = standardize_rgi(raw_df)
    elif tool_name == "AMRFinderPlus":
        std_df = standardize_amrfinder(raw_df)
    elif tool_name == "ABRicate-ResFinder":
        std_df = standardize_abricate(raw_df)
    else:
        std_df = pd.DataFrame()

    if std_df.empty:
        st.warning("No standardized AMR entries could be extracted.")
        return

    col1, col2, col3 = st.columns(3)

    with col1:
        st.metric(f"{tool_name} AMR genes", std_df["gene"].nunique())

    with col2:
        st.metric(f"{tool_name} drug classes", std_df["drug_class"].nunique())

    with col3:
        st.metric(f"{tool_name} mechanisms", std_df["mechanism"].nunique())

    drug_counts = (
        std_df.groupby("drug_class")["gene"]
        .nunique()
        .reset_index()
        .rename(columns={"drug_class": "Drug class", "gene": "Detected genes"})
    )

    if not drug_counts.empty:
        fig = px.bar(
            drug_counts,
            x="Drug class",
            y="Detected genes",
            title=f"{tool_name}: Drug-class distribution"
        )
        st.plotly_chart(fig, width="stretch")


def show_comparative_dashboard(rgi_df, amrfinder_df, abricate_df, sample_name):
    with st.expander("Debug: raw output dimensions and column names"):
        if rgi_df is not None:
            st.write("CARD-RGI shape:", rgi_df.shape)
            st.write("CARD-RGI columns:", list(rgi_df.columns))
        else:
            st.write("CARD-RGI: None")

        if amrfinder_df is not None:
            st.write("AMRFinderPlus shape:", amrfinder_df.shape)
            st.write("AMRFinderPlus columns:", list(amrfinder_df.columns))
        else:
            st.write("AMRFinderPlus: None")

        if abricate_df is not None:
            st.write("ABRicate shape:", abricate_df.shape)
            st.write("ABRicate columns:", list(abricate_df.columns))
        else:
            st.write("ABRicate: None")

    combined_df = create_combined_df(rgi_df, amrfinder_df, abricate_df)
    consensus_df = create_consensus_table(combined_df)
    tool_counts = create_tool_count_table(combined_df)
    drug_class_df = create_drug_class_table(combined_df)
    mechanism_df = create_mechanism_table(combined_df)

    st.subheader("Combined Standardized AMR Evidence Table")
    st.dataframe(combined_df, width="stretch")

    st.subheader("Consensus Gene-Level Comparison")
    st.dataframe(consensus_df, width="stretch")

    high_count = 0
    moderate_count = 0
    low_count = 0

    if not consensus_df.empty and "Consensus confidence" in consensus_df.columns:
        high_count = (consensus_df["Consensus confidence"] == "High").sum()
        moderate_count = (consensus_df["Consensus confidence"] == "Moderate").sum()
        low_count = (consensus_df["Consensus confidence"] == "Tool-specific / Low").sum()

    col1, col2, col3, col4 = st.columns(4)

    with col1:
        st.metric(
            "Total unique AMR genes",
            consensus_df["gene"].nunique() if not consensus_df.empty and "gene" in consensus_df.columns else 0
        )

    with col2:
        st.metric("High-confidence genes", high_count)

    with col3:
        st.metric("Moderate-confidence genes", moderate_count)

    with col4:
        st.metric("Tool-specific genes", low_count)

    st.subheader("Tool-wise AMR Gene Count")

    if not tool_counts.empty:
        fig_tool = px.bar(
            tool_counts,
            x="Tool",
            y="Detected AMR entries",
            title="Number of unique AMR genes detected by each tool"
        )
        st.plotly_chart(fig_tool, width="stretch")
        st.dataframe(tool_counts, width="stretch")
    else:
        st.info("No tool-wise AMR gene count available.")

    st.subheader("Drug-Class Comparison")

    if not drug_class_df.empty:
        fig_drug = px.bar(
            drug_class_df,
            x="Drug class",
            y="Detected genes",
            color="Tool",
            barmode="group",
            title="Drug-class distribution across AMR tools"
        )
        st.plotly_chart(fig_drug, width="stretch")
        st.dataframe(drug_class_df, width="stretch")
    else:
        st.info("No drug-class comparison available.")

    st.subheader("Mechanism / Element-Type Comparison")

    if not mechanism_df.empty:
        fig_mech = px.bar(
            mechanism_df,
            x="Mechanism",
            y="Detected genes",
            color="Tool",
            barmode="group",
            title="Mechanism or element-type distribution across AMR tools"
        )
        st.plotly_chart(fig_mech, width="stretch")
        st.dataframe(mechanism_df, width="stretch")
    else:
        st.info("No mechanism comparison available.")

    st.subheader("Gene Presence/Absence Heatmap")

    if not combined_df.empty:
        presence_matrix = combined_df.pivot_table(
            index="gene",
            columns="tool",
            values="present",
            aggfunc="max",
            fill_value=0
        )

        if not presence_matrix.empty:
            fig_heatmap = px.imshow(
                presence_matrix.T,
                x=presence_matrix.index,
                y=presence_matrix.columns,
                aspect="auto",
                title="AMR gene presence/absence heatmap"
            )
            st.plotly_chart(fig_heatmap, width="stretch")
    else:
        st.info("No gene presence/absence heatmap available.")

    st.subheader("Download Combined Excel Report")

    report_path = RESULT_DIR / f"{sample_name}_AMR_Comparative_Report.xlsx"

    with pd.ExcelWriter(report_path, engine="openpyxl") as writer:
        if rgi_df is not None:
            rgi_df.to_excel(writer, sheet_name="CARD_RGI_raw", index=False)
        if amrfinder_df is not None:
            amrfinder_df.to_excel(writer, sheet_name="AMRFinderPlus_raw", index=False)
        if abricate_df is not None:
            abricate_df.to_excel(writer, sheet_name="ABRicate_raw", index=False)

        combined_df.to_excel(writer, sheet_name="Combined_standardized", index=False)
        consensus_df.to_excel(writer, sheet_name="Consensus_comparison", index=False)
        tool_counts.to_excel(writer, sheet_name="Tool_counts", index=False)
        drug_class_df.to_excel(writer, sheet_name="Drug_class_summary", index=False)
        mechanism_df.to_excel(writer, sheet_name="Mechanism_summary", index=False)

    with open(report_path, "rb") as f:
        st.download_button(
            label="Download Comparative Excel Report",
            data=f,
            file_name=f"{sample_name}_AMR_Comparative_Report.xlsx",
            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
        )


# =========================
# Core AMR pipeline
# =========================

def run_full_amr_pipeline(input_fasta, sample_name):
    sample_name = safe_sample_name(sample_name)

    original_fasta = Path(input_fasta)
    cleaned_fasta = UPLOAD_DIR / f"{sample_name}_cleaned.fasta"

    protein_faa = RESULT_DIR / f"{sample_name}_proteins.faa"
    protein_clean_faa = RESULT_DIR / f"{sample_name}_proteins_clean.faa"

    rgi_output_prefix = RESULT_DIR / f"{sample_name}_rgi"
    rgi_txt = RESULT_DIR / f"{sample_name}_rgi.txt"

    amrfinder_output = RESULT_DIR / f"{sample_name}_amrfinder.tsv"
    abricate_output = RESULT_DIR / f"{sample_name}_abricate_resfinder.tsv"

    try:
        record_count, total_bases = clean_genome_fasta(
            original_fasta,
            cleaned_fasta
        )

        st.subheader("Cleaned FASTA Summary")
        st.write(f"Number of sequence records: {record_count}")
        st.write(f"Total nucleotide bases: {total_bases}")

        if record_count == 0 or total_bases == 0:
            st.error("The selected file does not contain valid nucleotide FASTA sequence.")
            st.stop()

    except Exception as e:
        st.error(f"FASTA cleaning failed: {e}")
        st.stop()

    wsl_cleaned_fasta = win_to_wsl(cleaned_fasta)
    wsl_protein_faa = win_to_wsl(protein_faa)
    wsl_protein_clean_faa = win_to_wsl(protein_clean_faa)
    wsl_rgi_output_prefix = win_to_wsl(rgi_output_prefix)
    wsl_card_db_dir = win_to_wsl(CARD_DB_DIR)
    wsl_amrfinder_output = win_to_wsl(amrfinder_output)
    wsl_abricate_output = win_to_wsl(abricate_output)

    for pattern in [
        f"{sample_name}_rgi*",
        f"{sample_name}_amrfinder*",
        f"{sample_name}_abricate*",
        f"{sample_name}_proteins*"
    ]:
        for old_file in RESULT_DIR.glob(pattern):
            try:
                old_file.unlink()
            except Exception:
                pass

    bash_command = (
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
        f"--clean && "

        f"amrfinder "
        f"-p {wsl_protein_clean_faa} "
        f"-o {wsl_amrfinder_output} && "

        f"abricate "
        f"--db resfinder "
        f"{wsl_cleaned_fasta} "
        f"> {wsl_abricate_output}"
    )

    st.subheader("Running backend command")
    st.code(bash_command)

    with st.spinner("Running Prodigal, CARD-RGI, AMRFinderPlus, and ABRicate. Please wait..."):
        result, wsl_command = run_wsl_command(bash_command)

    if result.returncode != 0:
        st.error("AMR prediction failed.")
        st.subheader("Error message")
        st.code(result.stderr)
        st.subheader("Command output")
        st.code(result.stdout)
        st.stop()

    st.success("AMR prediction completed successfully.")

    rgi_df = None
    amrfinder_df = None
    abricate_df = None

    if rgi_txt.exists():
        rgi_df = read_rgi_table(rgi_txt)
    else:
        st.warning("CARD-RGI output file was not found.")

    if amrfinder_output.exists():
        amrfinder_df = read_amrfinder_table(amrfinder_output)
    else:
        st.warning("AMRFinderPlus output file was not found.")

    if abricate_output.exists():
        try:
            abricate_df = read_abricate_table(abricate_output)
        except Exception:
            abricate_df = pd.DataFrame()
            st.warning("ABRicate output file was found but could not be parsed.")
    else:
        st.warning("ABRicate output file was not found.")

    tabs = st.tabs([
        "CARD-RGI",
        "AMRFinderPlus",
        "ABRicate",
        "Comparative Summary"
    ])

    with tabs[0]:
        if rgi_df is not None:
            show_individual_summary(rgi_df, "CARD-RGI")
        else:
            st.warning("No CARD-RGI result available.")

    with tabs[1]:
        if amrfinder_df is not None:
            show_individual_summary(amrfinder_df, "AMRFinderPlus")
        else:
            st.warning("No AMRFinderPlus result available.")

    with tabs[2]:
        if abricate_df is not None:
            show_individual_summary(abricate_df, "ABRicate-ResFinder")
        else:
            st.warning("No ABRicate result available.")

    with tabs[3]:
        show_comparative_dashboard(rgi_df, amrfinder_df, abricate_df, sample_name)


# =========================
# NCBI accession downloader
# =========================

def download_ncbi_genome(accession, sample_name):
    sample_name = safe_sample_name(sample_name)

    accession_dir = UPLOAD_DIR / f"{sample_name}_ncbi_download"
    accession_zip = UPLOAD_DIR / f"{sample_name}_ncbi.zip"

    if accession_dir.exists():
        shutil.rmtree(accession_dir)

    if accession_zip.exists():
        accession_zip.unlink()

    accession_dir.mkdir(parents=True, exist_ok=True)

    wsl_accession_zip = win_to_wsl(accession_zip)

    bash_command = (
        f"datasets download genome accession {accession} "
        f"--include genome "
        f"--filename {wsl_accession_zip}"
    )

    st.subheader("NCBI genome download command")
    st.code(bash_command)

    with st.spinner("Downloading genome from NCBI..."):
        result, wsl_command = run_wsl_command(bash_command)

    if result.returncode != 0:
        st.error("NCBI genome download failed.")
        st.subheader("Error message")
        st.code(result.stderr)
        st.subheader("Command output")
        st.code(result.stdout)
        st.stop()

    if not accession_zip.exists():
        st.error("NCBI download finished, but ZIP file was not found.")
        st.stop()

    try:
        with zipfile.ZipFile(accession_zip, "r") as zip_ref:
            zip_ref.extractall(accession_dir)
    except Exception as e:
        st.error(f"Could not extract NCBI ZIP file: {e}")
        st.stop()

    fna_files = list(accession_dir.rglob("*.fna"))

    if not fna_files:
        st.error("No .fna genome file was found after NCBI download.")
        st.stop()

    selected_fasta = fna_files[0]

    st.success(f"Genome downloaded successfully: {selected_fasta.name}")

    return selected_fasta


# =========================
# WGS assembly through Shovill
# =========================

def assemble_wgs_reads(r1_path, r2_path, sample_name):
    sample_name = safe_sample_name(sample_name)

    assembly_dir = RESULT_DIR / f"{sample_name}_shovill"
    contigs_path = assembly_dir / "contigs.fa"

    if assembly_dir.exists():
        shutil.rmtree(assembly_dir)

    wsl_r1 = win_to_wsl(r1_path)
    wsl_r2 = win_to_wsl(r2_path)
    wsl_assembly_dir = win_to_wsl(assembly_dir)

    bash_command = (
        f"shovill "
        f"--R1 {wsl_r1} "
        f"--R2 {wsl_r2} "
        f"--outdir {wsl_assembly_dir} "
        f"--force"
    )

    st.subheader("WGS assembly command")
    st.code(bash_command)

    with st.spinner("Assembling WGS reads using Shovill. This may take time..."):
        result, wsl_command = run_wsl_command(bash_command)

    if result.returncode != 0:
        st.error("WGS assembly failed.")
        st.subheader("Error message")
        st.code(result.stderr)
        st.subheader("Command output")
        st.code(result.stdout)
        st.stop()

    if not contigs_path.exists():
        st.error("Assembly completed, but contigs.fa was not found.")
        st.stop()

    st.success("WGS assembly completed successfully.")
    return contigs_path


# =========================
# Sidebar/status information
# =========================

st.sidebar.header("AMR-Intel Cloud Setup")
st.sidebar.write("Required tools inside the cloud container/Linux server:")
st.sidebar.write("1. Prodigal, CARD-RGI, AMRFinderPlus, NCBI Datasets CLI and Shovill")
st.sidebar.write("2. ABRicate")
st.sidebar.write("3. CARD database loaded in `card_database/localDB`")
st.sidebar.write("4. AMRFinderPlus database updated with `amrfinder -u`")
st.sidebar.write("5. ABRicate databases prepared with `abricate --setupdb`")


# =========================
# Analysis mode
# =========================

analysis_mode = st.radio(
    "Select analysis mode",
    [
        "Upload genome FASTA and run AMR prediction",
        "Enter NCBI genome accession and run AMR prediction",
        "Upload WGS FASTQ reads and run AMR prediction",
        "Upload existing RGI, AMRFinderPlus and ABRicate result files",
        "Upload existing RGI result file only"
    ]
)


# =========================
# Mode 1: Upload genome FASTA
# =========================

if analysis_mode == "Upload genome FASTA and run AMR prediction":

    st.info(
        """
        Upload a bacterial genome FASTA file. The app will run Prodigal, CARD-RGI,
        AMRFinderPlus and ABRicate-ResFinder, followed by comparative visualization.
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

        with open(original_fasta, "wb") as f:
            f.write(uploaded_fasta.getbuffer())

        st.success("Genome uploaded successfully.")

        run_full_amr_pipeline(original_fasta, sample_name)


# =========================
# Mode 2: NCBI accession
# =========================

if analysis_mode == "Enter NCBI genome accession and run AMR prediction":

    st.info(
        """
        Enter an NCBI assembly accession such as `GCF_000006945.2` or `GCA_...`.
        The app will download the genome FASTA using NCBI Datasets CLI and run AMR analysis.
        """
    )

    accession = st.text_input(
        "Enter NCBI assembly accession",
        value="GCF_000006945.2"
    )

    sample_name_input = st.text_input("Sample name", value="ncbi_sample")

    run_accession = st.button("Download Genome and Run AMR Prediction")

    if run_accession:
        sample_name = safe_sample_name(sample_name_input)
        fasta_path = download_ncbi_genome(accession, sample_name)
        run_full_amr_pipeline(fasta_path, sample_name)


# =========================
# Mode 3: WGS FASTQ input
# =========================

if analysis_mode == "Upload WGS FASTQ reads and run AMR prediction":

    st.warning(
        """
        This mode is for paired-end WGS reads. It assembles reads using Shovill,
        then runs the same AMR pipeline on assembled contigs.
        This can take substantial time for large datasets.
        """
    )

    uploaded_r1 = st.file_uploader(
        "Upload R1 FASTQ / FASTQ.GZ",
        type=["fastq", "fq", "gz"],
        key="r1_fastq"
    )

    uploaded_r2 = st.file_uploader(
        "Upload R2 FASTQ / FASTQ.GZ",
        type=["fastq", "fq", "gz"],
        key="r2_fastq"
    )

    sample_name_input = st.text_input("Sample name", value="wgs_sample")

    run_wgs = st.button("Assemble WGS Reads and Run AMR Prediction")

    if uploaded_r1 is not None and uploaded_r2 is not None and run_wgs:
        sample_name = safe_sample_name(sample_name_input)

        r1_path = UPLOAD_DIR / f"{sample_name}_R1.fastq.gz"
        r2_path = UPLOAD_DIR / f"{sample_name}_R2.fastq.gz"

        with open(r1_path, "wb") as f:
            f.write(uploaded_r1.getbuffer())

        with open(r2_path, "wb") as f:
            f.write(uploaded_r2.getbuffer())

        st.success("WGS FASTQ files uploaded successfully.")

        contigs_path = assemble_wgs_reads(r1_path, r2_path, sample_name)

        run_full_amr_pipeline(contigs_path, sample_name)


# =========================
# Mode 4: Upload existing outputs
# =========================

if analysis_mode == "Upload existing RGI, AMRFinderPlus and ABRicate result files":

    st.subheader("Upload Existing AMR Tool Outputs")

    uploaded_rgi = st.file_uploader(
        "Upload CARD-RGI tab-separated result file",
        type=["txt", "tsv", "csv"],
        key="uploaded_rgi_existing"
    )

    uploaded_amrfinder = st.file_uploader(
        "Upload AMRFinderPlus tab-separated result file",
        type=["txt", "tsv", "csv"],
        key="uploaded_amrfinder_existing"
    )

    uploaded_abricate = st.file_uploader(
        "Upload ABRicate-ResFinder tab-separated result file",
        type=["txt", "tsv", "csv"],
        key="uploaded_abricate_existing"
    )

    sample_name_existing = st.text_input("Sample name", value="sample_existing")

    if uploaded_rgi is not None and uploaded_amrfinder is not None and uploaded_abricate is not None:

        sample_name_existing = safe_sample_name(sample_name_existing)

        rgi_path = RESULT_DIR / uploaded_rgi.name
        amrfinder_path = RESULT_DIR / uploaded_amrfinder.name
        abricate_path = RESULT_DIR / uploaded_abricate.name

        with open(rgi_path, "wb") as f:
            f.write(uploaded_rgi.getbuffer())

        with open(amrfinder_path, "wb") as f:
            f.write(uploaded_amrfinder.getbuffer())

        with open(abricate_path, "wb") as f:
            f.write(uploaded_abricate.getbuffer())

        try:
            rgi_df = read_rgi_table(rgi_path)
            amrfinder_df = read_amrfinder_table(amrfinder_path)
            abricate_df = read_abricate_table(abricate_path)

            tabs = st.tabs([
                "CARD-RGI",
                "AMRFinderPlus",
                "ABRicate",
                "Comparative Summary"
            ])

            with tabs[0]:
                show_individual_summary(rgi_df, "CARD-RGI")

            with tabs[1]:
                show_individual_summary(amrfinder_df, "AMRFinderPlus")

            with tabs[2]:
                show_individual_summary(abricate_df, "ABRicate-ResFinder")

            with tabs[3]:
                show_comparative_dashboard(
                    rgi_df,
                    amrfinder_df,
                    abricate_df,
                    sample_name_existing
                )

        except Exception as e:
            st.error(f"Could not process uploaded result files: {e}")


# =========================
# Mode 5: Upload existing RGI only
# =========================

if analysis_mode == "Upload existing RGI result file only":

    uploaded_result = st.file_uploader(
        "Upload CARD-RGI tab-separated result file",
        type=["txt", "tsv", "csv"],
        key="uploaded_rgi_only"
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

            show_individual_summary(df, "CARD-RGI")

        except Exception as e:
            st.error(f"Could not read uploaded RGI file: {e}")
