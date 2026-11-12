"""
DNA Sequence Analyzer - Streamlit web app
Run locally:  streamlit run streamlit_app.py
Requires fasta_toolkit.py in the same folder.
"""

import io

import numpy as np
import pandas as pd
import streamlit as st
from Bio import Entrez, SeqIO
from Bio.Seq import Seq

from fasta_toolkit import clean_seq, find_orfs, sequence_stats, translate_dna

st.set_page_config(page_title="DNA Sequence Analyzer", page_icon="🧬", layout="wide")
st.title("🧬 DNA Sequence Analyzer & FASTA Toolkit")
st.caption("Inspect sequences, find ORFs in all six frames, and translate DNA to protein.")


# ---------------------------------------------------------------- helpers
@st.cache_data(show_spinner="Downloading from NCBI GenBank...")
def fetch_ncbi(accession: str, email: str) -> str:
    Entrez.email = email
    with Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text") as h:
        return h.read()


def parse_fasta(text: str):
    text = text.strip()
    if text and not text.startswith(">"):
        text = ">pasted_sequence\n" + text
    return list(SeqIO.parse(io.StringIO(text), "fasta"))


@st.cache_data(show_spinner="Scanning six reading frames...")
def cached_orfs(seq_str: str, min_aa: int, table: int):
    return find_orfs(Seq(seq_str), min_aa=min_aa, table=table)


def gc_windows(seq_str: str, window: int):
    arr = np.frombuffer(seq_str.encode(), dtype=np.uint8)
    gc = ((arr == ord("G")) | (arr == ord("C"))).astype(np.int64)
    csum = np.concatenate([[0], np.cumsum(gc)])
    step = max(window // 2, len(arr) // 3000, 1)
    starts = np.arange(0, max(len(arr) - window, 0) + 1, step)
    vals = 100 * (csum[starts + window] - csum[starts]) / window if len(arr) >= window else np.array([])
    return pd.DataFrame({"position": starts + window // 2, "GC %": vals}).set_index("position")


# ---------------------------------------------------------------- sidebar input
with st.sidebar:
    st.header("1. Load a sequence")
    source = st.radio("Source", ["Upload FASTA", "Paste sequence", "Fetch from NCBI"])
    raw = ""
    if source == "Upload FASTA":
        up = st.file_uploader("FASTA file", type=["fasta", "fa", "fna", "txt"])
        if up:
            raw = up.getvalue().decode("utf-8", errors="ignore")
    elif source == "Paste sequence":
        raw = st.text_area("FASTA or raw DNA", height=200)
    else:
        acc = st.text_input("Accession", "NC_045512.2")
        email = st.text_input("Your email (required by NCBI)")
        if st.button("Fetch") and acc and email:
            try:
                st.session_state["ncbi_raw"] = fetch_ncbi(acc, email)
            except Exception as e:
                st.error(f"Download failed: {e}")
        raw = st.session_state.get("ncbi_raw", "")

    st.header("2. Settings")
    table = st.number_input("Genetic code table", 1, 33, 1, help="1 = standard, 11 = bacterial")

if not raw.strip():
    st.info("👈 Load a sequence from the sidebar. Try accession **NC_045512.2** (SARS-CoV-2) via *Fetch from NCBI*.")
    st.stop()

records = parse_fasta(raw)
if not records:
    st.error("No valid FASTA records found.")
    st.stop()

ids = [r.id for r in records]
rec = records[ids.index(st.sidebar.selectbox("Record", ids))] if len(ids) > 1 else records[0]
seq = clean_seq(rec)
seq_str = str(seq)

tab_stats, tab_orf, tab_tr = st.tabs(["📊 Statistics", "🔎 ORF finder", "🔤 Translation"])

# ---------------------------------------------------------------- stats
with tab_stats:
    s = sequence_stats(rec)
    st.subheader(s["id"])
    st.write(s["description"])
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Length (bp)", f"{s['length']:,}")
    c2.metric("GC content", f"{s['gc_percent']:.2f}%")
    c3.metric("GC skew", f"{s['gc_skew']:+.4f}")
    c4.metric("Ambiguous bases", f"{s['ambiguous']:,}")

    left, right = st.columns(2)
    with left:
        st.markdown("**Base composition**")
        st.bar_chart(pd.Series(s["composition"], name="count"))
    with right:
        window = st.slider("GC window (bp)", 100, 5000, 1000, 100)
        st.markdown("**GC content along the sequence**")
        df = gc_windows(seq_str, window)
        if df.empty:
            st.caption("Sequence shorter than the window.")
        else:
            st.line_chart(df)

    rc = seq.reverse_complement()
    st.download_button("Download reverse complement (FASTA)",
                       f">{rec.id}_revcomp\n{rc}\n", file_name=f"{rec.id}_revcomp.fasta")

# ---------------------------------------------------------------- ORFs
with tab_orf:
    min_aa = st.slider("Minimum ORF length (amino acids)", 30, 1000, 100, 10)
    orfs = cached_orfs(seq_str, min_aa, int(table))
    st.write(f"**{len(orfs)}** ORFs ≥ {min_aa} aa (ATG start, six frames)")
    if orfs:
        df = pd.DataFrame(orfs)
        df.insert(0, "#", range(1, len(df) + 1))
        st.dataframe(df.drop(columns="protein"), use_container_width=True, hide_index=True)

        fasta_out = "".join(
            f">{rec.id}_ORF{i} strand={o['strand']} frame={o['frame']} {o['start']}-{o['end']} {o['aa_length']}aa\n"
            f"{o['protein']}\n" for i, o in enumerate(orfs, 1))
        st.download_button("Download ORF proteins (FASTA)", fasta_out, file_name=f"{rec.id}_orfs.faa")

        pick = st.number_input("View protein #", 1, len(orfs), 1)
        o = orfs[pick - 1]
        st.code(o["protein"], language=None)

        # simple map of ORFs
        st.markdown("**ORF length by genome position**")
        st.scatter_chart(df, x="start", y="aa_length", color="strand")
    else:
        st.warning("No ORFs found at this threshold. Try lowering the minimum length.")

# ---------------------------------------------------------------- translate
with tab_tr:
    c1, c2, c3 = st.columns(3)
    start = c1.number_input("Start (1-based)", 1, len(seq), 1)
    end = c2.number_input("End (inclusive)", 1, len(seq), len(seq))
    frame = c3.selectbox("Frame", [1, 2, 3, -1, -2, -3])
    to_stop = st.checkbox("Stop at first stop codon")
    if end < start:
        st.error("End must be ≥ start.")
    else:
        protein = translate_dna(seq[start - 1:end], frame=frame, table=int(table), to_stop=to_stop)
        st.write(f"**{len(protein):,}** amino acids")
        st.code("\n".join(protein[i:i + 60] for i in range(0, len(protein), 60)), language=None)
        st.download_button("Download protein (FASTA)",
                           f">{rec.id}_{start}-{end}_frame{frame}\n{protein}\n",
                           file_name=f"{rec.id}_translation.faa")
