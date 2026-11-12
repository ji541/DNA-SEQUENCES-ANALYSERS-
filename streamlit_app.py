"""Interactive Streamlit interface for DNA sequence exploration."""

from __future__ import annotations

import io

import numpy as np
import pandas as pd
import streamlit as st
from Bio import Entrez, SeqIO
from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from fasta_toolkit import clean_seq, find_orfs, sequence_stats, translate_dna

st.set_page_config(page_title="DNA Sequence Explorer", page_icon="🧬", layout="wide")

st.title("🧬 DNA Sequence Explorer")
st.markdown(
    "Explore FASTA sequences with clear, reproducible summaries. Calculate base composition, "
    "inspect local GC variation, identify **candidate** ORFs, and translate selected regions."
)
st.caption("Educational bioinformatics tool · Results are computational predictions, not experimental validation.")


@st.cache_data(show_spinner="Retrieving sequence from NCBI…")
def fetch_ncbi(accession: str, email: str) -> str:
    Entrez.email = email
    with Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text") as handle:
        data = handle.read()
    if not data.lstrip().startswith(">"):
        raise ValueError("NCBI did not return a FASTA record. Check the accession and try again.")
    return data


def parse_fasta(text: str) -> list[SeqRecord]:
    """Parse FASTA text or wrap a raw nucleotide string as one record."""
    text = text.strip()
    if not text:
        return []
    if not text.startswith(">"):
        text = ">pasted_sequence\n" + text
    records = list(SeqIO.parse(io.StringIO(text), "fasta"))
    for record in records:
        clean_seq(record)
    return records


@st.cache_data(show_spinner="Scanning six reading frames…")
def cached_orfs(sequence: str, min_aa: int, table: int):
    return find_orfs(sequence, min_aa=min_aa, table=table)


def gc_windows(sequence: str, window: int) -> pd.DataFrame:
    """Calculate GC% in overlapping windows; ambiguous bases are excluded."""
    text = sequence.upper()
    if not text or window < 1 or len(text) < window:
        return pd.DataFrame(columns=["GC %"])
    arr = np.frombuffer(text.encode("ascii"), dtype=np.uint8)
    canonical = np.isin(arr, [ord("A"), ord("C"), ord("G"), ord("T")]).astype(np.int64)
    gc = ((arr == ord("G")) | (arr == ord("C"))).astype(np.int64)
    gc_sum = np.concatenate(([0], np.cumsum(gc)))
    base_sum = np.concatenate(([0], np.cumsum(canonical)))
    step = max(window // 2, len(arr) // 3000, 1)
    starts = np.arange(0, len(arr) - window + 1, step)
    gc_counts = gc_sum[starts + window] - gc_sum[starts]
    base_counts = base_sum[starts + window] - base_sum[starts]
    values = np.divide(100 * gc_counts, base_counts, out=np.zeros(len(starts), dtype=float), where=base_counts > 0)
    return pd.DataFrame({"position (bp)": starts + window // 2, "GC %": values}).set_index("position (bp)")


with st.sidebar:
    st.header("1 · Load sequence")
    source = st.radio("Choose a source", ["Upload FASTA", "Paste sequence", "NCBI accession"])
    raw = ""
    if source == "Upload FASTA":
        upload = st.file_uploader("FASTA file", type=["fasta", "fa", "fna", "fas", "txt"])
        if upload:
            raw = upload.getvalue().decode("utf-8-sig", errors="replace")
    elif source == "Paste sequence":
        raw = st.text_area("Paste FASTA or raw DNA/RNA", height=180, placeholder=">sequence_id\nATG...")
    else:
        accession = st.text_input("NCBI nucleotide accession", placeholder="e.g. NC_045512.2")
        email = st.text_input("Contact email (required by NCBI)", placeholder="you@example.org")
        if st.button("Fetch sequence", type="primary"):
            if not accession.strip() or not email.strip():
                st.error("Enter both an accession and contact email.")
            else:
                try:
                    st.session_state["ncbi_fasta"] = fetch_ncbi(accession.strip(), email.strip())
                    st.session_state["ncbi_accession"] = accession.strip()
                except Exception as exc:
                    st.error(f"Could not fetch sequence: {exc}")
        raw = st.session_state.get("ncbi_fasta", "")

    st.divider()
    st.header("2 · Analysis settings")
    table_ids = sorted(CodonTable.unambiguous_dna_by_id)
    table_names = {number: CodonTable.unambiguous_dna_by_id[number].names[0] for number in table_ids}
    selected_table = st.selectbox(
        "Genetic code",
        table_ids,
        index=table_ids.index(1),
        format_func=lambda number: f"{number} · {table_names[number]}",
        help="Choose the NCBI translation table appropriate for the organism/sequence.",
    )

if not raw.strip():
    st.info("Load a sequence from the sidebar to begin. For a first example, try an NCBI accession such as **NC_045512.2**.")
    st.stop()

try:
    records = parse_fasta(raw)
except (ValueError, UnicodeError) as exc:
    st.error(f"Sequence input could not be parsed: {exc}")
    st.stop()

if not records:
    st.error("No FASTA records were found. Check the input and try again.")
    st.stop()

if len(records) > 1:
    selected_id = st.selectbox("Select FASTA record", range(len(records)), format_func=lambda i: f"{records[i].id} — {len(records[i].seq):,} bp")
else:
    selected_id = 0
record = records[selected_id]
sequence = clean_seq(record)
sequence_text = str(sequence)
if not sequence_text:
    st.warning("The selected FASTA record contains no sequence bases.")
    st.stop()

with st.expander("Sequence metadata", expanded=False):
    st.write(f"**Record:** `{record.id}`")
    st.write(f"**Description:** {record.description or 'No description provided'}")
    st.write(f"**Input:** {len(records)} FASTA record(s) · selected record has {len(sequence):,} bp")

stats_tab, composition_tab, orf_tab, translation_tab = st.tabs(
    ["📊 Overview", "🧪 Composition & GC", "🔎 ORF candidates", "🔤 Translation"]
)

with stats_tab:
    summary = sequence_stats(record)
    st.subheader("Sequence summary")
    col1, col2, col3, col4 = st.columns(4)
    col1.metric("Length", f"{summary['length']:,} bp")
    col2.metric("GC content", f"{summary['gc_percent']:.2f}%", help="Calculated among unambiguous A/C/G/T bases.")
    col3.metric("GC skew", f"{summary['gc_skew']:+.4f}", help="(G − C) / (G + C); descriptive and not diagnostic by itself.")
    col4.metric("Ambiguous bases", f"{summary['ambiguous']:,}")
    st.write(f"AT content among unambiguous bases: **{summary['at_percent']:.2f}%**")
    st.caption("GC/AT percentages use only canonical A, C, G, and T bases; ambiguous symbols are counted separately.")
    reverse = SeqRecord(sequence.reverse_complement(), id=f"{record.id}_revcomp", description="reverse complement")
    rc_buffer = io.StringIO()
    SeqIO.write(reverse, rc_buffer, "fasta")
    st.download_button("Download reverse complement (FASTA)", rc_buffer.getvalue(), file_name=f"{record.id}_revcomp.fasta", mime="text/plain")

with composition_tab:
    st.subheader("Nucleotide composition")
    composition = pd.DataFrame({"Count": summary["composition"], "Percent of A/C/G/T": summary["composition_percent"]})
    left, right = st.columns(2)
    with left:
        st.bar_chart(composition["Count"], color="#19a974")
    with right:
        st.dataframe(composition.style.format({"Count": "{:,.0f}", "Percent of A/C/G/T": "{:.2f}%"}), use_container_width=True)
    max_window = min(5000, len(sequence))
    if max_window < 1:
        st.info("No bases available for a GC window plot.")
    else:
        default_window = min(500, max_window)
        window = st.slider("Sliding-window width (bp)", min_value=1, max_value=max_window, value=default_window)
        gc_data = gc_windows(sequence_text, window)
        if gc_data.empty:
            st.info("The selected sequence is shorter than the window.")
        else:
            st.markdown("**Local GC content**")
            st.line_chart(gc_data, y="GC %", color="#4c78a8")
            st.caption("Windows overlap by approximately 50%. Ambiguous bases are excluded from each window's denominator.")

with orf_tab:
    st.subheader("Six-frame ORF candidate search")
    st.warning("ORFs are sequence-based candidates, not proof of gene expression or function. This simple method uses ATG starts and the selected genetic code's stop codons; it does not model splice sites, promoters, or organism-specific biology.")
    min_aa = st.slider("Minimum ORF length (amino acids; stop codon excluded)", 1, min(1000, max(1, len(sequence) // 3)), min(100, max(1, len(sequence) // 3)), 1)
    with st.spinner("Searching forward and reverse-complement frames…"):
        orfs = cached_orfs(sequence_text, min_aa, selected_table)
    st.metric("Candidates found", len(orfs))
    if orfs:
        orf_frame = pd.DataFrame(orfs)
        display = orf_frame.drop(columns=["protein"]).copy()
        st.dataframe(display, use_container_width=True, hide_index=True)
        fasta_text = "".join(
            f">{record.id}_ORF{i} strand={hit['strand']} frame={hit['frame']} "
            f"coordinates={hit['start']}-{hit['end']} length={hit['aa_length']}aa\n{hit['protein']}\n"
            for i, hit in enumerate(orfs, 1)
        )
        st.download_button("Download candidate proteins (FASTA)", fasta_text, file_name=f"{record.id}_orfs.faa", mime="text/plain")
        st.markdown("**Candidate length by start coordinate**")
        st.scatter_chart(display, x="start", y="aa_length", color="strand")
        selected_orf = st.selectbox("Inspect predicted protein", range(len(orfs)), format_func=lambda i: f"ORF {i + 1}: {orfs[i]['strand']} strand, {orfs[i]['start']}–{orfs[i]['end']} ({orfs[i]['aa_length']} aa)")
        st.code("\n".join(orfs[selected_orf]["protein"][i:i + 60] for i in range(0, len(orfs[selected_orf]["protein"]), 60)), language=None)
    else:
        st.info("No ORF candidates met this threshold. Try a shorter minimum length or another genetic code.")

with translation_tab:
    st.subheader("Translate a selected region")
    region_col1, region_col2, region_col3 = st.columns(3)
    start = region_col1.number_input("Start coordinate (1-based)", min_value=1, max_value=len(sequence), value=1)
    end = region_col2.number_input("End coordinate (inclusive)", min_value=1, max_value=len(sequence), value=len(sequence))
    frame = region_col3.selectbox("Reading frame", [1, 2, 3, -1, -2, -3], format_func=lambda value: f"{value:+d}" if value > 0 else str(value))
    to_stop = st.checkbox("Stop translation at the first in-frame stop codon")
    if end < start:
        st.error("End coordinate must be greater than or equal to start coordinate.")
    else:
        protein = translate_dna(sequence[start - 1:end], frame=frame, table=selected_table, to_stop=to_stop)
        st.write(f"**Translated length:** {len(protein):,} amino acid(s)")
        if protein:
            st.code("\n".join(protein[i:i + 60] for i in range(0, len(protein), 60)), language=None)
        else:
            st.info("The selected region does not contain a complete codon in this frame.")
        protein_record = SeqRecord(Seq(protein), id=f"{record.id}_{start}-{end}_frame{frame}", description=f"genetic_code={selected_table}")
        protein_buffer = io.StringIO()
        SeqIO.write(protein_record, protein_buffer, "fasta")
        st.download_button("Download translated protein (FASTA)", protein_buffer.getvalue(), file_name=f"{record.id}_translation.faa", mime="text/plain")
