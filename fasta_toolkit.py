#!/usr/bin/env python3
"""
DNA Sequence Analyzer & FASTA Toolkit
=====================================
Inspect a FASTA file, find open reading frames (ORFs) on all six frames,
and translate DNA into amino acids, using Biopython.

Setup
-----
    pip install biopython

Usage
-----
    # 1. Download SARS-CoV-2 reference genome from NCBI (needs internet + your email)
    python fasta_toolkit.py download NC_045512.2 --email you@example.com -o NC_045512.2.fasta

    # 2. Basic sequence statistics
    python fasta_toolkit.py stats NC_045512.2.fasta

    # 3. Find ORFs of at least 100 amino acids (all 6 frames), save proteins to FASTA
    python fasta_toolkit.py orfs NC_045512.2.fasta --min-aa 100 --proteins-out orfs.faa

    # 4. Translate: whole sequence in frame 1, or a region (1-based, inclusive)
    python fasta_toolkit.py translate NC_045512.2.fasta --frame 1
    python fasta_toolkit.py translate NC_045512.2.fasta --start 21563 --end 25384 --to-stop

    # 5. Reverse complement
    python fasta_toolkit.py revcomp NC_045512.2.fasta -o revcomp.fasta
"""

import argparse
import sys
from collections import Counter

from Bio import Entrez, SeqIO
from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# --------------------------------------------------------------------------
# File I/O
# --------------------------------------------------------------------------
def download_fasta(accession: str, email: str, out_path: str) -> str:
    """Download a nucleotide record from NCBI GenBank as FASTA."""
    Entrez.email = email  # NCBI requires a contact email
    with Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text") as handle:
        data = handle.read()
    if not data.startswith(">"):
        sys.exit(f"Unexpected response from NCBI for '{accession}':\n{data[:200]}")
    with open(out_path, "w") as fh:
        fh.write(data)
    return out_path


def load_records(path: str) -> list:
    """Parse a FASTA file and return all records (validates the format)."""
    records = list(SeqIO.parse(path, "fasta"))
    if not records:
        sys.exit(f"No FASTA records found in '{path}'. Is it a valid FASTA file?")
    return records


def clean_seq(record: SeqRecord) -> Seq:
    """Uppercase the sequence and convert RNA 'U' to 'T'."""
    return Seq(str(record.seq).upper().replace("U", "T"))


# --------------------------------------------------------------------------
# Sequence statistics
# --------------------------------------------------------------------------
def gc_content(seq: Seq) -> float:
    """GC percentage over unambiguous bases (A, C, G, T)."""
    s = str(seq)
    acgt = sum(s.count(b) for b in "ACGT")
    return 100.0 * (s.count("G") + s.count("C")) / acgt if acgt else 0.0


def sequence_stats(record: SeqRecord) -> dict:
    seq = clean_seq(record)
    counts = Counter(str(seq))
    length = len(seq)
    g, c = counts.get("G", 0), counts.get("C", 0)
    return {
        "id": record.id,
        "description": record.description,
        "length": length,
        "gc_percent": gc_content(seq),
        "at_percent": 100.0 - gc_content(seq) if length else 0.0,
        "gc_skew": (g - c) / (g + c) if (g + c) else 0.0,
        "composition": {b: counts.get(b, 0) for b in "ACGT"},
        "ambiguous": length - sum(counts.get(b, 0) for b in "ACGT"),
        "starts_with_atg": str(seq).startswith("ATG"),
    }


# --------------------------------------------------------------------------
# ORF finding
# --------------------------------------------------------------------------
def find_orfs(seq: Seq, min_aa: int = 100, table: int = 1) -> list:
    """
    Find ORFs (ATG ... stop) in all six reading frames.

    Returns a list of dicts sorted by length (longest first). Coordinates are
    1-based, inclusive, and always relative to the forward strand.
    Only the outermost (longest) ORF per stop codon is reported.
    """
    stops = set(CodonTable.unambiguous_dna_by_id[table].stop_codons)
    n = len(seq)
    orfs = []

    for strand, s_seq in (("+", seq), ("-", seq.reverse_complement())):
        s = str(s_seq)
        for frame in range(3):
            i = frame
            while i + 3 <= n:
                if s[i:i + 3] == "ATG":
                    j = i
                    while j + 3 <= n and s[j:j + 3] not in stops:
                        j += 3
                    if j + 3 <= n:  # a stop codon was reached
                        aa_len = (j - i) // 3
                        if aa_len >= min_aa:
                            if strand == "+":
                                start, end = i + 1, j + 3
                            else:
                                start, end = n - (j + 3) + 1, n - i
                            orfs.append({
                                "strand": strand,
                                "frame": frame + 1,
                                "start": start,
                                "end": end,
                                "nt_length": j + 3 - i,
                                "aa_length": aa_len,
                                "protein": str(Seq(s[i:j]).translate(table=table)),
                            })
                        i = j + 3  # continue after the stop; skips nested ATGs
                        continue
                    break  # ran off the end of the sequence without a stop
                i += 3

    return sorted(orfs, key=lambda o: o["aa_length"], reverse=True)


# --------------------------------------------------------------------------
# Translation
# --------------------------------------------------------------------------
def translate_dna(seq: Seq, frame: int = 1, table: int = 1, to_stop: bool = False) -> str:
    """Translate DNA in frame 1-3 (forward) or -1 to -3 (reverse strand)."""
    if frame not in (1, 2, 3, -1, -2, -3):
        raise ValueError("frame must be one of 1, 2, 3, -1, -2, -3")
    if frame < 0:
        seq = seq.reverse_complement()
    seq = seq[abs(frame) - 1:]
    seq = seq[: len(seq) - len(seq) % 3]  # drop trailing partial codon
    return str(seq.translate(table=table, to_stop=to_stop))


# --------------------------------------------------------------------------
# CLI commands
# --------------------------------------------------------------------------
def cmd_download(a):
    path = download_fasta(a.accession, a.email, a.output)
    rec = load_records(path)[0]
    print(f"Saved {rec.id} ({len(rec.seq):,} bp) to {path}")


def cmd_stats(a):
    for rec in load_records(a.fasta):
        st = sequence_stats(rec)
        print(f"\n>{st['id']}  {st['description']}")
        print(f"  Length        : {st['length']:,} bp")
        print(f"  GC content    : {st['gc_percent']:.2f}%")
        print(f"  AT content    : {st['at_percent']:.2f}%")
        print(f"  GC skew       : {st['gc_skew']:+.4f}")
        print(f"  Ambiguous (N) : {st['ambiguous']:,}")
        print(f"  Starts w/ ATG : {st['starts_with_atg']}")
        print("  Composition   : " + ", ".join(
            f"{b}={n:,} ({100 * n / st['length']:.1f}%)" for b, n in st["composition"].items()))


def cmd_orfs(a):
    rec = load_records(a.fasta)[0]
    seq = clean_seq(rec)
    orfs = find_orfs(seq, min_aa=a.min_aa, table=a.table)
    print(f"{rec.id}: {len(orfs)} ORFs >= {a.min_aa} aa (6 frames, ATG start)\n")
    print(f"{'#':>3} {'strand':>6} {'frame':>5} {'start':>8} {'end':>8} {'nt':>7} {'aa':>6}  protein (first 30 aa)")
    for k, o in enumerate(orfs[: a.top], 1):
        print(f"{k:>3} {o['strand']:>6} {o['frame']:>5} {o['start']:>8} {o['end']:>8} "
              f"{o['nt_length']:>7} {o['aa_length']:>6}  {o['protein'][:30]}")
    if len(orfs) > a.top:
        print(f"... {len(orfs) - a.top} more (raise --top to show them)")

    if a.proteins_out:
        out = [
            SeqRecord(Seq(o["protein"]),
                      id=f"{rec.id}_ORF{k}",
                      description=f"strand={o['strand']} frame={o['frame']} {o['start']}-{o['end']} {o['aa_length']}aa")
            for k, o in enumerate(orfs, 1)
        ]
        SeqIO.write(out, a.proteins_out, "fasta")
        print(f"\nWrote {len(out)} protein sequences to {a.proteins_out}")


def cmd_translate(a):
    rec = load_records(a.fasta)[0]
    seq = clean_seq(rec)
    if a.start or a.end:
        start = (a.start or 1) - 1
        end = a.end or len(seq)
        seq = seq[start:end]
        # for region translation the region itself starts at its own frame 1
    protein = translate_dna(seq, frame=a.frame, table=a.table, to_stop=a.to_stop)
    print(f">{rec.id}_translated frame={a.frame} length={len(protein)}aa")
    for i in range(0, len(protein), 60):
        print(protein[i:i + 60])


def cmd_revcomp(a):
    recs = load_records(a.fasta)
    out = [SeqRecord(clean_seq(r).reverse_complement(), id=r.id + "_revcomp", description="reverse complement")
           for r in recs]
    if a.output:
        SeqIO.write(out, a.output, "fasta")
        print(f"Wrote {len(out)} record(s) to {a.output}")
    else:
        SeqIO.write(out, sys.stdout, "fasta")


def build_parser():
    p = argparse.ArgumentParser(description="DNA Sequence Analyzer & FASTA Toolkit")
    sub = p.add_subparsers(dest="command", required=True)

    d = sub.add_parser("download", help="Download a record from NCBI GenBank as FASTA")
    d.add_argument("accession", help="e.g. NC_045512.2")
    d.add_argument("--email", required=True, help="Your email (required by NCBI)")
    d.add_argument("-o", "--output", default=None, help="Output FASTA path")
    d.set_defaults(func=cmd_download)

    s = sub.add_parser("stats", help="Length, GC content, base composition")
    s.add_argument("fasta")
    s.set_defaults(func=cmd_stats)

    o = sub.add_parser("orfs", help="Find ORFs in all six frames")
    o.add_argument("fasta")
    o.add_argument("--min-aa", type=int, default=100, help="Minimum ORF length in amino acids (default 100)")
    o.add_argument("--table", type=int, default=1, help="NCBI genetic code table (default 1)")
    o.add_argument("--top", type=int, default=20, help="Rows to print (default 20)")
    o.add_argument("--proteins-out", help="Write ORF proteins to this FASTA file")
    o.set_defaults(func=cmd_orfs)

    t = sub.add_parser("translate", help="Translate DNA to protein")
    t.add_argument("fasta")
    t.add_argument("--frame", type=int, default=1, choices=[1, 2, 3, -1, -2, -3])
    t.add_argument("--start", type=int, help="Region start (1-based)")
    t.add_argument("--end", type=int, help="Region end (1-based, inclusive)")
    t.add_argument("--table", type=int, default=1)
    t.add_argument("--to-stop", action="store_true", help="Stop at the first stop codon")
    t.set_defaults(func=cmd_translate)

    r = sub.add_parser("revcomp", help="Reverse complement")
    r.add_argument("fasta")
    r.add_argument("-o", "--output")
    r.set_defaults(func=cmd_revcomp)
    return p


def main():
    args = build_parser().parse_args()
    if args.command == "download" and args.output is None:
        args.output = f"{args.accession}.fasta"
    args.func(args)


if __name__ == "__main__":
    main()
