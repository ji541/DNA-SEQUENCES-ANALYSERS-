#!/usr/bin/env python3
"""FASTA sequence analysis utilities and command-line interface.

This toolkit provides descriptive sequence statistics, six-frame ORF candidate
search, translation, reverse-complement generation, and NCBI nucleotide fetches.
Coordinates reported by the ORF command are 1-based and inclusive on the input
sequence, including for candidates found on the reverse-complement strand.
"""

from __future__ import annotations

import argparse
import sys
from collections import Counter
from pathlib import Path
from typing import Any

from Bio import Entrez, SeqIO
from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

DNA_IUPAC = set("ACGTRYSWKMBDHVN")


def clean_seq(record: SeqRecord) -> Seq:
    """Normalize a record to uppercase DNA, converting RNA U to T."""
    sequence = str(record.seq).upper().replace("U", "T").replace(" ", "").replace("\n", "")
    invalid = sorted(set(sequence) - DNA_IUPAC)
    if invalid:
        raise ValueError(f"Sequence contains unsupported character(s): {', '.join(invalid)}")
    return Seq(sequence)


def load_records(path: str | Path) -> list[SeqRecord]:
    """Read all records from a FASTA file or raise a helpful error."""
    try:
        records = list(SeqIO.parse(str(path), "fasta"))
    except (OSError, UnicodeError) as exc:
        raise ValueError(f"Could not read FASTA file '{path}': {exc}") from exc
    if not records:
        raise ValueError(f"No FASTA records found in '{path}'. Check that it is valid FASTA.")
    for record in records:
        clean_seq(record)  # validate all records before downstream analysis
    return records


def gc_content(seq: Seq | str) -> float:
    """Return GC percentage among unambiguous A/C/G/T bases (0-100)."""
    sequence = str(seq).upper().replace("U", "T")
    canonical = sum(sequence.count(base) for base in "ACGT")
    return 100.0 * (sequence.count("G") + sequence.count("C")) / canonical if canonical else 0.0


def sequence_stats(record: SeqRecord) -> dict[str, Any]:
    """Calculate length, base composition, GC/AT content, and GC skew."""
    seq = clean_seq(record)
    counts = Counter(str(seq))
    length = len(seq)
    canonical = sum(counts.get(base, 0) for base in "ACGT")
    g, c = counts.get("G", 0), counts.get("C", 0)
    composition = {base: counts.get(base, 0) for base in "ACGT"}
    percentages = {
        base: (100.0 * count / canonical if canonical else 0.0)
        for base, count in composition.items()
    }
    return {
        "id": record.id,
        "description": record.description,
        "length": length,
        "gc_percent": gc_content(seq),
        "at_percent": (100.0 * (counts.get("A", 0) + counts.get("T", 0)) / canonical
                       if canonical else 0.0),
        "gc_skew": (g - c) / (g + c) if g + c else 0.0,
        "composition": composition,
        "composition_percent": percentages,
        "ambiguous": length - canonical,
        "starts_with_atg": str(seq).startswith("ATG"),
    }


def find_orfs(seq: Seq | str, min_aa: int = 100, table: int = 1) -> list[dict[str, Any]]:
    """Find ATG-to-stop ORF candidates in all six frames.

    Every qualifying ATG is considered, including nested starts. Protein length
    excludes the terminal stop codon; nucleotide length and end coordinate
    include it. Coordinates are 1-based, inclusive, on the original sequence.
    """
    if min_aa < 1:
        raise ValueError("min_aa must be at least 1")
    if table not in CodonTable.unambiguous_dna_by_id:
        raise ValueError(f"Unknown NCBI genetic code table: {table}")
    sequence = Seq(str(seq).upper().replace("U", "T"))
    n = len(sequence)
    stops = set(CodonTable.unambiguous_dna_by_id[table].stop_codons)
    candidates: list[dict[str, Any]] = []

    for strand, oriented in (("+", sequence), ("-", sequence.reverse_complement())):
        text = str(oriented)
        for frame in range(3):
            for start0 in range(frame, n - 2, 3):
                if text[start0:start0 + 3] != "ATG":
                    continue
                stop0 = None
                for pos in range(start0 + 3, n - 2, 3):
                    if text[pos:pos + 3] in stops:
                        stop0 = pos
                        break
                if stop0 is None:
                    continue
                aa_length = (stop0 - start0) // 3
                if aa_length < min_aa:
                    continue
                if strand == "+":
                    start, end = start0 + 1, stop0 + 3
                else:
                    start, end = n - (stop0 + 3) + 1, n - start0
                candidates.append({
                    "strand": strand,
                    "frame": frame + 1,
                    "start": start,
                    "end": end,
                    "nt_length": stop0 + 3 - start0,
                    "aa_length": aa_length,
                    "start_codon": "ATG",
                    "stop_codon": text[stop0:stop0 + 3],
                    "protein": str(Seq(text[start0:stop0]).translate(table=table)),
                })
    return sorted(candidates, key=lambda item: (-item["aa_length"], item["start"], item["strand"]))


def translate_dna(seq: Seq | str, frame: int = 1, table: int = 1, to_stop: bool = False) -> str:
    """Translate DNA in frame 1-3 or reverse-complement frame -1 to -3."""
    if frame not in (1, 2, 3, -1, -2, -3):
        raise ValueError("frame must be one of 1, 2, 3, -1, 2, or -3")
    if table not in CodonTable.unambiguous_dna_by_id:
        raise ValueError(f"Unknown NCBI genetic code table: {table}")
    dna = Seq(str(seq).upper().replace("U", "T"))
    if frame < 0:
        dna = dna.reverse_complement()
    dna = dna[abs(frame) - 1:]
    dna = dna[:len(dna) - len(dna) % 3]
    return str(dna.translate(table=table, to_stop=to_stop))


def download_fasta(accession: str, email: str, out_path: str) -> str:
    """Fetch one nucleotide record from NCBI and save it as FASTA."""
    if not email.strip():
        raise ValueError("An email address is required by NCBI E-utilities.")
    Entrez.email = email.strip()
    try:
        with Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text") as handle:
            data = handle.read()
    except Exception as exc:
        raise RuntimeError(f"NCBI download failed: {exc}") from exc
    if not data.lstrip().startswith(">"):
        raise RuntimeError(f"NCBI returned no FASTA record for accession '{accession}'.")
    Path(out_path).write_text(data, encoding="utf-8")
    return out_path


def _write_fasta(records: list[SeqRecord], path: str | None) -> None:
    if path:
        SeqIO.write(records, path, "fasta")
        print(f"Wrote {len(records)} record(s) to {path}")
    else:
        SeqIO.write(records, sys.stdout, "fasta")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="DNA Sequence Analyzer & FASTA Toolkit")
    sub = parser.add_subparsers(dest="command", required=True)

    download = sub.add_parser("download", help="Download a nucleotide FASTA record from NCBI")
    download.add_argument("accession", help="NCBI accession, for example NC_045512.2")
    download.add_argument("--email", required=True, help="Contact email required by NCBI")
    download.add_argument("-o", "--output", help="Output path (default: <accession>.fasta)")

    stats = sub.add_parser("stats", help="Summarize every FASTA record")
    stats.add_argument("fasta")

    orfs = sub.add_parser("orfs", help="Find ATG-to-stop candidates in six reading frames")
    orfs.add_argument("fasta")
    orfs.add_argument("--min-aa", type=int, default=100, help="Minimum protein length (default: 100)")
    orfs.add_argument("--table", type=int, default=1, help="NCBI genetic code table (default: 1)")
    orfs.add_argument("--top", type=int, default=20, help="Maximum rows to print (default: 20)")
    orfs.add_argument("--proteins-out", help="Optional output FASTA for predicted proteins")

    translate = sub.add_parser("translate", help="Translate a sequence or 1-based inclusive region")
    translate.add_argument("fasta")
    translate.add_argument("--record", type=int, default=1, help="1-based FASTA record index")
    translate.add_argument("--frame", type=int, default=1, choices=[1, 2, 3, -1, -2, -3])
    translate.add_argument("--start", type=int, help="Region start, 1-based")
    translate.add_argument("--end", type=int, help="Region end, 1-based inclusive")
    translate.add_argument("--table", type=int, default=1)
    translate.add_argument("--to-stop", action="store_true", help="Stop at the first in-frame stop codon")

    revcomp = sub.add_parser("revcomp", help="Write reverse complements for all FASTA records")
    revcomp.add_argument("fasta")
    revcomp.add_argument("-o", "--output")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        if args.command == "download":
            output = args.output or f"{args.accession}.fasta"
            path = download_fasta(args.accession, args.email, output)
            record = load_records(path)[0]
            print(f"Saved {record.id} ({len(record.seq):,} bp) to {path}")
        elif args.command == "stats":
            for record in load_records(args.fasta):
                result = sequence_stats(record)
                print(f"\n>{result['id']}  {result['description']}")
                print(f"  Length             : {result['length']:,} bp")
                print(f"  GC content         : {result['gc_percent']:.2f}% (of A/C/G/T)")
                print(f"  AT content         : {result['at_percent']:.2f}% (of A/C/G/T)")
                print(f"  GC skew            : {result['gc_skew']:+.4f}")
                print(f"  Ambiguous bases    : {result['ambiguous']:,}")
                print(f"  Starts with ATG    : {result['starts_with_atg']}")
                print("  Composition        : " + ", ".join(
                    f"{base}={count:,} ({result['composition_percent'][base]:.1f}%)"
                    for base, count in result["composition"].items()))
        elif args.command == "orfs":
            if args.top < 1:
                raise ValueError("--top must be at least 1")
            proteins: list[SeqRecord] = []
            for record in load_records(args.fasta):
                hits = find_orfs(clean_seq(record), min_aa=args.min_aa, table=args.table)
                print(f"\n{record.id}: {len(hits)} ORF candidate(s) >= {args.min_aa} aa")
                print(f"{'#':>3} {'strand':>6} {'frame':>5} {'start':>8} {'end':>8} {'nt':>7} {'aa':>6}  protein (first 30 aa)")
                for index, hit in enumerate(hits[:args.top], 1):
                    print(f"{index:>3} {hit['strand']:>6} {hit['frame']:>5} {hit['start']:>8} {hit['end']:>8} "
                          f"{hit['nt_length']:>7} {hit['aa_length']:>6}  {hit['protein'][:30]}")
                proteins.extend(SeqRecord(Seq(hit["protein"]), id=f"{record.id}_ORF{index}",
                    description=f"strand={hit['strand']} frame={hit['frame']} {hit['start']}-{hit['end']} {hit['aa_length']}aa")
                    for index, hit in enumerate(hits, 1))
            if args.proteins_out:
                _write_fasta(proteins, args.proteins_out)
        elif args.command == "translate":
            records = load_records(args.fasta)
            if not 1 <= args.record <= len(records):
                raise ValueError(f"--record must be between 1 and {len(records)}")
            record = records[args.record - 1]
            sequence = clean_seq(record)
            start = args.start if args.start is not None else 1
            end = args.end if args.end is not None else len(sequence)
            if start < 1 or end < start or end > len(sequence):
                raise ValueError(f"Region must satisfy 1 <= start <= end <= {len(sequence)}")
            protein = translate_dna(sequence[start - 1:end], frame=args.frame, table=args.table, to_stop=args.to_stop)
            print(f">{record.id}_translated region={start}-{end} frame={args.frame} length={len(protein)}aa")
            for offset in range(0, len(protein), 60):
                print(protein[offset:offset + 60])
        elif args.command == "revcomp":
            records = load_records(args.fasta)
            output = [SeqRecord(clean_seq(record).reverse_complement(), id=f"{record.id}_revcomp",
                                description="reverse complement") for record in records]
            _write_fasta(output, args.output)
    except (ValueError, RuntimeError, OSError) as exc:
        parser.error(str(exc))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
