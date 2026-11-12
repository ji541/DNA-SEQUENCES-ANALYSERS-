# DNA Sequence Analyzer & FASTA Toolkit

A small educational bioinformatics project for inspecting DNA/RNA FASTA records. It combines a Python command-line toolkit with an interactive Streamlit interface for sequence summaries, base composition, local GC content, six-frame ORF candidate searches, translation, reverse complements, and optional NCBI sequence retrieval.

> **Scientific note:** Outputs are computational summaries and candidate predictions. An ORF is not, by itself, evidence that a sequence is expressed or encodes a functional protein. Interpret results in the context of the organism, genetic code, annotation, and experimental evidence.

## What it can do

- Read one or more FASTA records and normalize sequence letters to uppercase DNA (RNA `U` is converted to `T`). IUPAC ambiguity symbols are accepted.
- Summarize sequence length, A/C/G/T counts and percentages, GC/AT content, GC skew, and ambiguous-base count.
- Plot base composition and GC content in overlapping sequence windows.
- Search all three frames on both strands for ATG-to-stop ORF candidates, including nested start codons.
- Translate a full sequence or a selected region in any of six reading frames and choose an NCBI genetic code table.
- Export reverse complements, ORF candidate proteins, and translated proteins as FASTA.
- Retrieve nucleotide FASTA records from NCBI using an accession and contact email.

## Project files

- `fasta_toolkit.py` — reusable analysis functions and command-line interface.
- `streamlit_app.py` — interactive sequence exploration app.
- `requirements.txt` — Python package dependencies.

## Setup

Python 3.9 or newer is recommended.

```bash
python -m venv .venv
```

Activate the environment, then install dependencies:

```bash
# macOS/Linux
source .venv/bin/activate

# Windows PowerShell
.venv\Scripts\Activate.ps1

python -m pip install -r requirements.txt
```

## Run the Streamlit app

```bash
streamlit run streamlit_app.py
```

Use the sidebar to upload a FASTA file, paste a FASTA/raw sequence, or fetch a record from NCBI. Select a record and genetic code, then explore the Overview, Composition & GC, ORF candidates, and Translation tabs.

For NCBI retrieval, provide an accession (for example `NC_045512.2`) and a contact email as requested by NCBI E-utilities. Network access is required.

## Command-line examples

### Download a nucleotide FASTA record

```bash
python fasta_toolkit.py download NC_045512.2 --email you@example.org -o reference.fasta
```

### Summarize every record in a FASTA file

```bash
python fasta_toolkit.py stats reference.fasta
```

### Search six frames for ORF candidates

```bash
python fasta_toolkit.py orfs reference.fasta --min-aa 100 --table 1 --proteins-out candidates.faa
```

Useful options include `--min-aa` (minimum amino-acid length), `--table` (NCBI genetic code table), and `--top` (maximum number of candidates printed per record). The optional protein FASTA contains all matching candidates.

### Translate a sequence or region

```bash
# Translate the first record in forward frame 1
python fasta_toolkit.py translate reference.fasta --frame 1

# Coordinates are 1-based and inclusive
python fasta_toolkit.py translate reference.fasta --start 21563 --end 25384 --frame 1 --to-stop

# Select a different FASTA record or reverse frame
python fasta_toolkit.py translate reference.fasta --record 2 --frame -1 --table 11
```

### Generate reverse complements

```bash
python fasta_toolkit.py revcomp reference.fasta -o reverse_complements.fasta
```

Run `python fasta_toolkit.py --help` or append `--help` to any subcommand for the complete options.

## Method notes and interpretation

- GC and AT percentages are calculated among canonical A/C/G/T bases; ambiguous bases are reported separately and excluded from those percentages.
- GC skew is calculated as `(G - C) / (G + C)`; it is descriptive and should not be treated as a stand-alone biological conclusion.
- ORFs use ATG as the start codon and the selected genetic code's stop codons. The terminal stop is included in nucleotide coordinates/length but excluded from protein length.
- ORF coordinates are 1-based and inclusive relative to the original input sequence, including reverse-strand candidates.
- ORF detection is a simple sequence scan. It does not account for introns, alternative start codons, regulatory elements, transcript evidence, or protein function.
- Partial trailing codons are ignored during translation. Ambiguous codons may translate to `X` according to Biopython's translation behavior.
- The chosen genetic code should match the organism and genomic compartment being studied.

## Dependencies

Biopython, Streamlit, NumPy, and Pandas. Install the versions listed in `requirements.txt`.

## License

No license file is currently included. Please add a license before redistributing or reusing this project publicly.
