# DNA Sequence Analyzer & FASTA Toolkit

A Python toolkit for analyzing DNA sequences, detecting open reading frames (ORFs), translating genetic code, and exploring sequence composition through a Streamlit web app.

This project includes:
- A reusable FASTA processing toolkit in `fasta_toolkit.py`
- A browser-based interface in `streamlit_app.py`
- A small dependency list in `requirements.txt`

## Features

- Parse and validate FASTA records
- Clean and normalize DNA sequence data
- Compute basic sequence statistics
  - length
  - GC content
  - GC skew
  - AT content
  - ambiguous base counts
  - base composition breakdown
- Detect ORFs across all six reading frames
- Translate DNA sequences into protein sequences
- Reverse complement sequences
- Download results as FASTA files
- Fetch sequences from NCBI GenBank using an accession number
- Explore results visually in a Streamlit app

## Project structure

- `fasta_toolkit.py` — command-line toolkit and core analysis functions
- `streamlit_app.py` — interactive web dashboard for sequence analysis
- `requirements.txt` — Python dependencies
- `README.md` — project documentation

## Requirements

Python 3.9+ is recommended.

Install dependencies:

```bash
pip install -r requirements.txt
```

## Quick start

### 1. Run the web app

```bash
streamlit run streamlit_app.py
```

The app lets you:
- upload a FASTA file
- paste a DNA sequence directly
- fetch a record from NCBI GenBank
- inspect GC content and sequence composition
- view six-frame ORF hits
- translate DNA regions in different frames

### 2. Use the CLI toolkit

#### Download a sequence from NCBI

```bash
python fasta_toolkit.py download NC_045512.2 --email you@example.com -o NC_045512.2.fasta
```

#### Show sequence statistics

```bash
python fasta_toolkit.py stats NC_045512.2.fasta
```

#### Find ORFs

```bash
python fasta_toolkit.py orfs NC_045512.2.fasta --min-aa 100 --proteins-out orfs.faa
```

#### Translate a sequence

```bash
python fasta_toolkit.py translate NC_045512.2.fasta --frame 1
```

Translate a specific region:

```bash
python fasta_toolkit.py translate NC_045512.2.fasta --start 21563 --end 25384 --to-stop
```

#### Reverse complement

```bash
python fasta_toolkit.py revcomp NC_045512.2.fasta -o revcomp.fasta
```

## Example use cases

- Analyze viral genomes or reference sequences
- Identify candidate open reading frames in a gene region
- Translate coding segments into amino acid sequences
- Compare sequence composition across records
- Explore a FASTA dataset interactively in the browser

## Notes

- ORF detection is based on ATG start codons and stops defined by the selected genetic code table.
- NCBI downloads require a valid email address.
- The app uses Biopython and Streamlit for sequence parsing, translation, and visualization.

## Dependencies

Key libraries used by this project:
- Biopython
- Streamlit
- NumPy
- Pandas

## License

This project does not currently include a repository license file. If you plan to distribute it publicly, consider adding an appropriate open-source license.
