# Codon Atlas

Codon Atlas is a Nim library for DNA and RNA sequence files. It computes
composition statistics. These include GC content, k-mer counts, and codon
usage tables. A command-line tool writes the results to a text report.

## Features

- Read FASTA files with a permissive parser.
- Compute GC content and per-base counts.
- Count overlapping k-mers.
- Build codon usage tables from the standard genetic code.
- Write a text report to the terminal or a file.
- Share one API across the library and the CLI.

## Requirements

- Nim 2.0 or newer.
- Nimble 0.14 or newer.
- A C compiler, such as gcc or MSVC. Nim compiles to C.

## Setup

Build the CLI:

```
nim c -d:release --path:src -o:build/codonatlas src/codonatlas_cli.nim
```

Run the test suite:

```
nim c -r --path:src --outdir:build/tests tests/test_codonatlas.nim
```

On a checked-out repository, `nimble build` and `nimble test` are the
equivalent shorthand.

## Usage

### Command line

Produce a report for a FASTA file:

```
./build/codonatlas tests/data/sample.fasta
```

Write the report to a file:

```
./build/codonatlas tests/data/sample.fasta --out=report.txt
```

Add a k-mer count section:

```
./build/codonatlas tests/data/sample.fasta --kmer=3
```

Show only the top 10 k-mers:

```
./build/codonatlas tests/data/sample.fasta --kmer=3 --top=10
```

### Library

Import the package and build a report in code:

```nim
import codonatlas

let records = readFasta("tests/data/sample.fasta")
echo buildReport(records, source = "sample.fasta")
```

Run the bundled example:

```
nim c -r --path:src --outdir:build/examples examples/demo.nim
```

## Sample output

The command below prints a report for the bundled sample file:

```
./build/codonatlas tests/data/sample.fasta --kmer=3 --top=3
```

```
Codon Atlas report
==================
Source: tests/data/sample.fasta
Generated: 2026-08-03T00:27:23Z
Records: 2
Total canonical bases: 181

Nucleotide composition
----------------------
Record 1: demo_orf
  Description: synthetic coding sequence, 144 bp, starts with ATG and ends with TAA
  Length: 144 canonical bases
  Bases: A 40  C 31  G 38  T 35  U 0  N 0
  GC content: 47.9%

Record 2: promoter_fragment
  Description: synthetic non-coding fragment
  Length: 37 canonical bases
  Bases: A 9  C 8  G 9  T 11  U 0  N 0
  GC content: 45.9%

Overall
-------
Total canonical bases: 181
GC content: 47.5%

K-mer counts (k = 3)
--------------------------------
Top 3 of 61 distinct k-mers.
CGT        8
GTA        7
GCG        6

Codon usage table (reading frame 0)
-----------------------------------
Total codons: 60

Count   Frequency   Per 1000  Codon AA  Amino acid
--------------------------------------------
    4      0.0667       66.7  CGU  R   Arginine
    3      0.0500       50.0  AUG  M   Methionine
    3      0.0500       50.0  GAC  D   Aspartic acid
    3      0.0500       50.0  GAU  D   Aspartic acid
    2      0.0333       33.3  AAG  K   Lysine
```

The report lists all 64 codons of the standard genetic code. The output above
is abbreviated. Codons use RNA letters (U instead of T). A stop codon shows as
`*`.

## Project layout

```
codonatlas.nimble          package metadata and tasks
src/codonatlas.nim         public API entry point
src/codonatlas/fasta.nim   FASTA parsing
src/codonatlas/composition.nim   GC content, base counts, k-mers
src/codonatlas/gencode.nim standard genetic code
src/codonatlas/codons.nim  codon usage tables
src/codonatlas/report.nim  text report builder
src/codonatlas_cli.nim     command-line interface
examples/demo.nim          runnable example
tests/                     unit tests and fixtures
```

## Testing

Run `nimble check` to typecheck every module. Run `nimble test` to execute
the unit tests. The suite covers parsing, statistics, the genetic code, the
report builder, and a golden report comparison.

## Limitations

- The parser keeps sequences as given. It does not reverse-complement.
- Codon counts use reading frame 0 on the given strands.
- Ambiguous bases (N and others) are excluded from composition ratios.
- The table normalizes codons to RNA form. It merges DNA and RNA input.

## Roadmap

Release 2 adds FASTQ support. It will summarize quality scores per record.

Release 3 adds translation. It will map codons to amino acids. It will also
detect open reading frames.

Release 4 adds reverse complement and configurable reading frames.

## License

Released under the MIT License. See `LICENSE`.
