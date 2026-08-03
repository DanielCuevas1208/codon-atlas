# Codon Atlas

Codon Atlas is a Nim library for DNA and RNA sequence files. It reads FASTA
and FASTQ files. It computes GC content, k-mer counts, and codon usage tables.
It summarizes FASTQ quality scores. A command-line tool writes the results to
a text report.

## Features

- Read FASTA files with a permissive parser.
- Read FASTQ files with their quality lines.
- Compute GC content and per-base counts.
- Count overlapping k-mers.
- Build codon usage tables from the standard genetic code.
- Summarize FASTQ quality scores per record.
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
nim c -r --path:src --outdir:build/tests tests/test_fastq.nim
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

Produce a report for a FASTQ file. The tool detects the format from the file
extension:

```
./build/codonatlas tests/data/sample.fastq
```

Add a quality threshold line to a FASTQ report:

```
./build/codonatlas tests/data/sample.fastq --quality=20
```

Force a format when the file extension is not standard:

```
./build/codonatlas reads.data --format=fastq
```

### Library

Import the package and build a report in code:

```nim
import codonatlas

let records = readFasta("tests/data/sample.fasta")
echo buildReport(records, source = "sample.fasta")
```

Build a report from a FASTQ file:

```nim
import codonatlas

let records = readFastq("tests/data/sample.fastq")
echo buildFastqReport(records, source = "sample.fastq",
  options = ReportOptions(qualityThreshold = 20))
```

Run the bundled examples:

```
nim c -r --path:src --outdir:build/examples examples/demo.nim
nim c -r --path:src --outdir:build/examples examples/demo_fastq.nim
```

## Sample output

The command below prints a FASTQ report for the bundled sample file:

```
./build/codonatlas tests/data/sample.fastq --quality=20
```

```
Codon Atlas report
==================
Source: tests/data/sample.fastq
Generated: 2026-08-03T02:15:15Z
Records: 2
Total canonical bases: 35

Nucleotide composition
----------------------
Record 1: read1
  Description: synthetic high-quality read, 23 bp
  Length: 23 canonical bases
  Bases: A 5  C 5  G 7  T 6  U 0  N 0
  GC content: 52.2%
  Quality: mean 38.3, min 37, max 40
  Bases below Q20: 0

Record 2: read2
  Description: synthetic low-quality read, 12 bp
  Length: 12 canonical bases
  Bases: A 3  C 4  G 2  T 3  U 0  N 0
  GC content: 50.0%
  Quality: mean 7.0, min 0, max 13
  Bases below Q20: 12

Overall
-------
Total canonical bases: 35
GC content: 51.4%
```

Quality scores use the Phred+33 encoding. A higher score means a more
confident base call. A score of 20 is a common quality threshold.

The FASTA report adds a codon usage table. It lists all 64 codons of the
standard genetic code. Codons use RNA letters (U instead of T). A stop codon
shows as `*`.

## Project layout

```
codonatlas.nimble          package metadata and tasks
src/codonatlas.nim         public API entry point
src/codonatlas/fasta.nim   FASTA parsing
src/codonatlas/fastq.nim   FASTQ parsing
src/codonatlas/quality.nim quality score statistics
src/codonatlas/composition.nim   GC content, base counts, k-mers
src/codonatlas/gencode.nim standard genetic code
src/codonatlas/codons.nim  codon usage tables
src/codonatlas/report.nim  text report builders
src/codonatlas_cli.nim     command-line interface
examples/demo.nim          runnable FASTA example
examples/demo_fastq.nim    runnable FASTQ example
tests/                     unit tests and fixtures
```

## Testing

Run `nimble check` to typecheck every module. Run `nimble test` to execute
the unit tests. The suite covers parsing, statistics, the genetic code, the
report builders, quality scores, and golden report comparisons.

## Limitations

- The parser keeps sequences as given. It does not reverse-complement.
- Codon counts use reading frame 0 on the given strands.
- Ambiguous bases (N and others) are excluded from composition ratios.
- The table normalizes codons to RNA form. It merges DNA and RNA input.
- The FASTQ parser assumes Sanger quality scores. It rejects sequences and
  quality lines with different lengths.
- FASTQ reports do not include codon usage tables. Reads are usually short
  unassembled fragments.

## Roadmap

Release 2 is complete. Codon Atlas reads FASTQ files and summarizes quality
scores per record.

Release 3 adds translation. It will map codons to amino acids. It will also
detect open reading frames.

Release 4 adds reverse complement and configurable reading frames.

## License

Released under the MIT License. See `LICENSE`.
