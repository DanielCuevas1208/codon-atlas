# Changelog

This project follows semantic versioning.

## [0.2.0] - 2026-08-03

Second coherent release. This release completes the FASTQ milestone from the
roadmap.

Added:

- FASTQ parsing with headers, descriptions, wrapped lines, and quality lines.
- Phred+33 quality score statistics per record.
- A FASTQ text report with per-record quality summaries.
- A `--format` option and automatic format detection by file extension.
- A `--quality` option to add a below-threshold count line.
- A bundled FASTQ fixture and a golden report comparison.
- Unit tests for the FASTQ parser and the quality statistics.

The package uses only the Nim standard library. It has no external
dependencies, so it needs no lockfile.

## [0.1.0] - 2026-08-02

First coherent release.

Added:

- FASTA parsing with headers, descriptions, and multi-line sequences.
- GC content and per-base nucleotide counts.
- Overlapping k-mer counting.
- Codon usage tables built from the standard genetic code.
- A text report generator.
- A command-line interface.
- Unit tests for known sequences and a golden report comparison.
- A continuous integration workflow for Linux, macOS, and Windows.
