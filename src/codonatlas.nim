## Codon Atlas: composition statistics for DNA and RNA sequence files.
##
## The package reads FASTA and FASTQ files, computes GC content and k-mer
## counts, and builds codon usage tables. It summarizes FASTQ quality scores.
## A CLI writes the results to a plain-text report.

import codonatlas/[fasta, fastq, quality, composition, gencode, codons, report]

export fasta, fastq, quality, composition, gencode, codons, report
