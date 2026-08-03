## Codon Atlas: composition statistics for DNA and RNA sequence files.
##
## The package reads FASTA files, computes GC content and k-mer counts, and
## builds codon usage tables. A CLI writes the results to a plain-text report.

import codonatlas/[fasta, composition, gencode, codons, report]

export fasta, composition, gencode, codons, report
