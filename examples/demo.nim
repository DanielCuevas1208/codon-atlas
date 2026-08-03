## Runnable example: build a Codon Atlas report for the bundled sample file.

import std/os
import codonatlas

let samplePath = currentSourcePath().parentDir() / ".." / "tests" / "data" / "sample.fasta"
let records = readFasta(samplePath)

stdout.write buildReport(records, source = "sample.fasta",
  options = ReportOptions(kmer: 3, kmerTop: 5),
  generatedAt = "2026-08-02T12:00:00Z")
