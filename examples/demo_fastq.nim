## Runnable example: build a Codon Atlas report for the bundled FASTQ file.

import std/os
import codonatlas

let samplePath = currentSourcePath().parentDir() / ".." / "tests" / "data" / "sample.fastq"
let records = readFastq(samplePath)

stdout.write buildFastqReport(records, source = "sample.fastq",
  options = ReportOptions(kmer: 3, qualityThreshold: 20),
  generatedAt = "2026-08-02T12:00:00Z")
