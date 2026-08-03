## Plain-text report generation.

import std/algorithm
import std/sequtils
import std/strutils
import std/tables
import std/times

import codonatlas/codons
import codonatlas/composition
import codonatlas/fasta
import codonatlas/fastq
import codonatlas/quality

type
  ReportOptions* = object
    kmer*: int
      ## When greater than 0, add a k-mer count section for this k.
    kmerTop*: int
      ## When greater than 0, show only the top k-mers in that section.
    qualityThreshold*: int
      ## When greater than 0, add a below-threshold count for each FASTQ
      ## record's quality line.

proc cmpKmer(a, b: tuple[key: string, val: int]): int =
  if a.val != b.val:
    cmp(b.val, a.val)
  else:
    cmp(a.key, b.key)

proc fmtDecimal(value: float, precision: int): string =
  formatFloat(value, ffDecimal, precision)

proc formatCounts(counts: NucleotideCounts): string =
  "A " & $counts.a & "  C " & $counts.c & "  G " & $counts.g &
    "  T " & $counts.t & "  U " & $counts.u & "  N " & $counts.n

proc kmerSection(sequences: openArray[string], k, top: int): seq[string] =
  var combined = initCountTable[string]()
  for sequence in sequences:
    for key, value in kmerCounts(sequence, k):
      combined.inc(key, value)
  var items = combined.pairs.toSeq
  items.sort(cmpKmer)
  let shown = if top > 0: min(top, items.len) else: items.len
  result.add ""
  result.add "K-mer counts (k = " & $k & ")"
  result.add "--------------------------------"
  if top > 0:
    result.add "Top " & $shown & " of " & $items.len & " distinct k-mers."
  for i in 0 ..< shown:
    result.add alignLeft(items[i][0], 4) & "  " & align($items[i][1], 6)

proc gcPercentOf(counts: NucleotideCounts): float =
  let total = totalBases(counts)
  if total == 0:
    0.0
  else:
    (counts.g + counts.c).float * 100.0 / total.float

proc buildReport*(records: seq[Record], source: string,
                  options = ReportOptions(), generatedAt = ""): string =
  ## Build a text report for a list of FASTA records.
  ##
  ## `source` names the origin of the data, such as a file path. When
  ## `generatedAt` is empty, the current UTC time is used.
  let whenGenerated =
    if generatedAt.len > 0: generatedAt
    else: now().utc.format("yyyy-MM-dd'T'HH:mm:ss'Z'")

  var
    outLines: seq[string]
    allCounts = NucleotideCounts()

  for record in records:
    let counts = nucleotideCounts(record.sequence)
    allCounts.a += counts.a
    allCounts.c += counts.c
    allCounts.g += counts.g
    allCounts.t += counts.t
    allCounts.u += counts.u
    allCounts.n += counts.n

  outLines.add "Codon Atlas report"
  outLines.add "=================="
  outLines.add "Source: " & source
  outLines.add "Generated: " & whenGenerated
  outLines.add "Records: " & $records.len
  outLines.add "Total canonical bases: " & $totalBases(allCounts)
  outLines.add ""

  outLines.add "Nucleotide composition"
  outLines.add "----------------------"
  for index, record in records:
    let counts = nucleotideCounts(record.sequence)
    outLines.add "Record " & $(index + 1) & ": " & record.id
    if record.description.len > 0:
      outLines.add "  Description: " & record.description
    outLines.add "  Length: " & $totalBases(counts) & " canonical bases"
    outLines.add "  Bases: " & formatCounts(counts)
    outLines.add "  GC content: " & fmtDecimal(gcPercent(record.sequence), 1) & "%"
    outLines.add ""

  outLines.add "Overall"
  outLines.add "-------"
  outLines.add "Total canonical bases: " & $totalBases(allCounts)
  outLines.add "GC content: " & fmtDecimal(gcPercentOf(allCounts), 1) & "%"

  if options.kmer > 0:
    outLines.add kmerSection(records.mapIt(it.sequence), options.kmer, options.kmerTop)

  let usage = countCodonUsageMany(records.mapIt(it.sequence))
  outLines.add ""
  outLines.add "Codon usage table (reading frame 0)"
  outLines.add "-----------------------------------"
  outLines.add "Total codons: " & $usage.totalCodons
  outLines.add ""
  var countWidth = 5
  for entry in usage.entries:
    countWidth = max(countWidth, len($entry.count))
  outLines.add align("Count", countWidth) & "  " &
    align("Frequency", 10) & "  " & align("Per 1000", 9) & "  " &
    alignLeft("Codon", 4) & " " & alignLeft("AA", 2) & "  Amino acid"
  outLines.add repeat("-", 44)
  for entry in usage.entries:
    outLines.add align($entry.count, countWidth) & "  " &
      align(fmtDecimal(entry.frequency, 4), 10) & "  " &
      align(fmtDecimal(entry.perMille, 1), 9) & "  " &
      alignLeft(entry.codon, 4) & " " & alignLeft($entry.aminoAcid, 2) & "  " &
      entry.aminoAcidName
  outLines.add ""
  outLines.add "Codons use RNA letters (U instead of T). AA uses one-letter"
  outLines.add "amino acid codes. A '*' marks a stop codon."

  outLines.join("\n") & "\n"

proc buildFastqReport*(records: seq[FastqRecord], source: string,
                       options = ReportOptions(), generatedAt = ""): string =
  ## Build a text report for a list of FASTQ records.
  ##
  ## `source` names the origin of the data, such as a file path. When
  ## `generatedAt` is empty, the current UTC time is used.
  let whenGenerated =
    if generatedAt.len > 0: generatedAt
    else: now().utc.format("yyyy-MM-dd'T'HH:mm:ss'Z'")

  var
    outLines: seq[string]
    allCounts = NucleotideCounts()

  for record in records:
    let counts = nucleotideCounts(record.sequence)
    allCounts.a += counts.a
    allCounts.c += counts.c
    allCounts.g += counts.g
    allCounts.t += counts.t
    allCounts.u += counts.u
    allCounts.n += counts.n

  outLines.add "Codon Atlas report"
  outLines.add "=================="
  outLines.add "Source: " & source
  outLines.add "Generated: " & whenGenerated
  outLines.add "Records: " & $records.len
  outLines.add "Total canonical bases: " & $totalBases(allCounts)
  outLines.add ""

  outLines.add "Nucleotide composition"
  outLines.add "----------------------"
  for index, record in records:
    let counts = nucleotideCounts(record.sequence)
    outLines.add "Record " & $(index + 1) & ": " & record.id
    if record.description.len > 0:
      outLines.add "  Description: " & record.description
    outLines.add "  Length: " & $totalBases(counts) & " canonical bases"
    outLines.add "  Bases: " & formatCounts(counts)
    outLines.add "  GC content: " & fmtDecimal(gcPercent(record.sequence), 1) & "%"
    let summary = qualitySummary(record.quality)
    outLines.add "  Quality: mean " & fmtDecimal(summary.mean, 1) &
      ", min " & $summary.min & ", max " & $summary.max
    if options.qualityThreshold > 0:
      let below = basesBelowQuality(record.quality, options.qualityThreshold)
      outLines.add "  Bases below Q" & $options.qualityThreshold & ": " & $below
    outLines.add ""

  outLines.add "Overall"
  outLines.add "-------"
  outLines.add "Total canonical bases: " & $totalBases(allCounts)
  outLines.add "GC content: " & fmtDecimal(gcPercentOf(allCounts), 1) & "%"

  if options.kmer > 0:
    outLines.add kmerSection(records.mapIt(it.sequence), options.kmer, options.kmerTop)

  outLines.add ""
  outLines.add "Quality scores use the Phred+33 encoding. A higher score"
  outLines.add "means a more confident base call."

  outLines.join("\n") & "\n"
