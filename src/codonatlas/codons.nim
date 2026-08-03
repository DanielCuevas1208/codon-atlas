## Codon usage tables.
##
## Codons are counted in reading frame `start`, which is 0 by default. Codons
## that contain ambiguous bases are skipped. Each counted codon is normalized
## to its RNA form, so the table is always consistent. The returned table
## always lists all 64 codons of the standard genetic code.

import std/algorithm
import std/strutils
import std/tables

import codonatlas/composition
import codonatlas/gencode

type
  CodonEntry* = object
    codon*: string
      ## RNA form of the codon (T replaced with U).
    aminoAcid*: char
      ## One-letter amino acid code. `*` marks a stop codon.
    aminoAcidName*: string
    count*: int
    frequency*: float
      ## Fraction of all counted codons, in the range 0.0 to 1.0.
    perMille*: float
      ## Count per 1000 counted codons.

  CodonUsage* = object
    totalCodons*: int
      ## Number of codons counted, excluding skipped codons.
    entries*: seq[CodonEntry]
      ## All 64 codons, sorted by count descending and then by codon.

proc cmpEntry(a, b: CodonEntry): int =
  if a.count != b.count:
    cmp(b.count, a.count)
  else:
    cmp(a.codon, b.codon)

proc buildEntries(usage: var CodonUsage, counts: CountTable[string]) =
  for codon in sortedCodons:
    let
      n = counts.getOrDefault(codon)
      total = usage.totalCodons
      freq = if total > 0: n.float / total.float else: 0.0
      aa = codonAminoAcid(codon)
    usage.entries.add CodonEntry(
      codon: codon,
      aminoAcid: aa,
      aminoAcidName: if aa == '*': "Stop" else: AminoAcidNames.getOrDefault(aa, "Unknown"),
      count: n,
      frequency: freq,
      perMille: freq * 1000.0)
  usage.entries.sort(cmpEntry)

proc countCodonUsageMany*(sequences: openArray[string], start = 0): CodonUsage =
  ## Count codons across many sequences. Each sequence starts in frame
  ## `start`. Counts are summed across all sequences.
  var counts = initCountTable[string]()
  for sequence in sequences:
    var i = start
    while i + 2 < sequence.len:
      let raw = toUpperAscii(sequence[i ..< i + 3])
      if isCanonical(raw):
        counts.inc(dnaToRna(raw))
        inc result.totalCodons
      i += 3
  buildEntries(result, counts)

proc countCodonUsage*(sequence: string, start = 0): CodonUsage =
  ## Count codons in a single sequence, starting in frame `start`.
  countCodonUsageMany([sequence], start)
