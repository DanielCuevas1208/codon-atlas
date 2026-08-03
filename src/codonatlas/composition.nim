## Composition statistics for nucleotide sequences.
##
## All statistics treat bases case-insensitively. A base is canonical when it
## is A, C, G, T, or U. Every other base is ambiguous and is excluded from
## composition ratios.

import std/strutils
import std/tables

type
  NucleotideCounts* = object
    a*: int
    c*: int
    g*: int
    t*: int
    u*: int
    n*: int
      ## Ambiguous bases (any base outside A, C, G, T, U).

proc nucleotideCounts*(sequence: string): NucleotideCounts =
  ## Count the bases in a sequence.
  for base in sequence:
    case base
    of 'A', 'a': result.a += 1
    of 'C', 'c': result.c += 1
    of 'G', 'g': result.g += 1
    of 'T', 't': result.t += 1
    of 'U', 'u': result.u += 1
    else: result.n += 1

proc totalBases*(counts: NucleotideCounts): int =
  ## Number of canonical bases in `counts`. Ambiguous bases are not counted.
  counts.a + counts.c + counts.g + counts.t + counts.u

proc isCanonical*(sequence: string): bool =
  ## Return true when every base in `sequence` is canonical.
  for base in sequence:
    case base
    of 'A', 'a', 'C', 'c', 'G', 'g', 'T', 't', 'U', 'u':
      discard
    else:
      return false
  true

proc gcContent*(sequence: string): float =
  ## Fraction of canonical bases that are G or C, in the range 0.0 to 1.0.
  ##
  ## Returns 0.0 for a sequence with no canonical bases.
  let counts = nucleotideCounts(sequence)
  let total = totalBases(counts)
  if total == 0:
    return 0.0
  (counts.g + counts.c).float / total.float

proc gcPercent*(sequence: string): float =
  ## GC content as a percentage, in the range 0.0 to 100.0.
  gcContent(sequence) * 100.0

proc kmerCounts*(sequence: string, k: int): CountTable[string] =
  ## Count every overlapping window of length `k` that contains only
  ## canonical bases.
  ##
  ## Windows with ambiguous bases are skipped. A window is reported with
  ## uppercase bases. Returns an empty table when `k` is not positive or when
  ## the sequence is shorter than `k`.
  if k <= 0:
    return
  for i in 0 .. sequence.len - k:
    let window = sequence[i ..< i + k]
    if not isCanonical(window):
      continue
    result.inc(toUpperAscii(window))
