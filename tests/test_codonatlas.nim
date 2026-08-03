import std/os
import std/sequtils
import std/strutils
import std/tables
import std/unittest

import codonatlas

proc fixturePath(): string =
  currentSourcePath().parentDir() / "data" / "sample.fasta"

proc goldenPath(): string =
  currentSourcePath().parentDir() / "data" / "expected_report.txt"

proc entry(usage: CodonUsage, codon: string): CodonEntry =
  for item in usage.entries:
    if item.codon == codon:
      return item
  raise newException(ValueError, "codon not in table: " & codon)

suite "GC content":
  test "hand-checked fractions":
    check gcContent("ATGC") == 0.5
    check gcContent("AAAA") == 0.0
    check gcContent("CCGG") == 1.0
    check gcContent("") == 0.0
    check gcContent("atgcc") == 0.6
    check gcContent("AUGUU") == 0.2
    check gcContent("NNNAT") == 0.0
    check gcContent("NNCGG") == 1.0

  test "percentage":
    check gcPercent("ATGC") == 50.0
    check gcPercent("CCCC") == 100.0

suite "nucleotide counts":
  test "canonical and ambiguous bases":
    let counts = nucleotideCounts("AACCGGTTUUNN")
    check counts.a == 2
    check counts.c == 2
    check counts.g == 2
    check counts.t == 2
    check counts.u == 2
    check counts.n == 2
    check totalBases(counts) == 10

  test "case insensitivity":
    check nucleotideCounts("ATgc") == NucleotideCounts(a: 1, c: 1, g: 1, t: 1)

suite "k-mer counts":
  test "overlapping windows":
    let four = kmerCounts("AAAA", 2)
    check four.len == 1
    check four["AA"] == 3
    let eight = kmerCounts("ATGCATGC", 3)
    check eight.len == 4
    check eight["ATG"] == 2
    check eight["TGC"] == 2
    check eight["GCA"] == 1
    check eight["CAT"] == 1

  test "windows with ambiguous bases are skipped":
    let counts = kmerCounts("AATGNN", 2)
    check counts.len == 3
    check counts["AA"] == 1
    check counts["AT"] == 1
    check counts["TG"] == 1

  test "case insensitivity":
    let counts = kmerCounts("atgca", 3)
    check counts.len == 3
    check counts["ATG"] == 1
    check counts["TGC"] == 1
    check counts["GCA"] == 1

  test "edge cases":
    check kmerCounts("", 3).len == 0
    check kmerCounts("AT", 3).len == 0
    check kmerCounts("ACGT", 0).len == 0
    check kmerCounts("ACGT", 5).len == 0

suite "genetic code":
  test "known codons":
    check aminoAcid("AUG") == 'M'
    check aminoAcid("UUU") == 'F'
    check aminoAcid("UGG") == 'W'
    check aminoAcid("UAA") == '*'
    check aminoAcid("XYZ") == '?'
    check codonAminoAcid("ATG") == 'M'
    check codonAminoAcid("TAA") == '*'

  test "full table is complete":
    check StandardCode.len == 64
    var byAminoAcid = initCountTable[char]()
    for aa in StandardCode.values:
      byAminoAcid.inc(aa)
    check byAminoAcid['M'] == 1
    check byAminoAcid['W'] == 1
    check byAminoAcid['I'] == 3
    check byAminoAcid['L'] == 6
    check byAminoAcid['S'] == 6
    check byAminoAcid['R'] == 6
    check byAminoAcid['F'] == 2
    check byAminoAcid['A'] == 4
    check byAminoAcid['*'] == 3

suite "codon usage":
  test "simple sequence":
    let usage = countCodonUsage("ATGAAATAA")
    check usage.totalCodons == 3
    check usage.entries.len == 64
    let aug = entry(usage, "AUG")
    check aug.count == 1
    check aug.aminoAcid == 'M'
    check aug.aminoAcidName == "Methionine"
    let uaa = entry(usage, "UAA")
    check uaa.aminoAcid == '*'
    check uaa.aminoAcidName == "Stop"

  test "unused codons have zero counts":
    let usage = countCodonUsage("ATG")
    check usage.totalCodons == 1
    check usage.entries.filterIt(it.count == 0).len == 63

  test "reading frame":
    let usage = countCodonUsage("ATGCAT", start = 1)
    check usage.totalCodons == 1
    check entry(usage, "UGC").count == 1

  test "DNA codons are normalized to RNA":
    let usage = countCodonUsage("ATGAAATAA")
    check entry(usage, "AUG").aminoAcid == 'M'
    check entry(usage, "AAA").count == 1
    check entry(usage, "UAA").count == 1

  test "ambiguous codons are skipped":
    let usage = countCodonUsage("ATGNNN")
    check usage.totalCodons == 1

  test "aggregation across sequences":
    let usage = countCodonUsageMany(["ATGAAA", "ATG"])
    check usage.totalCodons == 3
    check entry(usage, "AUG").count == 2
    check entry(usage, "AAA").count == 1

  test "empty and short sequences":
    check countCodonUsage("").totalCodons == 0
    check countCodonUsage("AT").totalCodons == 0

suite "FASTA parsing":
  test "multiple records with descriptions":
    let records = parseFasta(">a\nATGC\n>b desc here\nGG\nCC\n")
    check records.len == 2
    check records[0].id == "a"
    check records[0].description == ""
    check records[0].sequence == "ATGC"
    check records[1].id == "b"
    check records[1].description == "desc here"
    check records[1].sequence == "GGCC"

  test "whitespace inside sequence lines is removed":
    let records = parseFasta(">a\nA T\nGC\n")
    check records[0].sequence == "ATGC"

  test "comments and blank lines are ignored":
    let records = parseFasta("; a note\n\n>a\nAT\n")
    check records.len == 1
    check records[0].sequence == "AT"

  test "lowercase is preserved":
    check parseFasta(">a\natgc\n")[0].sequence == "atgc"

  test "empty input produces no records":
    check parseFasta("").len == 0

  test "sequence before a header is rejected":
    expect FastaError:
      discard parseFasta("ATGC")

  test "empty header is rejected":
    expect FastaError:
      discard parseFasta(">\nAT\n")

  test "record without sequence is rejected":
    expect FastaError:
      discard parseFasta(">a\n")

  test "readFasta reads the bundled fixture":
    let records = readFasta(fixturePath())
    check records.len == 2
    check records[0].id == "demo_orf"
    check records[0].sequence.len == 144
    check records[1].id == "promoter_fragment"
    check records[1].description.len > 0

suite "report":
  test "required sections are present":
    let records = parseFasta(">demo\nATGGCCTAA\n")
    let report = buildReport(records, source = "demo.fasta",
      generatedAt = "2026-08-02T12:00:00Z")
    check report.contains("Codon Atlas report")
    check report.contains("Source: demo.fasta")
    check report.contains("Generated: 2026-08-02T12:00:00Z")
    check report.contains("Records: 1")
    check report.contains("GC content: 44.4%")
    check report.contains("Total codons: 3")
    check report.contains("Methionine")
    check report.contains("Stop")

  test "output is deterministic":
    let records = parseFasta(">a\nATG\n")
    let first = buildReport(records, source = "a.fasta",
      generatedAt = "2026-08-02T12:00:00Z")
    let second = buildReport(records, source = "a.fasta",
      generatedAt = "2026-08-02T12:00:00Z")
    check first == second

  test "a generated timestamp is added by default":
    let report = buildReport(@[], source = "empty.fasta")
    check report.contains("Generated: ")

  test "k-mer section appears on request":
    let records = parseFasta(">a\nATGCAT\n")
    let report = buildReport(records, source = "a.fasta",
      options = ReportOptions(kmer: 2),
      generatedAt = "2026-08-02T12:00:00Z")
    check report.contains("K-mer counts (k = 2)")
    check report.contains("AT")

  test "k-mer section stays hidden by default":
    let records = parseFasta(">a\nATGCAT\n")
    let report = buildReport(records, source = "a.fasta",
      generatedAt = "2026-08-02T12:00:00Z")
    check not report.contains("K-mer counts")

  test "golden report matches the bundled fixture":
    let records = readFasta(fixturePath())
    let report = buildReport(records, source = "sample.fasta",
      options = ReportOptions(kmer: 3, kmerTop: 5),
      generatedAt = "2026-08-02T12:00:00Z")
    check report == readFile(goldenPath())
