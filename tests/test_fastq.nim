import std/os
import std/strutils
import std/unittest

import codonatlas

proc fixturePath(): string =
  currentSourcePath().parentDir() / "data" / "sample.fastq"

proc goldenPath(): string =
  currentSourcePath().parentDir() / "data" / "expected_fastq_report.txt"

suite "quality scores":
  test "Phred+33 encoding":
    check qualityScore('!') == 0
    check qualityScore('I') == 40
    check qualityScore('~') == 93
    check qualityScore('@', offset = 64) == 0

  test "out-of-range scores are rejected":
    expect ValueError:
      discard qualityScore(' ')

suite "quality summary":
  test "known strings":
    let summary = qualitySummary("IHH")
    check summary.count == 3
    check summary.min == 39
    check summary.max == 40
    check abs(summary.mean - 39.3333) < 0.001

  test "empty string yields zeros":
    let summary = qualitySummary("")
    check summary.count == 0
    check summary.min == 0
    check summary.max == 0
    check summary.mean == 0.0

suite "bases below threshold":
  test "strictly below the threshold":
    check basesBelowQuality("IHH", 40) == 2
    check basesBelowQuality("IHH", 39) == 0
    check basesBelowQuality("", 20) == 0

suite "FASTQ parsing":
  test "multiple records with descriptions":
    let records = parseFastq("@a\nAT\n+\nII\n@b note here\nGC\n+\n!!\n")
    check records.len == 2
    check records[0].id == "a"
    check records[0].description == ""
    check records[0].sequence == "AT"
    check records[0].quality == "II"
    check records[1].id == "b"
    check records[1].description == "note here"
    check records[1].sequence == "GC"
    check records[1].quality == "!!"

  test "wrapped sequence and quality lines":
    let records = parseFastq("@a\nAT\nGC\n+\nII\n!!\n")
    check records[0].sequence == "ATGC"
    check records[0].quality == "II!!"

  test "blank lines between records are ignored":
    let records = parseFastq("@a\nAT\n+\nII\n\n@b\nGC\n+\n!!\n")
    check records.len == 2

  test "whitespace inside lines is removed":
    let records = parseFastq("@a\nA T\nGC\n+\nI I\n!!\n")
    check records[0].sequence == "ATGC"
    check records[0].quality == "II!!"

  test "lowercase is preserved":
    let records = parseFastq("@a\natgc\n+\nIIII\n")
    check records[0].sequence == "atgc"

  test "the plus line may carry text":
    let records = parseFastq("@a\nAT\n+ some text\nII\n")
    check records[0].sequence == "AT"
    check records[0].quality == "II"

  test "empty input produces no records":
    check parseFastq("").len == 0

  test "sequence before a header is rejected":
    expect FastqError:
      discard parseFastq("ATGC\n+\nIIII\n")

  test "empty header is rejected":
    expect FastqError:
      discard parseFastq("@\nAT\n+\nII\n")

  test "record without sequence is rejected":
    expect FastqError:
      discard parseFastq("@a\n+\nII\n")

  test "missing plus separator is rejected":
    expect FastqError:
      discard parseFastq("@a\nAT\nII\n")

  test "missing quality line is rejected":
    expect FastqError:
      discard parseFastq("@a\nAT\n+\n")

  test "quality shorter than sequence is rejected":
    expect FastqError:
      discard parseFastq("@a\nAT\n+\nI\n")

  test "quality longer than sequence is rejected":
    expect FastqError:
      discard parseFastq("@a\nAT\n+\nIII\n")

  test "readFastq reads the bundled fixture":
    let records = readFastq(fixturePath())
    check records.len == 2
    check records[0].id == "read1"
    check records[0].sequence.len == 23
    check records[0].quality.len == 23
    check records[1].id == "read2"
    check records[1].description.len > 0

suite "FASTQ report":
  test "required sections are present":
    let records = readFastq(fixturePath())
    let report = buildFastqReport(records, source = "sample.fastq",
      generatedAt = "2026-08-02T12:00:00Z")
    check report.contains("Codon Atlas report")
    check report.contains("Source: sample.fastq")
    check report.contains("Generated: 2026-08-02T12:00:00Z")
    check report.contains("Records: 2")
    check report.contains("Total canonical bases: 35")
    check report.contains("Quality: mean 38.3, min 37, max 40")

  test "the below-threshold line appears on request":
    let records = readFastq(fixturePath())
    let report = buildFastqReport(records, source = "sample.fastq",
      options = ReportOptions(qualityThreshold: 20),
      generatedAt = "2026-08-02T12:00:00Z")
    check report.contains("Bases below Q20: 12")

  test "the below-threshold line stays hidden by default":
    let records = readFastq(fixturePath())
    let report = buildFastqReport(records, source = "sample.fastq",
      generatedAt = "2026-08-02T12:00:00Z")
    check not report.contains("Bases below Q20")

  test "k-mer section appears on request":
    let records = readFastq(fixturePath())
    let report = buildFastqReport(records, source = "sample.fastq",
      options = ReportOptions(kmer: 3),
      generatedAt = "2026-08-02T12:00:00Z")
    check report.contains("K-mer counts (k = 3)")

  test "output is deterministic":
    let records = readFastq(fixturePath())
    let first = buildFastqReport(records, source = "sample.fastq",
      generatedAt = "2026-08-02T12:00:00Z")
    let second = buildFastqReport(records, source = "sample.fastq",
      generatedAt = "2026-08-02T12:00:00Z")
    check first == second

  test "golden report matches the bundled fixture":
    let records = readFastq(fixturePath())
    let report = buildFastqReport(records, source = "sample.fastq",
      options = ReportOptions(kmer: 3, kmerTop: 3, qualityThreshold: 20),
      generatedAt = "2026-08-02T12:00:00Z")
    check report == readFile(goldenPath())
