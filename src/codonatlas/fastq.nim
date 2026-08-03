## FASTQ reading and parsing.
##
## The parser accepts the standard Sanger FASTQ format. Each record has a
## header line, a sequence, a `+` separator, and a quality line. Sequence and
## quality lines may wrap over multiple lines. Their combined lengths must
## match. Quality scores use the Phred+33 encoding.

import std/strutils

type
  FastqError* = object of CatchableError
    ## Raised when the FASTQ text is malformed.

  FastqRecord* = object
    id*: string
      ## First token of the header line, after the `@`.
    description*: string
      ## Remainder of the header line, if any.
    sequence*: string
      ## Concatenated sequence bases, with whitespace removed.
    quality*: string
      ## Concatenated quality characters, with whitespace removed.

proc stripWhitespace(line: string): string =
  result = newStringOfCap(line.len)
  for ch in line:
    if ch notin {' ', '\t', '\r', '\n'}:
      result.add ch

proc parseHeader(line: string, lineNo: int): tuple[id, description: string] =
  let header = line[1 .. ^1].strip
  if header.len == 0:
    raise newException(FastqError, "empty sequence header on line " & $lineNo)
  let fields = header.splitWhitespace()
  result.id = fields[0]
  result.description = fields[1 .. ^1].join(" ")

proc parseFastq*(content: string): seq[FastqRecord] =
  ## Parse FASTQ text into a list of records.
  ##
  ## Raises `FastqError` when the text is malformed. An empty input produces
  ## an empty list.
  let lines = splitLines(content)
  var index = 0
  while index < lines.len:
    if lines[index].strip.len == 0:
      inc index
      continue
    if lines[index][0] != '@':
      raise newException(FastqError,
        "expected a header starting with '@' (line " & $(index + 1) & ")")
    var record = FastqRecord()
    let (id, description) = parseHeader(lines[index], index + 1)
    record.id = id
    record.description = description
    inc index

    var sequenceParts: seq[string]
    while index < lines.len:
      let line = lines[index]
      if line.len > 0 and line[0] == '+':
        break
      if line.len > 0:
        sequenceParts.add stripWhitespace(line)
      inc index
    if sequenceParts.len == 0:
      raise newException(FastqError,
        "record '" & record.id & "' has no sequence")
    if index >= lines.len:
      raise newException(FastqError,
        "record '" & record.id & "' is missing the '+' separator")
    record.sequence = sequenceParts.join("")
    inc index

    if index >= lines.len:
      raise newException(FastqError,
        "record '" & record.id & "' has no quality line")
    var quality = ""
    while index < lines.len and quality.len < record.sequence.len:
      if lines[index].len > 0:
        quality.add stripWhitespace(lines[index])
      inc index
    if quality.len < record.sequence.len:
      raise newException(FastqError,
        "quality is shorter than the sequence in record '" & record.id & "'")
    if quality.len > record.sequence.len:
      raise newException(FastqError,
        "quality is longer than the sequence in record '" & record.id & "'")
    record.quality = quality
    result.add record

proc readFastq*(path: string): seq[FastqRecord] =
  ## Read a FASTQ file and parse its records.
  ##
  ## Raises `IOError` when the file cannot be read. Raises `FastqError` when
  ## the file content is malformed.
  parseFastq(readFile(path))
