## FASTA reading and parsing.
##
## The parser accepts standard FASTA text. Blank lines and lines that start
## with `;` are ignored. Sequence lines may wrap over many lines. Spaces and
## tabs inside a sequence line are removed.

import std/strutils

type
  FastaError* = object of CatchableError
    ## Raised when the FASTA text is malformed.

  Record* = object
    id*: string
      ## First token of the header line, after the `>`.
    description*: string
      ## Remainder of the header line, if any.
    sequence*: string
      ## Concatenated sequence bases, with whitespace removed.

proc cleanSequence(line: string): string =
  result = newStringOfCap(line.len)
  for base in line:
    if base notin {' ', '\t', '\r', '\n'}:
      result.add base

proc parseFasta*(content: string): seq[Record] =
  ## Parse FASTA text into a list of records.
  ##
  ## Raises `FastaError` when the text is malformed. An empty input produces
  ## an empty list.
  var
    lineNo = 0
    inRecord = false
  for rawLine in splitLines(content):
    inc lineNo
    if rawLine.len == 0 or rawLine[0] == ';':
      continue
    if rawLine[0] == '>':
      if inRecord and result.len > 0 and result[^1].sequence.len == 0:
        raise newException(FastaError,
          "record '" & result[^1].id & "' has no sequence (line " & $lineNo & ")")
      let header = rawLine[1 .. ^1].strip
      if header.len == 0:
        raise newException(FastaError,
          "empty sequence header on line " & $lineNo)
      let fields = header.splitWhitespace()
      result.add Record(
        id: fields[0],
        description: fields[1 .. ^1].join(" "))
      inRecord = true
    else:
      if not inRecord:
        raise newException(FastaError,
          "sequence data before the first header (line " & $lineNo & ")")
      result[^1].sequence.add cleanSequence(rawLine)
  if inRecord and result.len > 0 and result[^1].sequence.len == 0:
    raise newException(FastaError,
      "record '" & result[^1].id & "' has no sequence")

proc readFasta*(path: string): seq[Record] =
  ## Read a FASTA file and parse its records.
  ##
  ## Raises `IOError` when the file cannot be read. Raises `FastaError` when
  ## the file content is malformed.
  parseFasta(readFile(path))
