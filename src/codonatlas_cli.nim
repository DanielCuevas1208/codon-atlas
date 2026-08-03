## Command-line interface for Codon Atlas.

import std/os
import std/parseopt
import std/strutils

import codonatlas

const
  Version = "0.1.0"
  ToolName = "codonatlas_cli"

proc printUsage() =
  echo ToolName & " " & Version & " - composition statistics for sequence files"
  echo ""
  echo "Usage:"
  echo "  " & ToolName & " <file.fasta> [options]"
  echo ""
  echo "Options:"
  echo "  -o, --out=<path>   Write the report to a file."
  echo "  -k, --kmer=<n>     Add a k-mer count section for k = n."
  echo "  -t, --top=<n>      Show only the top n k-mers."
  echo "  -h, --help         Show this help."
  echo "  -v, --version      Show the version."

proc fail(message: string) =
  stderr.writeLine("Error: " & message)
  quit(1)

proc main() =
  var
    input: seq[string]
    outFile = ""
    kmer = 0
    top = 0
    showHelp = false
    showVersion = false

  var parser = initOptParser()
  for kind, key, val in parser.getopt():
    case kind
    of cmdArgument:
      input.add key
    of cmdLongOption, cmdShortOption:
      case key
      of "out", "o":
        outFile = val
      of "kmer", "k":
        try:
          kmer = parseInt(val)
        except ValueError:
          fail("invalid k-mer size: " & val)
      of "top", "t":
        try:
          top = parseInt(val)
        except ValueError:
          fail("invalid top value: " & val)
      of "help", "h":
        showHelp = true
      of "version", "v":
        showVersion = true
      else:
        fail("unknown option: " & key)
    of cmdEnd:
      discard

  if showHelp:
    printUsage()
    quit(0)
  if showVersion:
    echo Version
    quit(0)
  if input.len != 1:
    stderr.writeLine("Error: provide exactly one FASTA file.")
    printUsage()
    quit(1)

  let inFile = input[0]
  if not fileExists(inFile):
    fail("file not found: " & inFile)

  var records: seq[Record]
  try:
    records = readFasta(inFile)
  except FastaError as exc:
    fail("invalid FASTA in " & inFile & ": " & exc.msg)

  let options = ReportOptions(kmer: kmer, kmerTop: top)
  let report = buildReport(records, source = inFile, options = options)

  if outFile.len > 0:
    writeFile(outFile, report)
    echo "Report written to " & outFile
  else:
    stdout.write(report)

when isMainModule:
  main()
