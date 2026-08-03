## The standard genetic code.
##
## Codons use RNA bases (U instead of T). The table maps all 64 RNA codons to
## a one-letter amino acid code. A stop codon maps to `*`.

import std/algorithm
import std/strutils
import std/tables

let StandardCode* = {
  "UUU": 'F', "UUC": 'F', "UUA": 'L', "UUG": 'L',
  "UCU": 'S', "UCC": 'S', "UCA": 'S', "UCG": 'S',
  "UAU": 'Y', "UAC": 'Y', "UAA": '*', "UAG": '*',
  "UGU": 'C', "UGC": 'C', "UGA": '*', "UGG": 'W',
  "CUU": 'L', "CUC": 'L', "CUA": 'L', "CUG": 'L',
  "CCU": 'P', "CCC": 'P', "CCA": 'P', "CCG": 'P',
  "CAU": 'H', "CAC": 'H', "CAA": 'Q', "CAG": 'Q',
  "CGU": 'R', "CGC": 'R', "CGA": 'R', "CGG": 'R',
  "AUU": 'I', "AUC": 'I', "AUA": 'I', "AUG": 'M',
  "ACU": 'T', "ACC": 'T', "ACA": 'T', "ACG": 'T',
  "AAU": 'N', "AAC": 'N', "AAA": 'K', "AAG": 'K',
  "AGU": 'S', "AGC": 'S', "AGA": 'R', "AGG": 'R',
  "GUU": 'V', "GUC": 'V', "GUA": 'V', "GUG": 'V',
  "GCU": 'A', "GCC": 'A', "GCA": 'A', "GCG": 'A',
  "GAU": 'D', "GAC": 'D', "GAA": 'E', "GAG": 'E',
  "GGU": 'G', "GGC": 'G', "GGA": 'G', "GGG": 'G'
}.toTable

const AminoAcidNames* = {
  'A': "Alanine", 'R': "Arginine", 'N': "Asparagine", 'D': "Aspartic acid",
  'C': "Cysteine", 'Q': "Glutamine", 'E': "Glutamic acid", 'G': "Glycine",
  'H': "Histidine", 'I': "Isoleucine", 'L': "Leucine", 'K': "Lysine",
  'M': "Methionine", 'F': "Phenylalanine", 'P': "Proline", 'S': "Serine",
  'T': "Threonine", 'W': "Tryptophan", 'Y': "Tyrosine", 'V': "Valine"
}.toTable

var sortedCodons*: seq[string]
for codon, _ in StandardCode:
  sortedCodons.add codon
sortedCodons.sort()

proc dnaToRna*(codon: string): string =
  ## Convert a codon to RNA by replacing every T with U.
  result = toUpperAscii(codon)
  result = result.replace('T', 'U')

proc aminoAcid*(rnaCodon: string): char =
  ## Return the one-letter amino acid for an RNA codon.
  ##
  ## A stop codon returns `*`. An unknown codon returns `?`.
  StandardCode.getOrDefault(toUpperAscii(rnaCodon), '?')

proc codonAminoAcid*(codon: string): char =
  ## Return the one-letter amino acid for a DNA or RNA codon.
  ##
  ## The input may use T or U. A stop codon returns `*`.
  aminoAcid(dnaToRna(codon))
