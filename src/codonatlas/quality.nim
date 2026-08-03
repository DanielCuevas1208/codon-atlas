## Quality score statistics for FASTQ sequences.
##
## Scores use the Sanger Phred+33 encoding. Each quality character maps to a
## score of `ord(c) - 33`. A valid score is an integer from 0 to 93.

const StandardQualityOffset* = 33
  ## Phred offset for Sanger FASTQ quality strings.

type
  QualitySummary* = object
    count*: int
      ## Number of quality scores.
    min*: int
      ## Lowest score.
    max*: int
      ## Highest score.
    mean*: float
      ## Average score. Zero when there are no scores.

proc qualityScore*(base: char, offset = StandardQualityOffset): int =
  ## Return the Phred quality score for a base character.
  ##
  ## Raises `ValueError` when the score falls outside 0 to 93.
  result = ord(base) - offset
  if result < 0 or result > 93:
    raise newException(ValueError, "quality score out of range: " & $result)

proc qualitySummary*(quality: string, offset = StandardQualityOffset): QualitySummary =
  ## Summarize a quality string.
  ##
  ## An empty string produces a summary with all values at zero.
  if quality.len == 0:
    return
  result.count = quality.len
  result.min = high(int)
  result.max = low(int)
  var total = 0
  for base in quality:
    let score = qualityScore(base, offset)
    total += score
    result.min = min(result.min, score)
    result.max = max(result.max, score)
  result.mean = total.float / result.count.float

proc basesBelowQuality*(quality: string, threshold: int,
                        offset = StandardQualityOffset): int =
  ## Count positions with a score strictly below `threshold`.
  for base in quality:
    if qualityScore(base, offset) < threshold:
      inc result
