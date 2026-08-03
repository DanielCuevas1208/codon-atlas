version       = "0.2.0"
author        = "Codon Atlas contributors"
description   = "Composition statistics for DNA and RNA sequence files"
license       = "MIT"
srcDir        = "src"

bin           = @["codonatlas_cli"]

requires "nim >= 1.6.0"

task test, "Run the unit test suite":
  exec "nim c -r --path:src --outdir:build/tests tests/test_codonatlas.nim"
  exec "nim c -r --path:src --outdir:build/tests tests/test_fastq.nim"

task check, "Typecheck all modules":
  exec "nim check --path:src src/codonatlas.nim"
  exec "nim check --path:src src/codonatlas_cli.nim"
  exec "nim check --path:src tests/test_codonatlas.nim"
  exec "nim check --path:src tests/test_fastq.nim"
