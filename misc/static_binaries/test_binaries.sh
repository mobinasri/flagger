#!/bin/bash
# Checks prebuilt bam2cov and hmm_flagger binaries: runs bam2cov on a test BAM from
# programs/tests and compares its output with the expected output.
# Run from the repo root:
#   bash misc/static_binaries/test_binaries.sh <dir_with_binaries>
set -eu
BIN=$(cd "$1" && pwd)
# the annotation JSON lists BED paths relative to programs/
cd programs
T=tests/test_files/bam2cov
OUT=$(mktemp -d)
"${BIN}/hmm_flagger" -h > /dev/null
"${BIN}/bam2cov" -h > /dev/null
"${BIN}/bam2cov" \
  --bam ${T}/ptBlock_bam2cov_test.bam \
  --annotationJson ${T}/ptBlock_bam2cov_test.json \
  --output ${OUT}/test.cov \
  --runBiasDetection --minBiasCoverage 1 --minBiasLength 1 \
  --minAlignmentLength 0 --baselineAnnotation annot_1
# a checksum instead of cmp, which some minimal images lack; macOS has shasum instead of sha256sum
if command -v sha256sum > /dev/null; then sum() { sha256sum; }; else sum() { shasum -a 256; }; fi
[ "$(sum < ${OUT}/test.cov)" = "$(sum < ${T}/ptBlock_bam2cov_test.truth.cov)" ] ||
  { echo "::error::bam2cov output differs from the truth file"; exit 1; }
rm -rf "${OUT}"
echo "bam2cov output matches the truth file"
