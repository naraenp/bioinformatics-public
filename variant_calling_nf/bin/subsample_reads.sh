#!/usr/bin/env bash
# Deterministically subsample a paired-end FASTQ to <n> read pairs with seqtk.
# The SAME seed is used for both mates so pairing is preserved. A non-positive
# <n> (or seqtk being absent) passes the reads through unchanged.
#
# Usage: subsample_reads.sh <in_R1> <in_R2> <out_R1> <out_R2> <n> <seed>
set -euo pipefail

in1="$1"; in2="$2"; out1="$3"; out2="$4"; n="${5:-0}"; seed="${6:-42}"

if [[ "$n" -le 0 ]] || ! command -v seqtk >/dev/null 2>&1; then
    cp "$in1" "$out1"
    cp "$in2" "$out2"
    echo "subsample: pass-through ($(basename "$in1"))"
    exit 0
fi

seqtk sample -s"$seed" "$in1" "$n" | gzip > "$out1"
seqtk sample -s"$seed" "$in2" "$n" | gzip > "$out2"
echo "subsample: $n pairs @ seed=$seed ($(basename "$in1"))"
