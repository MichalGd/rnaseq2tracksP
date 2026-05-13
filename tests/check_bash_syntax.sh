#!/usr/bin/env bash
# ORIGIN: VERBATIM COPY — fastq2tracks/tests/check_bash_syntax.sh
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PASS=0; FAIL=0
while IFS= read -r -d '' f; do
  if bash -n "$f" 2>/dev/null; then echo "PASS: $f"; PASS=$((PASS+1))
  else echo "FAIL: $f"; FAIL=$((FAIL+1)); fi
done < <(find "$REPO/scripts" -name "*.sh" -print0)
echo "Results: $PASS passed, $FAIL failed"; [[ $FAIL -eq 0 ]]
