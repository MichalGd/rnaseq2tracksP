#!/usr/bin/env bash
# ORIGIN: VERBATIM — syntax check all .sh scripts
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
FAIL=0
while IFS= read -r -d '' f; do
  if bash -n "$f" 2>/dev/null; then echo "  [OK] $(basename "$f")"
  else echo "  [FAIL] $(basename "$f")"; FAIL=$((FAIL+1)); fi
done < <(find "$REPO/scripts" -name "*.sh" -print0)
[[ $FAIL -eq 0 ]] && echo "All scripts pass." || { echo "$FAIL script(s) failed."; exit 1; }
