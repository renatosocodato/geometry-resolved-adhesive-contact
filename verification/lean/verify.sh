#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "$0")"
lake env lean Contact.lean > check.log 2>&1 || { cat check.log; exit 1; }
python3 - <<'PY'
from pathlib import Path
import re
source = Path("Contact.lean").read_text()
assert not re.search(r"\b(?:axiom|sorry|admit|unsafe|native_decide)\b", source)
expected = set(re.findall(r"^theorem (\w+)", source, re.M))
log = Path("check.log").read_text()
assert "error:" not in log
lines = re.findall(r"'MRC\.(\w+)' depends on axioms: \[([^\]]*)\]", log)
allow = {"propext", "Classical.choice", "Quot.sound"}
assert expected == {theorem for theorem, _ in lines}
for theorem, printed in lines:
    assert set(printed.split(", ")) <= allow, (theorem, printed)
print(f"Lean kernel accepted all {len(lines)} public theorems with only standard axioms.")
PY
