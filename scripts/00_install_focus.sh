#!/bin/bash
# ============================================================
# FOCUS (pyfocus 0.802) for step 113, in its own virtual environment so it cannot disturb MetaXcan's.
#   FOCUS predates pandas 2 / numpy 1.24 / SQLAlchemy 2, so versions are pinned and two lines patched:
#   np.warnings -> warnings (exprref.py, ldref.py) and DataFrame.pivot keyword arguments (finemap.py).
#   Usage: bash scripts/00_install_focus.sh [python3.11 interpreter; default $FOCUS_PY]  ->  software/focus_venv/bin/focus
# ============================================================
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"
set -euo pipefail
PY=${1:-$FOCUS_PY}
"$PY" -c 'import sys; sys.exit(sys.version_info[:2] != (3, 11))' || {
  echo "FOCUS needs Python 3.11 (pandas 1.5.3 has no wheels for later versions): set FOCUS_PY in scripts/config.local.sh" >&2
  exit 1; }
V="$FECD_BASE/software/focus_venv"
"$PY" -m venv "$V"
"$V/bin/pip" install -q --upgrade pip
"$V/bin/pip" install -q "pyfocus==0.802" mygene "pandas==1.5.3" "numpy==1.26.4" "scipy==1.11.4" \
  "xarray==2023.1.0" "sqlalchemy==1.4.54"
P=$("$V/bin/python" -c "import pyfocus, os; print(os.path.dirname(pyfocus.__file__))")
for f in data/exprref.py data/ldref.py; do
  sed -i.bak 's/np\.warnings\./warnings./g' "$P/$f"
  grep -q '^import warnings' "$P/$f" || sed -i.bak '1s/^/import warnings\n/' "$P/$f"
done
sed -i.bak 's/attr_tmp.pivot("model_id", "attr_name", "value")/attr_tmp.pivot(index="model_id", columns="attr_name", values="value")/' "$P/finemap.py"
"$V/bin/focus" --help > /dev/null && echo "FOCUS installed: $V/bin/focus"
