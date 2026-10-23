#!/bin/sh
set -eu
structure_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
cd "$structure_dir"
if [ -z "${CHIMERAX_BIN:-}" ]; then
  if command -v ChimeraX >/dev/null 2>&1; then
    CHIMERAX_BIN=$(command -v ChimeraX)
  elif [ -x /Applications/ChimeraX-1.12.app/Contents/MacOS/ChimeraX ]; then
    CHIMERAX_BIN=/Applications/ChimeraX-1.12.app/Contents/MacOS/ChimeraX
  else
    printf "%s\n" "Set CHIMERAX_BIN to your ChimeraX executable." >&2
    exit 1
  fi
fi
FIGURE_PYTHON=${FIGURE_PYTHON:-python3}
"$CHIMERAX_BIN" --nogui --exit --script "$structure_dir/calculate_rmsd_sata_2qir.py"
"$CHIMERAX_BIN" --exit --script "$structure_dir/render_chimerax.py"
"$FIGURE_PYTHON" "$structure_dir/compose_figures.py"
