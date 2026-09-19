#!/usr/bin/env bash
set -euo pipefail
mkdir -p pseudo
base='https://pseudopotentials.quantum-espresso.org/upf_files'
files=(
  'Ti.pbe-spn-rrkjus_psl.1.0.0.UPF'
  'N.pbe-n-rrkjus_psl.1.0.0.UPF'
  'Hf.pbe-spn-rrkjus_psl.1.0.0.UPF'
  'O.pbe-n-rrkjus_psl.1.0.0.UPF'
)
for f in "${files[@]}"; do
  echo "Downloading $f"
  if command -v wget >/dev/null 2>&1; then
    wget -q --show-progress "$base/$f" -O "pseudo/$f"
  else
    curl -L "$base/$f" -o "pseudo/$f"
  fi
done
ls -lh pseudo/*.UPF
