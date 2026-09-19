#!/usr/bin/env bash
set -euo pipefail
NP=${NP:-1}
runexe () {
  exe=$1; infile=$2; outfile=$3
  mkdir -p "$(dirname "$outfile")"
  echo "[$(date)] $exe < $infile > $outfile (NP=$NP)"
  if [ "$NP" -gt 1 ]; then
    if ! command -v mpirun >/dev/null 2>&1; then echo 'mpirun not found; set NP=1 or install OpenMPI'; exit 2; fi
    mpirun -np "$NP" "$exe" -in "$infile" > "$outfile"
  else
    "$exe" -in "$infile" > "$outfile"
  fi
  if grep -q 'JOB DONE' "$outfile"; then echo "DONE $outfile"; else echo "WARNING: JOB DONE not found in $outfile"; tail -30 "$outfile"; return 1; fi
}
