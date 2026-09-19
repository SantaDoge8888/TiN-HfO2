#!/usr/bin/env bash
set -euo pipefail
source ./run_common.sh
mkdir -p tmp/interface outputs
runexe pw.x 05_interface_relax.in outputs/05_interface_relax.out
python3 06_extract_relaxed.py
