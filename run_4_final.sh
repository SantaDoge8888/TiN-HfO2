#!/usr/bin/env bash
set -euo pipefail
source ./run_common.sh
if [ ! -f 06_interface_scf.in ]; then echo 'Missing 06_interface_scf.in. Run run_3_relax.sh first.'; exit 2; fi
runexe pw.x 06_interface_scf.in outputs/06_interface_scf.out
runexe projwfc.x 07_interface_projwfc.in outputs/07_interface_projwfc.out
runexe pp.x 10_pp_interface.in outputs/10_pp_interface.out
python3 11_analyze.py
