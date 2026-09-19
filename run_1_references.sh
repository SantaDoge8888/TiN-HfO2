#!/usr/bin/env bash
set -euo pipefail
source ./run_common.sh
mkdir -p tmp/bulk_tin tmp/bulk_hfo2 outputs
runexe pw.x 02_bulk_tin_strained_scf.in outputs/02_bulk_tin_strained_scf.out
runexe pw.x 03_bulk_hfo2_scf.in outputs/03_bulk_hfo2_scf.out
runexe pp.x 08_pp_bulk_tin.in outputs/08_pp_bulk_tin.out
runexe pp.x 09_pp_bulk_hfo2.in outputs/09_pp_bulk_hfo2.out
