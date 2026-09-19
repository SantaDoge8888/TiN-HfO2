#!/usr/bin/env bash
set -euo pipefail
source ./run_common.sh
mkdir -p tmp/interface_smoke outputs
runexe pw.x 04_interface_smoke.in outputs/04_interface_smoke.out
