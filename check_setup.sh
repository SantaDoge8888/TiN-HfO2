#!/usr/bin/env bash
set -u
bad=0
for x in pw.x pp.x projwfc.x python3; do
  if command -v "$x" >/dev/null 2>&1; then echo "OK   $x -> $(command -v "$x")"; else echo "MISS $x"; bad=1; fi
done
for f in \
 pseudo/Ti.pbe-spn-rrkjus_psl.1.0.0.UPF \
 pseudo/N.pbe-n-rrkjus_psl.1.0.0.UPF \
 pseudo/Hf.pbe-spn-rrkjus_psl.1.0.0.UPF \
 pseudo/O.pbe-n-rrkjus_psl.1.0.0.UPF; do
  if [ -s "$f" ]; then echo "OK   $f ($(du -h "$f" | cut -f1))"; else echo "MISS $f"; bad=1; fi
done
python3 - <<'PY' || bad=1
try:
 import numpy
 print('OK   numpy',numpy.__version__)
except Exception as e:
 print('MISS numpy:',e)
 raise
PY
exit $bad
