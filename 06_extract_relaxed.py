#!/usr/bin/env python3
"""Extract the last ATOMIC_POSITIONS block from a QE relax output and
insert it into 06_interface_scf.template.in.
"""
from pathlib import Path
import json, re, sys
ROOT=Path(__file__).resolve().parent
out = ROOT/'outputs/05_interface_relax.out' if len(sys.argv)<2 else Path(sys.argv[1])
template = ROOT/'06_interface_scf.template.in'
dest = ROOT/'06_interface_scf.in'
meta=json.loads((ROOT/'interface_metadata.json').read_text())
nat=int(meta['nat'])
if not out.exists():
    raise SystemExit(f"Missing {out}")
lines=out.read_text(errors='ignore').splitlines()
starts=[i for i,l in enumerate(lines) if l.strip().startswith('ATOMIC_POSITIONS')]
if not starts:
    raise SystemExit('No ATOMIC_POSITIONS block found in relax output.')
# Find the last complete block with nat atom lines.
block=None
for i in reversed(starts):
    cand=[]
    for l in lines[i+1:i+1+nat]:
        p=l.split()
        if len(p)<4 or p[0] not in {'Ti','N','Hf','O'}:
            cand=[]; break
        try:
            float(p[1]);float(p[2]);float(p[3])
        except ValueError:
            cand=[];break
        cand.append((p[0],float(p[1]),float(p[2]),float(p[3])))
    if len(cand)==nat:
        block=cand; break
if block is None:
    raise SystemExit('Could not find a complete final coordinate block.')
text='ATOMIC_POSITIONS angstrom\n'+'\n'.join(
    f"{sp:2s} {x:14.8f} {y:14.8f} {z:14.8f}" for sp,x,y,z in block)
t=template.read_text()
if '__ATOMIC_POSITIONS__' not in t:
    raise SystemExit('Template placeholder not found.')
dest.write_text(t.replace('__ATOMIC_POSITIONS__',text))
print(f'Wrote {dest.name} with {nat} relaxed coordinates.')
