#!/usr/bin/env python3
"""Analyze TiN/HfO2 band alignment and effective work function.

Reports two literature-style VBO estimates when data are available:
  1) WM-like planar-potential lineup using bulk references.
  2) PDOS estimate from central HfO2 atoms, following the paper's preferred
     approach for short dielectric slabs.

W_eff = BG_exp(HfO2) + EA_exp(HfO2) - VBO
with 5.7 eV and 2.9 eV as in Fonseca & Knizhnik (2006).

Also reports the paper-style heuristic VBO scaling by Eg_exp/Eg_DFT.
This is a heuristic comparison, not a GW calculation.
"""
from __future__ import annotations
from pathlib import Path
import json, re, glob, math
import numpy as np

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'outputs'
RY_TO_EV=13.605693122994
BOHR_TO_A=0.529177210903
BG_EXP=5.7
EA_EXP=2.9


def last_float(pattern, text):
    vals=re.findall(pattern,text,re.I|re.M)
    if not vals: return None
    v=vals[-1]
    if isinstance(v,tuple): v=v[0]
    return float(v)


def parse_fermi(path):
    txt=Path(path).read_text(errors='ignore')
    return last_float(r'the Fermi energy is\s*([-+0-9.eEdD]+)\s*ev',txt)


def parse_edges(path):
    txt=Path(path).read_text(errors='ignore')
    matches=re.findall(r'highest occupied, lowest unoccupied level \(ev\):\s*([-+0-9.eEdD]+)\s+([-+0-9.eEdD]+)',txt,re.I)
    if matches:
        v,c=matches[-1]
        return float(v.replace('D','E').replace('d','e')),float(c.replace('D','E').replace('d','e'))
    v=last_float(r'highest occupied level \(ev\):\s*([-+0-9.eEdD]+)',txt)
    return (v,None)


def read_cube(path):
    lines=Path(path).read_text(errors='ignore').splitlines()
    if len(lines)<10: raise ValueError(f'{path} is not a valid cube file')
    p=lines[2].split(); nat=abs(int(p[0])); origin=np.array(list(map(float,p[1:4])))
    dims=[]; axes=[]
    for i in range(3):
        q=lines[3+i].split(); dims.append(abs(int(q[0]))); axes.append(np.array(list(map(float,q[1:4]))))
    start=6+nat
    vals=np.fromstring(' '.join(lines[start:]),sep=' ')
    n=np.prod(dims)
    if vals.size<n: raise ValueError(f'Cube has {vals.size} values, expected {n}')
    arr=vals[:n].reshape(tuple(dims))
    return origin,np.array(axes),arr


def cube_mean_ev(path):
    _,_,a=read_cube(path)
    return float(a.mean()*RY_TO_EV)


def cube_planar_z_ev(path):
    origin,axes,a=read_cube(path)
    planar=a.mean(axis=(0,1))*RY_TO_EV
    # For our interface cell the third FFT vector is along Cartesian z.
    z0=origin[2]*BOHR_TO_A
    dz=axes[2,2]*BOHR_TO_A
    z=z0+np.arange(a.shape[2])*dz
    return z,planar


def region_mean(z,v,lo,hi):
    m=(z>=lo)&(z<=hi)
    if m.sum()<3: raise ValueError(f'Only {m.sum()} grid points in region {lo}:{hi} A')
    return float(v[m].mean()),int(m.sum())


def load_atom_pdos(atom_ids):
    energies=None; total=None; files_used=[]
    for idx in atom_ids:
        pats=[str(OUT/f'interface_pdos.pdos_atm#{idx}(*)_wfc#*')]
        found=[]
        for p in pats: found.extend(glob.glob(p))
        for fn in sorted(set(found)):
            rows=[]
            for l in Path(fn).read_text(errors='ignore').splitlines():
                if not l.strip() or l.lstrip().startswith('#'): continue
                q=l.split()
                try: rows.append((float(q[0]),float(q[1])))
                except (ValueError,IndexError): pass
            if not rows: continue
            e=np.array([r[0] for r in rows]); d=np.array([r[1] for r in rows])
            if energies is None:
                energies=e; total=np.zeros_like(d)
            if len(e)!=len(energies) or np.max(np.abs(e-energies))>1e-6: continue
            total+=d; files_used.append(fn)
    return energies,total,files_used


def pdos_vbm(energies,dos,ef,frac):
    if energies is None: return None
    occ=energies < ef-0.03
    if occ.sum()<10: return None
    mx=float(np.max(dos[occ])); threshold=max(frac*mx,1e-8)
    # Search below EF for the highest robust oxide DOS. Requiring 3 nearby
    # grid points above threshold suppresses isolated numerical spikes.
    inds=np.where(occ & (dos>threshold))[0]
    good=[]
    S=set(int(i) for i in inds)
    for i in inds:
        if sum((j in S) for j in range(max(0,i-2),min(len(dos),i+3)))>=3:
            good.append(i)
    if not good: return None
    return float(energies[max(good)]),threshold


def main():
    meta=json.loads((ROOT/'interface_metadata.json').read_text())
    report=[]
    report.append('TiN(111) / m-HfO2(001) effective-work-function analysis')
    report.append('='*68)
    report.append(f"Model: {meta['model']}")
    report.append('')

    ef_tin=parse_fermi(OUT/'02_bulk_tin_strained_scf.out') if (OUT/'02_bulk_tin_strained_scf.out').exists() else None
    ev_hf=ec_hf=None
    if (OUT/'03_bulk_hfo2_scf.out').exists(): ev_hf,ec_hf=parse_edges(OUT/'03_bulk_hfo2_scf.out')
    eg_dft=(ec_hf-ev_hf) if ev_hf is not None and ec_hf is not None else None
    ef_int=parse_fermi(OUT/'06_interface_scf.out') if (OUT/'06_interface_scf.out').exists() else None

    report.append('Bulk/interface electronic references (QE energy zero is arbitrary):')
    report.append(f'  strained TiN bulk EF      = {ef_tin:.6f} eV' if ef_tin is not None else '  strained TiN bulk EF      = MISSING')
    report.append(f'  bulk HfO2 VBM             = {ev_hf:.6f} eV' if ev_hf is not None else '  bulk HfO2 VBM             = MISSING')
    report.append(f'  bulk HfO2 CBM             = {ec_hf:.6f} eV' if ec_hf is not None else '  bulk HfO2 CBM             = MISSING')
    report.append(f'  bulk HfO2 PBE gap         = {eg_dft:.6f} eV' if eg_dft is not None else '  bulk HfO2 PBE gap         = MISSING')
    report.append(f'  interface EF              = {ef_int:.6f} eV' if ef_int is not None else '  interface EF              = MISSING')
    report.append('')

    # WM-like lineup
    needed=[OUT/'bulk_tin_potential.cube',OUT/'bulk_hfo2_potential.cube',OUT/'interface_potential.cube']
    if all(p.exists() for p in needed) and ef_tin is not None and ev_hf is not None:
        vt_bulk=cube_mean_ev(needed[0]); vh_bulk=cube_mean_ev(needed[1])
        z,vi=cube_planar_z_ev(needed[2])
        tr=meta['potential_regions_A']['tin']; hr=meta['potential_regions_A']['hfo2']
        vt_int,nt=region_mean(z,vi,*tr); vh_int,nh=region_mean(z,vi,*hr)
        d_t=ef_tin-vt_bulk
        d_h=ev_hf-vh_bulk
        vbo_wm=d_t-d_h+(vt_int-vh_int)
        weff_wm=BG_EXP+EA_EXP-vbo_wm
        report.append('WM-like planar-potential lineup:')
        report.append(f'  <V> TiN bulk              = {vt_bulk:.6f} eV')
        report.append(f'  <V> HfO2 bulk             = {vh_bulk:.6f} eV')
        report.append(f'  <V> TiN interface center  = {vt_int:.6f} eV  ({nt} z-grid points)')
        report.append(f'  <V> HfO2 interface center = {vh_int:.6f} eV  ({nh} z-grid points)')
        report.append(f'  VBO_WM = EF(TiN)-VBM(HfO2)= {vbo_wm:.4f} eV')
        report.append(f'  W_eff_WM (Eg=5.7, EA=2.9)= {weff_wm:.4f} eV')
        if eg_dft and eg_dft>0:
            vbo_scaled=vbo_wm*(BG_EXP/eg_dft)
            report.append(f'  paper-style scaled VBO     = {vbo_scaled:.4f} eV  [heuristic]')
            report.append(f'  paper-style scaled W_eff   = {BG_EXP+EA_EXP-vbo_scaled:.4f} eV  [heuristic]')
        report.append('')
    else:
        report.append('WM-like lineup: not evaluated (missing cube/reference files).\n')

    # PDOS method
    if ef_int is not None:
        e,d,used=load_atom_pdos(meta['central_hfo2_atom_ids'])
        if e is not None and used:
            report.append(f'PDOS method using {len(used)} central-HfO2 orbital files:')
            for frac in (0.005,0.01,0.02):
                res=pdos_vbm(e,d,ef_int,frac)
                if res:
                    vbm,thr=res; vbo=ef_int-vbm; weff=BG_EXP+EA_EXP-vbo
                    report.append(f'  threshold {100*frac:4.1f}%: VBM={vbm:.4f} eV, VBO={vbo:.4f} eV, W_eff={weff:.4f} eV')
                    if eg_dft and eg_dft>0:
                        sv=vbo*(BG_EXP/eg_dft); sw=BG_EXP+EA_EXP-sv
                        report.append(f'                   scaled VBO={sv:.4f} eV, scaled W_eff={sw:.4f} eV [heuristic]')
            # export combined central oxide PDOS relative to EF
            np.savetxt(OUT/'central_hfo2_pdos.dat',np.column_stack([e-ef_int,d]),header='E-EF_eV  summed_central_HfO2_PDOS_states_per_eV')
            report.append('  Combined central-oxide PDOS written to outputs/central_hfo2_pdos.dat')
            report.append('')
        else:
            report.append('PDOS method: no interface_pdos files found.\n')

    report.append('Reference values from the 2006 paper (comparison only):')
    report.append('  measured TiN/HfO2 effective WF ~ 4.7 eV')
    report.append('  MIGS estimate in paper          ~ 4.6 eV')
    report.append('  raw GGA interface values can be substantially higher because VBO is underestimated')
    report.append('  do NOT force a raw PBE result to equal 4.7 eV; interface chemistry and band-gap error matter.')

    txt='\n'.join(report)+'\n'
    (OUT/'effective_work_function_report.txt').write_text(txt)
    print(txt)

if __name__=='__main__': main()
