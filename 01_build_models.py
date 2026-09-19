#!/usr/bin/env python3
"""Build a literature-anchored TiN(111)/m-HfO2(001) QE benchmark.

This is NOT an atom-for-atom reproduction of Fonseca & Knizhnik (PRB 74,
195304 (2006)); the paper does not publish machine-readable interface
coordinates. It follows their geometry scale/orientations:
  * TiN(111), ~10 A thick
  * monoclinic HfO2(001), ~10 A thick
  * interface cell ~15.78 x 5.29 A^2
  * double-interface periodic cell (no vacuum)

The default model uses Ti-terminated TiN on both sides and an O-terminated
sharp HfO2 slab. The two interfaces are allowed to relax while interior atoms
are fixed by default to retain bulk-like reference regions.
"""
from __future__ import annotations
import json, math, itertools
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent

# ---------- Literature-anchored geometry ----------
# Rocksalt TiN. 4.279 A reproduces the paper's quoted unstrained TiN(111)
# rectangular surface cell (~15.13 x 5.24 A for 5x1).
A_TIN = 4.279

# Interface lateral cell quoted for the TiN/HfO2 matching in the PRB paper.
LX = 15.78
LY = 5.29

# Experimental monoclinic HfO2 (P21/c) structural data; close to GGA values.
A_HFO2 = 5.1145
B_HFO2 = 5.1682
C_HFO2 = 5.2900
BETA_HFO2_DEG = 99.211
HFO2_SITES = {
    "Hf": (0.2751, 0.0396, 0.2073),
    "O1": (0.0660, 0.3260, 0.3390),
    "O2": (0.4560, 0.7570, 0.4870),
}

# 9 alternating (111) atomic planes -> ~9.8 A, Ti termination on both sides.
N_TIN_PLANES = 9
# 3 x 1 in-plane m-HfO2 cells and 2 cells along [001] -> ~10.4 A.
HFO2_REP_X = 3
HFO2_REP_Z = 2
# Initial vertical Ti-O separation at both periodic interfaces.
INTERFACE_GAP = 2.05
# Registry chosen by a small grid search so both interfaces begin with
# closest Ti-O distances ~2.15 A, avoiding overlaps.
HFO2_SHIFT_X = 0.9740740741
HFO2_SHIFT_Y = 0.9031707317

# ---------- DFT settings ----------
# PSlibrary PBE ultrasoft pseudopotentials below recommend up to ~52/575 Ry.
ECUTWFC = 60.0
ECUTRHO = 600.0
DEGAUSS = 0.01

PSEUDOS = {
    "Ti": "Ti.pbe-spn-rrkjus_psl.1.0.0.UPF",
    "N":  "N.pbe-n-rrkjus_psl.1.0.0.UPF",
    "Hf": "Hf.pbe-spn-rrkjus_psl.1.0.0.UPF",
    "O":  "O.pbe-n-rrkjus_psl.1.0.0.UPF",
}

MASS = {"Ti":47.867, "N":14.0067, "Hf":178.49, "O":15.999}


def p21c_equiv(frac):
    """General positions of P21/c (#14), unique axis b, origin choice 1."""
    x, y, z = frac
    ops = [
        (x, y, z),
        (-x, y + 0.5, -z + 0.5),
        (-x, -y, -z),
        (x, -y + 0.5, z + 0.5),
    ]
    return [tuple(float(v % 1.0) for v in p) for p in ops]


def build_tin111():
    """Return 5x1 rectangular TiN(111) slab, stretched to LX x LY."""
    a = A_TIN
    ex = np.array([1.0, -1.0, 0.0]) / math.sqrt(2.0)
    ey = np.array([1.0, 1.0, -2.0]) / math.sqrt(6.0)
    ez = np.array([1.0, 1.0, 1.0]) / math.sqrt(3.0)
    lx_nat = 5.0 * a / math.sqrt(2.0)
    ly_nat = a * math.sqrt(3.0 / 2.0)

    ti_basis = np.array([[0,0,0],[0,.5,.5],[.5,0,.5],[.5,.5,0]], float)
    n_basis  = np.array([[.5,0,0],[0,.5,0],[0,0,.5],[.5,.5,.5]], float)

    unique = {}
    for sp, basis in (("Ti", ti_basis), ("N", n_basis)):
        for t in itertools.product(range(-6, 7), repeat=3):
            t = np.array(t, float)
            for f in basis:
                r = a * (f + t)
                x, y, z = float(r @ ex), float(r @ ey), float(r @ ez)
                fx = round((x / lx_nat) % 1.0, 8) % 1.0
                fy = round((y / ly_nat) % 1.0, 8) % 1.0
                key = (sp, fx, fy, round(z, 6))
                unique[key] = (sp, fx * LX, fy * LY, z)

    zlevels = sorted({round(v[3], 6) for v in unique.values() if v[3] >= -1e-5})
    chosen = set(zlevels[:N_TIN_PLANES])
    atoms = [(sp, x, y, float(round(z,6))) for sp,x,y,z in unique.values()
             if round(z,6) in chosen]
    atoms.sort(key=lambda q: (q[3], q[0], q[1], q[2]))
    return atoms, max(z for _,_,_,z in atoms), lx_nat, ly_nat


def build_hfo2_001():
    """Return 3x1x2 m-HfO2(001) slab, stretched laterally to LX x LY."""
    beta = math.radians(BETA_HFO2_DEG)
    cperp = C_HFO2 * math.sin(beta)
    cshift = C_HFO2 * math.cos(beta)
    lx_nat = HFO2_REP_X * A_HFO2
    ly_nat = B_HFO2

    base = []
    for label, fr in HFO2_SITES.items():
        sp = "Hf" if label == "Hf" else "O"
        for x,y,z in p21c_equiv(fr):
            base.append((sp,x,y,z))

    atoms = []
    for ix in range(HFO2_REP_X):
        for iz in range(HFO2_REP_Z):
            for sp,x,y,z in base:
                X = (x + ix) * A_HFO2 + (z + iz) * cshift
                Y = y * B_HFO2
                Z = (z + iz) * cperp
                atoms.append((sp,X,Y,Z))

    zmin = min(z for *_,z in atoms)
    out = []
    for sp,x,y,z in atoms:
        x = ((x % lx_nat) * LX / lx_nat + HFO2_SHIFT_X) % LX
        y = ((y % ly_nat) * LY / ly_nat + HFO2_SHIFT_Y) % LY
        out.append((sp,x,y,z-zmin))
    out.sort(key=lambda q: (q[3], q[0], q[1], q[2]))
    return out, max(z for _,_,_,z in out), cperp


def atomic_species_block(species):
    order = [s for s in ("Ti","N","Hf","O") if s in species]
    lines = ["ATOMIC_SPECIES"]
    for s in order:
        lines.append(f"{s:2s} {MASS[s]:10.5f}  {PSEUDOS[s]}")
    return "\n".join(lines)


def positions_block(atoms, flags=None):
    lines = ["ATOMIC_POSITIONS angstrom"]
    for i,(sp,x,y,z) in enumerate(atoms):
        if flags is None:
            lines.append(f"{sp:2s} {x:14.8f} {y:14.8f} {z:14.8f}")
        else:
            fx,fy,fz = flags[i]
            lines.append(f"{sp:2s} {x:14.8f} {y:14.8f} {z:14.8f}  {fx:d} {fy:d} {fz:d}")
    return "\n".join(lines)


def cell_block(lx,ly,lz):
    return ("CELL_PARAMETERS angstrom\n"
            f"{lx:14.8f} 0.00000000 0.00000000\n"
            f"0.00000000 {ly:14.8f} 0.00000000\n"
            f"0.00000000 0.00000000 {lz:14.8f}")


def pw_header(calc,prefix,outdir,nat,ntyp,metal=True,conv=1e-8,forces=False):
    ctrl_extra = ""
    if forces:
        ctrl_extra = "\n  tprnfor=.true.,\n  tstress=.true.,\n  forc_conv_thr=2.5d-3,"
    sys_occ = ("\n  occupations='smearing',\n  smearing='mv',\n"
               f"  degauss={DEGAUSS:.4f},") if metal else "\n  occupations='fixed',"
    ions = "\n&IONS\n  ion_dynamics='bfgs',\n/" if calc == "relax" else ""
    return f"""&CONTROL
  calculation='{calc}',
  restart_mode='from_scratch',
  prefix='{prefix}',
  pseudo_dir='./pseudo/',
  outdir='{outdir}',
  verbosity='high',{ctrl_extra}
/
&SYSTEM
  ibrav=0,
  nat={nat},
  ntyp={ntyp},
  input_dft='PBE',
  ecutwfc={ECUTWFC:.1f},
  ecutrho={ECUTRHO:.1f},{sys_occ}
/
&ELECTRONS
  conv_thr={conv:.1e},
  electron_maxstep=300,
  mixing_mode='plain',
  mixing_beta=0.20,
  diagonalization='david',
/{ions}
"""


def write_bulk_tin():
    # Periodic strained TiN(111) reference using same lateral strain as interface.
    # One Ti plane + one N plane; repeat length between equivalent Ti planes.
    slab, _, lx_nat, ly_nat = build_tin111()
    # get first Ti plane and following N plane only
    zvals = sorted({round(z,6) for *_,z in slab})[:2]
    atoms = [(sp,x,y,z) for sp,x,y,z in slab if round(z,6) in zvals]
    z0 = min(z for *_,z in atoms)
    atoms = [(sp,x,y,z-z0) for sp,x,y,z in atoms]
    lz = 2.0 * (A_TIN / (2.0*math.sqrt(3.0)))  # Ti-to-next-Ti spacing
    text = pw_header('scf','bulk_tin_strained','./tmp/bulk_tin/',len(atoms),2,metal=True,conv=1e-9)
    text += atomic_species_block({"Ti","N"}) + "\n"
    text += cell_block(LX,LY,lz) + "\n"
    text += positions_block(atoms) + "\n"
    text += "K_POINTS automatic\n2 6 10 0 0 0\n"
    (ROOT/'02_bulk_tin_strained_scf.in').write_text(text)
    return atoms,lz


def write_bulk_hfo2():
    beta = math.radians(BETA_HFO2_DEG)
    avec = (A_HFO2,0.0,0.0)
    bvec = (0.0,B_HFO2,0.0)
    cvec = (C_HFO2*math.cos(beta),0.0,C_HFO2*math.sin(beta))
    atoms=[]
    for label,fr in HFO2_SITES.items():
        sp='Hf' if label=='Hf' else 'O'
        for x,y,z in p21c_equiv(fr):
            atoms.append((sp,x,y,z))
    atoms.sort(key=lambda q:(q[0],q[1],q[2],q[3]))
    text = pw_header('scf','bulk_hfo2','./tmp/bulk_hfo2/',len(atoms),2,metal=False,conv=1e-9)
    text = text.replace("  occupations='fixed',", "  occupations='fixed',\n  nbnd=64,")
    text += atomic_species_block({"Hf","O"}) + "\n"
    text += "CELL_PARAMETERS angstrom\n"
    text += f"{avec[0]:14.8f} {avec[1]:14.8f} {avec[2]:14.8f}\n"
    text += f"{bvec[0]:14.8f} {bvec[1]:14.8f} {bvec[2]:14.8f}\n"
    text += f"{cvec[0]:14.8f} {cvec[1]:14.8f} {cvec[2]:14.8f}\n"
    text += "ATOMIC_POSITIONS crystal\n"
    for sp,x,y,z in atoms:
        text += f"{sp:2s} {x:14.9f} {y:14.9f} {z:14.9f}\n"
    text += "K_POINTS automatic\n6 6 6 0 0 0\n"
    (ROOT/'03_bulk_hfo2_scf.in').write_text(text)
    return atoms


def write_interface():
    tin, tin_span, tin_lx_nat, tin_ly_nat = build_tin111()
    hfo, hfo_span, cperp = build_hfo2_001()
    # Place HfO2 after TiN; PBC supplies the second interface.
    hfo_z0 = tin_span + INTERFACE_GAP
    hfo = [(sp,x,y,z+hfo_z0) for sp,x,y,z in hfo]
    lz = tin_span + hfo_span + 2.0*INTERFACE_GAP
    atoms = tin + hfo

    # Fix bulk-like central portions; interface-adjacent atoms relax.
    flags=[]
    central_tin_ids=[]; central_hfo_ids=[]
    tin_center = tin_span/2
    hfo_center = hfo_z0 + hfo_span/2
    for i,(sp,x,y,z) in enumerate(atoms, start=1):
        is_tin = i <= len(tin)
        if is_tin:
            fixed = abs(z-tin_center) <= 1.35  # roughly central 3 planes
            if fixed: central_tin_ids.append(i)
        else:
            fixed = abs(z-hfo_center) <= 1.6
            if fixed: central_hfo_ids.append(i)
        flags.append((0,0,0) if fixed else (1,1,1))

    species=set(sp for sp,*_ in atoms)
    base = pw_header('relax','tinhfo2','./tmp/interface/',len(atoms),4,metal=True,conv=1e-7,forces=True)
    base += atomic_species_block(species) + "\n" + cell_block(LX,LY,lz) + "\n"
    base += positions_block(atoms,flags) + "\nK_POINTS automatic\n1 2 1 0 0 0\n"
    (ROOT/'05_interface_relax.in').write_text(base)

    # Smoke SCF: same structure, all coordinates fixed, Gamma-ish sparse k sampling.
    smoke = pw_header('scf','tinhfo2_smoke','./tmp/interface_smoke/',len(atoms),4,metal=True,conv=1e-6)
    smoke += atomic_species_block(species) + "\n" + cell_block(LX,LY,lz) + "\n"
    smoke += positions_block(atoms) + "\nK_POINTS automatic\n1 1 1 0 0 0\n"
    (ROOT/'04_interface_smoke.in').write_text(smoke)

    # Final SCF template populated by 06_extract_relaxed.py.
    final = pw_header('scf','tinhfo2','./tmp/interface/',len(atoms),4,metal=True,conv=1e-9)
    final += atomic_species_block(species) + "\n" + cell_block(LX,LY,lz) + "\n"
    final += "__ATOMIC_POSITIONS__\nK_POINTS automatic\n2 6 1 0 0 0\n"
    (ROOT/'06_interface_scf.template.in').write_text(final)

    metadata = {
      "model": "Ti-terminated TiN(111) / O-terminated monoclinic HfO2(001), periodic double interface",
      "nat": len(atoms), "tin_atoms": len(tin), "hfo2_atoms": len(hfo),
      "cell_A": [LX,LY,lz], "tin_span_A": tin_span, "hfo2_span_A": hfo_span,
      "hfo2_start_A": hfo_z0, "interface_gap_A": INTERFACE_GAP,
      "central_tin_atom_ids": central_tin_ids,
      "central_hfo2_atom_ids": central_hfo_ids,
      "potential_regions_A": {
          "tin": [tin_center-1.4, tin_center+1.4],
          "hfo2": [hfo_center-1.7, hfo_center+1.7],
      },
      "natural_tin_surface_A": [tin_lx_nat,tin_ly_nat],
      "tin_inplane_strain_percent": [100*(LX/tin_lx_nat-1),100*(LY/tin_ly_nat-1)],
      "hfo2_unstrained_inplane_A": [HFO2_REP_X*A_HFO2,B_HFO2],
      "hfo2_inplane_strain_percent": [100*(LX/(HFO2_REP_X*A_HFO2)-1),100*(LY/B_HFO2-1)],
      "notes": "Coordinates are a transparent literature-anchored model, not the unpublished interface coordinates of PRB 74 195304."
    }
    (ROOT/'interface_metadata.json').write_text(json.dumps(metadata,indent=2))
    return atoms,lz


def write_pp_inputs():
    jobs = [
      ('bulk_tin_strained','./tmp/bulk_tin/','outputs/bulk_tin_potential.cube','08_pp_bulk_tin.in'),
      ('bulk_hfo2','./tmp/bulk_hfo2/','outputs/bulk_hfo2_potential.cube','09_pp_bulk_hfo2.in'),
      ('tinhfo2','./tmp/interface/','outputs/interface_potential.cube','10_pp_interface.in'),
    ]
    for prefix,outdir,cube,name in jobs:
        txt=f"""&INPUTPP
  prefix='{prefix}',
  outdir='{outdir}',
  plot_num=11,
  filplot='outputs/{prefix}.pot.tmp',
/
&PLOT
  nfile=1,
  filepp(1)='outputs/{prefix}.pot.tmp',
  weight(1)=1.0,
  iflag=3,
  output_format=6,
  fileout='{cube}',
/
"""
        (ROOT/name).write_text(txt)

    proj="""&PROJWFC
  prefix='tinhfo2',
  outdir='./tmp/interface/',
  Emin=-15.0,
  Emax=15.0,
  DeltaE=0.02,
  ngauss=0,
  degauss=0.01,
  filpdos='outputs/interface_pdos',
/
"""
    (ROOT/'07_interface_projwfc.in').write_text(proj)


def main():
    for d in ('pseudo','tmp','outputs','tmp/bulk_tin','tmp/bulk_hfo2','tmp/interface','tmp/interface_smoke'):
        (ROOT/d).mkdir(parents=True,exist_ok=True)
    write_bulk_tin()
    write_bulk_hfo2()
    atoms,lz=write_interface()
    write_pp_inputs()
    print(f"Built TiN/HfO2 model: {len(atoms)} atoms; cell = {LX:.2f} x {LY:.2f} x {lz:.2f} A")
    print("Generated QE inputs 02-10 and interface_metadata.json")

if __name__ == '__main__':
    main()
