# TiN(111) / monoclinic HfO2(001) Quantum ESPRESSO benchmark

This project is the HfO2 benchmark we actually want for a DRAM-style electrode/high-k interface study. It is **not** the earlier bare TiN(100)/vacuum work-function model.

## What it models

A periodic **double-interface** heterostructure:

`... TiN(111) | m-HfO2(001) | TiN(111) | m-HfO2(001) ...`

The model is anchored to Fonseca & Knizhnik, *Phys. Rev. B* **74**, 195304 (2006):

- TiN(111) metal slab about 10 A thick.
- monoclinic HfO2(001) dielectric slab about 10 A thick.
- lateral matching close to the paper's quoted 15.78 x 5.29 A^2 HfO2 interface cell.
- periodic double-interface geometry (no vacuum) for the VBO/effective-work-function calculation.
- PBE/GGA.
- two analysis routes: projected DOS and a van-de-Walle/Martin-like planar-potential lineup.

The paper does **not** publish machine-readable atomic coordinates for its exact interface. Therefore this bundle uses a transparent, reproducible model: **Ti-terminated TiN(111) facing an O-terminated monoclinic HfO2(001) slab**. It is a literature-anchored benchmark, not an atom-for-atom reproduction.

The generated interface has 162 atoms: 90 Ti/N atoms plus 72 Hf/O atoms. The TiN slab has nine alternating (111) atomic planes (~9.88 A, Ti-terminated on both sides); the HfO2 slab is a 3 x 1 x 2 m-HfO2 supercell (~10.31 A). Interior atoms are fixed during the first relaxation so that bulk-like reference regions remain, while interface-adjacent atoms relax.

## Important expected numbers

Do **not** expect a raw PBE interface calculation automatically to equal 4.7 eV.

The paper reports:

- measured TiN/HfO2 effective work function: about **4.7 eV**;
- MIGS estimate: about **4.6 eV**;
- plain GGA interface VBOs were underestimated, so raw GGA-derived effective work functions could be considerably higher;
- the authors explicitly discuss band-gap/VBO corrections and report better agreement after correction/scaling.

The point of this workflow is to calculate the interface VBO and then obtain

`W_eff = Eg(HfO2) + EA(HfO2) - VBO`

using `Eg_exp = 5.7 eV` and `EA_exp = 2.9 eV`, exactly as in the paper's effective-work-function discussion.

## Pseudopotentials and cutoffs

The project uses one internally consistent PSlibrary PBE-USPP family:

- `Ti.pbe-spn-rrkjus_psl.1.0.0.UPF`
- `N.pbe-n-rrkjus_psl.1.0.0.UPF`
- `Hf.pbe-spn-rrkjus_psl.1.0.0.UPF`
- `O.pbe-n-rrkjus_psl.1.0.0.UPF`

The Ti and Hf pseudopotentials include semicore states, matching the paper's emphasis on semicore treatment. The chosen cutoffs are 60 Ry / 600 Ry, above the largest suggested minima in this PP set.

## Run order

From the project directory:

```bash
chmod +x *.sh *.py
./00_download_pseudos.sh
python3 01_build_models.py
./check_setup.sh
```

### Stage 1: bulk references

```bash
NP=4 bash run_1_references.sh
```

This computes a strained TiN(111) bulk reference and bulk monoclinic HfO2, then writes electrostatic-potential cube files.

### Optional smoke test of the 162-atom interface

```bash
NP=4 bash run_2_smoke.sh
```

This is a coarse 1x1x1 SCF on the unrelaxed interface. It only checks that the large model, pseudopotentials, bands, and memory requirements are viable.

### Stage 2: relax the interface

```bash
NP=4 bash run_3_relax.sh
```

The relaxation uses a 1x2x1 k mesh and a force threshold close to the paper's ~0.07 eV/A criterion. `06_extract_relaxed.py` automatically writes the final SCF input from the relaxed coordinates.

### Stage 3: final SCF, PDOS, electrostatic potential, W_eff

```bash
NP=4 bash run_4_final.sh
```

This runs:

1. final 2x6x1 interface SCF;
2. `projwfc.x` for the central HfO2 PDOS;
3. `pp.x` with `plot_num=11` for electrostatic-potential alignment;
4. `11_analyze.py`.

The main result is written to:

```text
outputs/effective_work_function_report.txt
```

The combined projected DOS from the central oxide region is written to:

```text
outputs/central_hfo2_pdos.dat
```

## How the analysis works

### PDOS method

The script sums the projected DOS of atoms in the middle of the HfO2 slab, identifies the approximate HfO2 valence-band edge below the metal Fermi level, and calculates

`VBO = E_F(TiN) - E_VBM(HfO2)`.

It reports several DOS thresholds because the exact edge in a finite broadened PDOS is not mathematically unique. This mirrors the paper's statement that PDOS was the more reliable method for its short dielectric slabs.

### Potential-lineup method

The script extracts `V_bare + V_H` using `pp.x`, reads the Gaussian cube files, averages the potential in bulk-like TiN and HfO2 regions, and combines those lineups with separate bulk references to obtain a WM-like VBO.

For both routes it reports:

- raw VBO;
- raw `W_eff` using experimental HfO2 gap and affinity;
- a clearly labelled **paper-style heuristic** in which the VBO is scaled by `Eg_exp / Eg_DFT`.

The scaled result is not a GW calculation and should not be presented as one.

## Computational cost

This is a much larger calculation than the earlier bare TiN test: 162 atoms and hundreds of occupied bands. A serial run can be very slow. Start with `NP=4` if the workstation has at least four cores; check `nproc` and available memory with `free -h`. If the machine has 8-16 physical/logical cores and enough RAM, `NP=8` may be reasonable.

Do not increase MPI processes blindly if memory becomes a problem.

## Files

- `00_download_pseudos.sh` — downloads the four matching PBE USPP files.
- `01_build_models.py` — builds bulk references and the interface model.
- `02_bulk_tin_strained_scf.in` — TiN reference with the interface in-plane strain.
- `03_bulk_hfo2_scf.in` — bulk monoclinic HfO2 reference.
- `04_interface_smoke.in` — optional cheap interface SCF.
- `05_interface_relax.in` — interface relaxation.
- `06_interface_scf.template.in` — final-SCF template.
- `06_extract_relaxed.py` — extracts relaxed coordinates and creates `06_interface_scf.in`.
- `07_interface_projwfc.in` — PDOS calculation.
- `08_pp_bulk_tin.in`, `09_pp_bulk_hfo2.in`, `10_pp_interface.in` — electrostatic-potential cube generation.
- `11_analyze.py` — VBO and effective-work-function analysis.
- `interface_metadata.json` — atom indices, cell/region definitions, strains.
