# PolyCleaver 2

**Non-polar, stoichiometric surface slabs of complex minerals — polyanion scaffolds,
counter-cations, hydration water and intercalated molecules — from a bulk CIF, ready
for DFT.**

PolyCleaver 2 is the successor to
[PolyCleaver](https://github.com/ericmates/PolyCleaver)
(Mates-Torres & Rimola, *J. Appl. Cryst.* **57**, 503–508, 2024,
[doi:10.1107/S1600576724001286](https://doi.org/10.1107/S1600576724001286)).
Where PolyCleaver handled simple silicates, sulfates, phosphates and carbonates by
trimming cations sequentially, PolyCleaver 2 rebuilds every slab chemically:

1. **Analyse the bulk** and classify each atom as a *scaffold centre* (Si, P, S, Al, …),
   *scaffold oxygen*, *counter-cation* (Mg²⁺, Ca²⁺, K⁺, …) or part of a *discrete
   molecule* (H₂O, OH⁻, NH₃, NH₄⁺ …). Structural hydroxyls stay attached to their
   polyhedra; hydration shells around cations are measured (6 H₂O per Mg²⁺ in struvite).
2. **Cut** with pymatgen's `SlabGenerator`, then **strip** the raw slab to its covalent
   scaffold — keeping every polyhedron intact for isolated-polyhedra minerals, keeping
   the broken surface tetrahedra for chain / sheet / framework silicates.
3. **Thin** the scaffold layer by layer and record every thickness state.
4. **Put the cations back** by a combinatorial search over the candidates closest to
   the slab centre, minimising the z-dipole.
5. **Put the molecules back** the same way — water sorted by distance to the nearest
   cation (PBC-correct), other species by distance to the slab centre.
6. **Accept** only slabs with |dipole<sub>z</sub>| ≤ tolerance (0.3 eÅ by default) and
   zero net charge, **rank** them by broken bonds per Å² (a DFT-validated proxy for the
   unrelaxed surface energy) and write one CIF per unique termination.

The result is the *minimal, non-polar, charge-neutral, stoichiometric* slab for each
Miller index — or, in range mode, **every** such termination between two thicknesses.

> **Status (September 2026):** public development snapshot of the code used in the
> forthcoming PolyCleaver 2 paper. The algorithm is complete and validated on the
> examples below; the API may still change before the first tagged release.

---

## What it handles

| Structure type | Example | How |
|---|---|---|
| Isolated polyhedra (nesosilicates, phosphates, sulfates) | forsterite Mg₂SiO₄, apatite | incomplete surface polyhedra removed; cations cancel the scaffold dipole |
| Chain silicates (pyroxenes) | diopside CaMgSi₂O₆, enstatite | layer-wise thinning; non-stoichiometric cuts rejected by the charge check |
| Sheet and framework silicates | micas, feldspars, zeolites | broken surface tetrahedra kept; scaffold-dipole filter before cation placement (implemented; fewer validated examples than for isolated polyhedra so far) |
| Structural hydroxyl | hydroxyapatite Ca₁₀(PO₄)₆(OH)₂, kaolinite | OH⁻ as a molecular species, or H attached to the Al/Si unit |
| Hydrated / intercalated molecules | struvite MgNH₄PO₄·6H₂O | H₂O and NH₄⁺ placed combinatorially, water near cations |
| Disordered CIFs (Fe/Mg, Al/Si, vacancies) | Al/Si on the tetrahedral sites of muscovite | ordered with pymatgen's Ewald-ranked `OrderDisorderedStructureTransformation`; the CLI warns if the ordered cell is not charge-neutral |
| Charged / non-stoichiometric terminations | reactivity studies | optional `allow_nonstoichiometric` (files tagged `_NS`) |

Not in scope: structures without a polyanion scaffold (simple oxides, sulfides,
metals) and surfaces that are polar by construction, such as the (001) face of 1:1
clays like kaolinite, whose neutral layer carries an intrinsic dipole.

---

## Installation

```bash
git clone https://github.com/ericmates/polycleaver2.git
cd PolyCleaver2
pip install -e .            # pymatgen, ase, numpy
pip install -e ".[webapp]"  # optional: FastAPI backend for the web interface
pip install -e ".[test]"    # optional: pytest
```

Python ≥ 3.10. Tested with pymatgen 2026.5 and ase 3.28.

---

## Quick start

```bash
polycleaver2                                          # reads ./config.json
polycleaver2 examples/configs/forsterite_110_range.json -o out/
python -m polycleaver2 my_run.json                    # same thing
```

Everything is set in a JSON file (keys beginning with `_` are comments):

| Key | Default | Meaning |
|---|---|---|
| `filename` | — | bulk structure (CIF, POSCAR, … anything pymatgen reads) |
| `target_thickness` | 15.0 | minimum slab thickness (Å); the thinnest valid state ≥ target is kept |
| `thickness_range` | `null` | `[min, max]` Å: return **all** valid terminations in the window (overrides `target_thickness`) |
| `vacuum_size` | 15.0 | vacuum on each side (Å) |
| `dipole_threshold` | 0.3 | max \|dipole<sub>z</sub>\| (eÅ) |
| `auto_select_stable` | `true` | rank by broken bonds / Å², then symmetry, then coordination saturation |
| `allow_nonstoichiometric` | `false` | keep charged slabs as `slab_hkl_N_NS.cif` |
| `mode` | `"single"` | `"single"` (`single_hkl`), `"max"` (all distinct indices ≤ `max_index`), `"list"` (`custom_list`) |

Output: `slab_{hkl}_{i}.cif` per accepted termination plus a terminal report with
formula, thickness, symmetry, residual dipole, surface cation hydration and the
stability rank.

### From Python

```python
from pymatgen.core import Structure
from polycleaver2 import ChemicalAnalyzer, ensure_oxidation_states, generate_slabs_for_hkl

bulk = ensure_oxidation_states(Structure.from_file("examples/structures/hydroxyapatite.cif"))
an = ChemicalAnalyzer(bulk)
an.print_detailed()                                   # scaffold / cations / molecules / hydration

slabs = generate_slabs_for_hkl(bulk, an, (0, 1, 0),
                               target_thickness=17.0, vacuum=15.0, dipole_tol=0.3)
for i, slab in enumerate(slabs):
    print(slab.formula, slab.is_symmetric(), float(slab.dipole[2]))
    slab.to(filename=f"HAP_010_{i}.cif")
```

---

## Examples and validation

All configs in `examples/configs/` run in seconds to a minute on a laptop.

| Config | System | Result |
|---|---|---|
| `forsterite_110_range.json` | Mg₂SiO₄ (110), 12–20 Å | 8 unique terminations; #1 is the symmetric Mg₃₂Si₁₆O₆₄ slab at 14.99 Å, dipole 0.000 eÅ — identical to the DFT-relaxed reference `examples/reference/POSCAR_forsterite_110_DFT` |
| `forsterite_010.json` | Mg₂SiO₄ (010), ≥ 10 Å | 2 symmetric terminations (11.9 and 14.9 Å). The bare SiO₄ scaffold is polar here; Mg²⁺ placement cancels it |
| `hydroxyapatite_010.json` | Ca₁₀(PO₄)₆(OH)₂ (010), ≥ 17 Å | Ca₂₀(PO₄)₁₂(OH)₄, 88 atoms, symmetric, dipole 0.000 eÅ — matches `POSCAR_hydroxyapatite_010_DFT`; the four cuts with n<sub>P</sub> ∤ 3 are correctly rejected as charged |
| `struvite_low_index.json` | MgNH₄PO₄·6H₂O, (001)/(010)/(100) | 7 terminations, all symmetric with dipole ≤ 0.12 eÅ and complete Mg(H₂O)₆ shells |
| `diopside_max1.json` | CaMgSi₂O₆, all indices ≤ 1 | 7 symmetric terminations on 5 planes; polar or odd-chain cuts discarded |

Stability ranking was checked against PBE surface energies of three forsterite (110)
terminations: broken bonds / Å² reproduces the DFT order (1.74 < 1.81 < 1.90 J m⁻²),
whereas the per-atom coordination-saturation average runs the wrong way. See
`docs/DESIGN_NOTES.md`.

Run the test-suite with `pytest tests/` (add `-m "not slow"` to skip the full
generation runs).

---

## Documentation

- [`docs/HOW_IT_WORKS.md`](docs/HOW_IT_WORKS.md) — plain-language walk-through of the pipeline and troubleshooting
- [`docs/ALGORITHM_FLOWCHART.md`](docs/ALGORITHM_FLOWCHART.md) — full Mermaid flowchart
- [`docs/DESIGN_NOTES.md`](docs/DESIGN_NOTES.md) — the non-obvious decisions: framework vs isolated polyhedra, why C(req+4, 4), the charge check, stability metric, DFT validation, limitations
- [`docs/WEBAPP.md`](docs/WEBAPP.md) — the Surface Lab web interface

## Web interface (Surface Lab)

`api/` (FastAPI) and `frontend/` (React + Three.js) provide an interactive front end:
upload a CIF, inspect the chemistry, generate and view slabs with coordination
polyhedra, download CIFs, and optionally relax bulk and slabs with MACE-MP to get
surface energies and a Wulff shape.

```bash
uvicorn api.main:app --reload --port 8000     # backend, from the repo root
cd frontend && npm install && npm run dev     # frontend on http://localhost:5173
```

![Surface Lab](docs/webapp_screenshot.png)

---

## Related tools

- **PolyCleaver** (v1): same goal for simple polyatomic-anion compounds; no molecular species, sequential cation trimming.
- **SALAMI** (Xie *et al.*, [arXiv:2608.28903](https://arxiv.org/abs/2608.28903), 2026): symmetric stoichiometric slabs of multinary compounds by combinatorial removal of surface atoms under coordination constraints, with Ewald/potential pre-screening. PolyCleaver 2 instead rebuilds the slab from a chemically classified scaffold, which is what makes hydrated and intercalated species tractable.
- **pymatgen** `get_tasker2_slabs`, **surfaxe**: layer- or ion-shuffling approaches without chemical classification.

## Citing

Until the PolyCleaver 2 paper is out, please cite the software (see `CITATION.cff`) together with the original PolyCleaver paper:

> E. Mates-Torres, A. Rimola, *J. Appl. Cryst.* **57**, 503–508 (2024). doi:10.1107/S1600576724001286

## License

GPL-3.0, as PolyCleaver. Example structures come from the Crystallography Open Database
(struvite COD 9007674, kaolinite COD 9009234) and from published DFT work; their headers
carry the original attribution.
