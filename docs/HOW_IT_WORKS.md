# How the Slab Generator Works

A plain-language guide for materials scientists who want to understand
what `polycleaver` does without reading every line of code.

---

## What problem does it solve?

When you cut a crystal along a crystallographic plane (a Miller index like (010)),
the exposed surface usually has unbalanced charge — positive charges on one side,
negative on the other. This creates a macroscopic electric dipole along the surface
normal (z-direction), which makes DFT calculations diverge.

This script finds surface cuts that are **non-polar**: the charges are distributed
so that the z-dipole moment is essentially zero.

---

## Quick start

1. Install the package: `pip install polycleaver` (or `pip install -e .` from the repository root).
2. Edit `config.json` — set `filename`, `mode`, and the Miller index you want —
   or give the same parameters as command-line options.
3. Run:

```bash
polycleaver                                  # reads ./config.json
polycleaver my_experiment.json               # a different config file
polycleaver my_experiment.json -o out/       # write CIFs somewhere else
polycleaver -s bulk.cif --hkl 0 1 0 --thickness 10 -o out/   # no config file
```

Output CIFs are written to the output directory: `slab_010_0.cif`, `slab_010_1.cif`, etc.,
together with `polycleaver_summary.json`, which lists the parameters, the library
versions and the properties of every slab. Ready-made configs for the bundled
examples live in `examples/configs/`.

---

## Key concepts

### The polyanion scaffold

Most complex minerals are built around a **scaffold**: a network of polyanion groups
held together by covalent bonds. Examples:

| Mineral class | Scaffold unit | Central atom |
|---------------|---------------|--------------|
| Phosphates (apatite) | PO₄³⁻ | P⁵⁺ |
| Sulfates (gypsum) | SO₄²⁻ | S⁶⁺ |
| Silicates (olivine) | SiO₄⁴⁻ | Si⁴⁺ |
| Framework silicates (zeolites) | –Si–O–Si– network | Si⁴⁺, Al³⁺ |

The script automatically identifies these from the bulk CIF by looking for
covalent-element atoms (P, S, Si, Al…) bonded to oxygen.

### Counter-cations

The scaffold's negative charge is balanced by **ionic cations** — atoms like
Ca²⁺, Mg²⁺, K⁺, Na⁺ that sit between scaffold units and are not covalently bonded.

### Molecular species

Some minerals also contain **discrete molecules** — water (H₂O), hydroxide (OH⁻),
ammonia (NH₃) — identified as connected clusters of atoms that contain no scaffold
centre atoms.

### The dipole moment

The dipole moment is a measure of charge separation. For a surface slab, we care
about the **z-component**: if it is large, the slab has a macroscopic electric
field through it that makes DFT calculations unreliable. The script only keeps
slabs where |dipole_z| ≤ `dipole_threshold` (default 0.3 eÅ).

---

## The generation pipeline, step by step

### Step 1 — Read and analyse the bulk

The script reads your CIF file and assigns oxidation states (Si⁴⁺, O²⁻, Ca²⁺, etc.)
using pymatgen's automatic guesser, with a fallback to a hardcoded table for
common minerals.

It then classifies every atom into one of four roles:
- **Scaffold centres** — the P, S, Si… atoms at the heart of each polyanion
- **Scaffold oxygens** — the O atoms covalently bonded to those centres
- **Counter-cations** — ionic species (Ca, Mg, K…)
- **Molecules** — discrete clusters like H₂O or OH⁻

### Step 2 — Generate raw surface cuts

For the chosen Miller index (e.g. (010)), pymatgen generates all
crystallographically distinct ways to cut the crystal (typically 5–15 "cuts").
Each cut exposes a different atomic termination. The script adds a 10 Å buffer
to the target thickness at this stage to ensure enough material is available
for trimming.

### Step 3 — Strip to bare scaffold

Each cut is stripped down to just the covalent scaffold (P + bonded O atoms,
or Si–O–Si network). All cations, water, and OH⁻ are removed.

For **isolated polyhedra** (phosphates, sulfates, orthosilicates):
incomplete polyhedral units at the surface — those missing one or more oxygens
compared to the bulk coordination number — are also removed. This ensures
only fully-coordinated scaffold groups remain.

For **framework silicates** (zeolites, feldspars, clays): broken surface
tetrahedra are kept, because the cleavage plane inherently exposes them.

### Step 4 — Thin the scaffold to the target thickness

The algorithm removes the topmost and bottommost scaffold units in pairs
(isolated polyhedra: one unit per side; framework/sheet silicates: one whole
crystallographic z-layer per side), one iteration at a time, until nothing is
left. Every intermediate state is recorded together with its thickness.

- **Single-target mode** (`target_thickness`): the *thinnest* recorded state
  that is still ≥ `target_thickness` is kept. If none is thick enough, the
  thickest valid state is used instead (the irreducible minimum period).
- **Range mode** (`thickness_range: [min, max]`): *every* recorded state whose
  thickness falls in the window is kept, and steps 5–7 run for each of them,
  so one Miller index can yield several distinct terminations.

Which states count as "valid" depends on the scaffold type:

- **Framework / sheet silicates** (M–O–M bridges): the bare scaffold is
  stoichiometrically self-contained, so a state is only valid if the
  scaffold's own |dipole_z| ≤ `dipole_threshold`.
- **Isolated polyhedra** (orthosilicates, phosphates, sulfates …): the bare
  SiO₄/PO₄ scaffold is *always* charged and polar — the counter-cations are
  what cancels the dipole. All states are therefore kept and the dipole
  filter is applied only to the finished slab after steps 5–7. Without this
  distinction forsterite (010) would never be found at 0.3 eÅ.

> **Important for HAP (hydroxyapatite) and other minerals with large unit cells:**
> If `target_thickness` is set too low, the algorithm thins past the natural
> crystallographic minimum needed for a balanced slab. For HAP (010), use
> `target_thickness = 17.0` (not 15.0). See [Troubleshooting](#troubleshooting).

### Step 5 — Restore cations (Phase 1)

The stoichiometric ratio of cations to scaffold centres is computed from the
bulk (e.g. Ca/P = 10/6 = 1.667 for HAP). Given the number of scaffold centres
in the thinned slab, the required number of cations is calculated.

From the full pool of available cation candidates (positions extracted from
the original raw cut), the script selects the `req + 4` candidates closest
to the slab's centre of mass. It then searches **all possible combinations**
of `req` cations from this pool — C(req+4, req) = C(req+4, 4) combinations
— evaluating the z-dipole for each. The combination with the lowest |dipole_z|
is kept; the search stops early as soon as one satisfies the threshold.

**Why C(req+4, 4) and not just take the closest `req`?**
Many minerals (especially HAP, with its two distinct Ca sites) require a
specific symmetry of cation placement that the simple "take the nearest"
heuristic misses. Searching all combinations within a small pool guarantees
the correct solution is found.

### Step 6 — Restore molecules (Phase 2)

The same combinatorial approach is applied to each molecular species (OH⁻,
H₂O, etc.). Candidates are pre-sorted so that the most physically reasonable
placements appear first — water/hydroxyl units are sorted by their distance
to the nearest cation; others by distance to the slab centre. This pre-sort
makes the early-exit fast when a valid solution exists.

### Step 7 — Final dipole and charge check

After all cations and molecules are placed, the total |dipole_z| of the
complete slab is evaluated one final time. Slabs above the threshold are
discarded, and so are slabs whose summed oxidation states are not zero
(unless `allow_nonstoichiometric` is set, in which case they are written with
an `_NS` suffix for reactivity studies where you add H atoms by hand).

### Step 8 — Stability ranking and de-duplication

Every accepted slab is scored with a broken-bond model: for each non-H atom
within 3.5 Å of either surface, the coordination number is compared with the
bulk value for that element, and the total deficit is divided by the surface
area. **Lower broken bonds / Å² = more stable** (it is proportional to the
unrelaxed surface energy). Slabs are sorted by this metric, then symmetric
ones first, and slabs with identical formula and thickness are merged. For
forsterite (110) the ranking reproduces the DFT surface-energy order of three
terminations; note that the simpler "coordination saturation" average goes
the *wrong* way and is only reported as a secondary number.

---

## The config.json reference

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `filename` | string | (required) | Path to the bulk structure file (CIF, POSCAR, …) |
| `target_thickness` | float | 15.0 | Minimum slab thickness in Å (single-target mode) |
| `thickness_range` | [min,max] | `null` | Enumerate *all* terminations between min and max Å; overrides `target_thickness` |
| `vacuum_size` | float | 15.0 | Vacuum gap in Å |
| `dipole_threshold` | float | 0.3 | Max allowed \|dipole_z\| in eÅ |
| `auto_select_stable` | bool | `true` | Rank accepted slabs by broken bonds / Å² (see below) |
| `allow_nonstoichiometric` | bool | `false` | Keep charged slabs (written as `slab_hkl_N_NS.cif`) |
| `mode` | string | `"single"` | `"single"` \| `"max"` \| `"list"` |
| `single_hkl` | [h,k,l] | [0,0,1] | Used when `mode = "single"` |
| `max_index` | int | 2 | Used when `mode = "max"` |
| `custom_list` | [[h,k,l]…] | — | Used when `mode = "list"` |

Keys beginning with `_` are treated as documentation comments and are ignored.

---

## Output files

For each accepted slab the script:
- Prints a summary table with formula, thickness, symmetry, dipole, and composition
- Writes a CIF file named `slab_{hkl}_{index}.cif`, e.g. `slab_010_0.cif`

A **symmetric** slab (pymatgen `is_symmetric() = True`) means the top and bottom
surfaces are related by inversion symmetry. This is desirable for DFT because
both surfaces are equivalent — a useful sanity check.

---

## Troubleshooting

### No surfaces found

The diagnostic table is printed when all cuts fail. Common fixes:

1. **Raise `dipole_threshold`** to 0.5 or 1.0 — some minerals cannot
   achieve 0.3 eÅ but are still usable.
2. **Raise `target_thickness`** — see the note below.
3. **Try different Miller indices** — some cuts are inherently polar
   regardless of cation placement.

### Wrong thickness / slabs too thin (HAP and large unit cells)

`step_5b_thin_to_target` keeps the thinnest scaffold state that is still
≥ `target_thickness`. If the target is set too low, that state has fewer
layers than the crystal chemistry requires for a balanced slab.

**Rule of thumb:** `target_thickness` must be at least as large as the
minimum crystallographic repeat along the surface normal that produces a
charge-balanced, non-polar slab.

| Mineral | Surface | Recommended `target_thickness` |
|---------|---------|-------------------------------|
| Hydroxyapatite (HAP) | (010) | **17.0 Å** (not 15.0) |
| Most silicates, sulfates, phosphates | any | 15.0 Å |

For HAP (010): the b-axis period is ≈ 9.4 Å. A valid slab needs at least
2 formula units ≈ 18.8 Å of crystal. Using target=15.0 trims past this,
leaving OH⁻ candidates that cannot balance the dipole in any combination.

### Slow generation

The combinatorial search grows as O(req⁴) but exits early when a solution
is found. If generation is slow, the mineral likely has `req_cats > 50`
(check the "Phase-1 cations: req=…" line in the output). In that case,
raising `target_thickness` to produce a thicker but simpler cut may help.
