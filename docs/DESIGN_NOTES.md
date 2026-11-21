# Design notes and lessons learned

Developer-facing notes on the non-obvious decisions in `polycleaver2/slabgen.py`.
Read `HOW_IT_WORKS.md` first for the plain-language pipeline.

## Atom classification (`ChemicalAnalyzer`)

| Category    | Rule |
|-------------|------|
| `centers`   | Covalent element (see `COVALENT_ELEMENTS`), positive oxidation state, ≥ 1 O neighbour within 2.3 Å. Si, P, S, Al, B, C, N, As … |
| `anions`    | O atoms in a covalent cluster that contains at least one centre |
| `cations`   | Every element **not** in `COVALENT_ELEMENTS` (Mg, Ca, K, Na, Fe, …) |
| `molecules` | Covalent clusters with **no** centre atom: H₂O, OH⁻, NH₃, NH₄⁺ … keyed by a sorted formula string (`H2O`, `HO`, `H4N`) |

Bond thresholds: X–H < 1.25 Å, other covalent bonds < 2.3 Å; neighbour search radius 2.8 Å.
Elements outside `OXI_MAP` that pymatgen cannot guess raise an actionable error.
Transition metals that form covalent polyhedra (e.g. Nb⁵⁺ in belkovite) are treated as
*cations* because they are not in `COVALENT_ELEMENTS`; their polyhedra are therefore not
protected during cleavage. Extend `COVALENT_ELEMENTS`/`OXI_MAP` to change that.

`cation_hydration_target` is the modal number of molecular O atoms within 2.8 Å of
each cation element in the bulk (6 for Mg²⁺ in struvite). It is reported in the bulk
summary; water placement itself is driven by the O–cation distance sort in Phase 2.

## Framework vs isolated polyhedra

`has_shared_oxygens()` — any O bonded to ≥ 2 centres — decides the code path:

| | Isolated polyhedra (olivine, apatite, gypsum, struvite) | Framework / chain / sheet (feldspar, pyroxene, mica, zeolite) |
|---|---|---|
| Stripping | Remove cations, molecules **and** any surface polyhedron with CN < bulk CN, then orphaned O | Remove cations and molecules only; broken surface tetrahedra are kept (the cleavage plane must cut them) |
| Thinning unit | One outermost centre per side | One whole crystallographic z-layer per side (`_get_center_layers`, gap threshold `max(1.5, 2.5 × median gap)`) |
| Scaffold-dipole filter | **None** — the bare scaffold is always polar; Mg²⁺/Ca²⁺ cancel it | Scaffold-only \|dipole_z\| ≤ tolerance |

**Why no scaffold filter for isolated polyhedra.** The first version filtered on the
scaffold dipole for every mineral. For forsterite (010) the SiO₄ sub-lattice never
reaches 0.3 eÅ at any thickness, so no state was ever recorded and Phase 1 never ran:
the surface "did not exist" unless the tolerance was raised to ~3 eÅ. Collecting all
thickness states and filtering only the finished slab fixed it.

## Combinatorial placement (Phases 1 and 2)

Candidates are the cations / molecules captured from the *full* raw cut before
stripping. They are sorted by proximity to the scaffold centre of mass (cations, and
molecules other than water) or by PBC-correct O–cation distance (water, via
`ase.geometry.get_distances`). The `req + 4` nearest form the pool and all
C(req + 4, req) = C(req + 4, 4) subsets are tried, exiting as soon as
\|dipole_z\| ≤ tolerance and always keeping the best subset seen.

The original pool (`req + 2`, with the innermost `req − 2` fixed) searched only 15
subsets for hydroxyapatite (010) and missed the Ca(I)/Ca(II) arrangement that cancels
the dipole; every cut failed with residual dipoles of 7–13 eÅ. `req + 4` fixed it.
Cost is O(req⁴) but early exit makes it milliseconds in practice
(HAP req = 17 → 5985 subsets; forsterite req = 44 → 194 580). A stagnation guard
(`_STAGNATION_PATIENCE`) aborts hopeless searches early.

Phase 2 originally placed molecules by distance only, which undid the Phase 1 balance;
it now runs the same combinatorial search per molecular species, with the distance
sort as a pre-sort.

## Charge check and rounding

`req = round(n_centres × bulk_ratio)` for each species. When a scaffold state has a
number of centres incompatible with the formula unit (HAP with n_P not divisible by 3,
odd numbers of chain tetrahedra in pyroxenes) the rounded counts leave a net charge.
Such slabs are discarded by the final `Σ oxi_state ≠ 0` check, or kept with an `_NS`
suffix when `allow_nonstoichiometric` is set. If the **bulk** itself is not neutral
with the assigned oxidation states (e.g. a hydroxyl-free mica model) the CLI warns
that no stoichiometric slab can pass.

## Minimum thickness is a crystallographic quantity

`target_thickness` is a lower bound and the thinnest valid state ≥ target is chosen.
For minerals with a long repeat along the surface normal the smallest neutral,
non-polar slab may need two formula units: HAP (010) needs ≥ 17 Å (b ≈ 9.4 Å),
forsterite (010) needs ≥ 10 Å (one b period). If a run yields nothing, raise the
target or use `thickness_range` to see every state.

## Stability ranking

`score_surface_stability` counts, for all non-H atoms within 3.5 Å of either surface,
`CN_bulk − CN_surf` with element-pair-specific rules (`_cn_of_site`):

| Atom | Counted neighbours |
|------|--------------------|
| O | non-O, non-H within 2.8 Å |
| Scaffold centre | O within **2.3 Å** (a 2.8 Å cutoff picks up ~2.5 Å non-bonded O and gives Si CN = 7) |
| Cation | O within 2.8 Å |

Bulk reference CNs for forsterite: Mg 6.0, Si 4.0, O 4.0.

Primary key = Σ deficit / area (`broken_per_area`, ∝ unrelaxed surface energy);
tiebreakers: symmetric first, then higher coordination saturation.

DFT validation, forsterite (110), PBE:

| Thickness (Å) | coord. saturation | broken / Å² | γ (J m⁻²) |
|---------------|-------------------|-------------|-----------|
| 12.64 | 0.653 | lowest | 1.74 |
| 13.85 | 0.719 | mid | 1.81 |
| 15.06 | 0.764 | highest | 1.90 |

Broken bonds / Å² orders these correctly; the per-atom *coordination saturation*
average is anti-correlated with γ and must not be used as the primary key.
The metric weights every bond equally (Si–O ≈ 460 kJ mol⁻¹ vs Mg–O ≈ 200 kJ mol⁻¹);
a bond-energy-weighted variant is future work. Surfaces terminated by intact molecules
(struvite) score 0 because no scaffold or ionic bond is cut.

`_deduplicate_slabs` merges slabs with the same formula and thickness (0.1 Å) after
sorting, so the best-ranked duplicate survives.

## Disordered input

CIFs with fractional occupancies (Fe/Mg solid solution, Al/Si on tetrahedral sites,
vacancy-ordered oxides) are ordered with pymatgen's
`OrderDisorderedStructureTransformation` (Ewald-ranked) before analysis. For phases
whose ordering pattern defines the symmetry (maghemite) supply an explicitly ordered
CIF instead.

## Known limitations

- A polyanion scaffold must exist; simple oxides, sulfides and metals are out of scope.
- Pool `req + 4` can miss the right cations in very large cells (req > ~50 is slow and
  the correct set may lie outside the nearest `req + 4`).
- Framework thinning removes whole z-layers; for chain and some sheet silicates this can
  leave non-stoichiometric centre counts on many cuts, which are then discarded.
- Fixed 10 Å generation buffer above the requested thickness.
- Stability metric is not bond-energy weighted.
