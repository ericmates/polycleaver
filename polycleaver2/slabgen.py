#!/usr/bin/env python3
"""
PolyCleaver 2 — non-polar, stoichiometric surface slab generator for complex
minerals (silicates, phosphates, sulfates, clays, hydrated and intercalated
structures).

This module is the computational core: bulk chemical analysis, slab stripping
and thinning, combinatorial cation/molecule placement, and stability scoring.
The command-line entry point lives in ``polycleaver2.cli``.

Given a bulk CIF structure, this script:
  1. Analyses the bulk to identify the polyanion scaffold (centers + O),
     counter-cations, and molecular species (H2O, NH3, OH⁻, etc.).
  2. Generates surface slabs using pymatgen's SlabGenerator.
  3. Strips away cations and molecules, leaving only the covalent scaffold.
  4. Thins the scaffold to the target thickness while minimising the
     z-dipole moment.
  5. Re-adds cations in the stoichiometric ratio, searching combinatorially
     for the configuration that yields the smallest z-dipole.
  6. Re-adds molecular species in the correct ratio, placing water molecules
     near surface cations and other species near the slab centre.
  7. Filters by a user-set dipole tolerance and writes the accepted slabs
     to CIF files.

Supports:
  - Isolated polyanion units (orthosilicates, phosphates, sulfates …)
  - 2-D framework / sheet silicates (phyllosilicates, e.g. muscovite,
    which are inherently polymerised but whose cleavage planes never expose
    under-coordinated metals).
  - 3-D framework silicates (zeolites, feldspars …)
  - Structures with intercalated molecules (water, ammonia, OH⁻ …)
  - Solvated / hydrated cations (e.g. Mg²⁺ in struvite).

Dependencies:  pymatgen, ase, numpy
"""

import sys
import copy
import json
import math
import argparse
import traceback
import numpy as np
import itertools
from collections import Counter, defaultdict
from pathlib import Path

from pymatgen.core import Structure, Species
from pymatgen.core.surface import (
    Slab,
    SlabGenerator,
    get_symmetrically_distinct_miller_indices,
)
from pymatgen.io.ase import AseAtomsAdaptor
from ase.geometry import get_distances


# =============================================================================
#  Utility helpers
# =============================================================================

# Rough covalent-element pool. Non-metal(oid) elements that commonly form
# covalent bonds in minerals.  Metals (K, Mg, Ca, Na, Fe …) are assumed
# to be ionic counter-cations.
COVALENT_ELEMENTS = frozenset({
    'H', 'B', 'C', 'N', 'O', 'F',
    'Al', 'Si', 'P', 'S', 'Cl',
    'Ge', 'Ga', 'As', 'Se', 'Br', 'I',
})

# Known (element → oxidation state) pairs for minerals where pymatgen's
# automatic guess fails (e.g. P1 structures without symmetry).
OXI_MAP = {
    'Si': 4, 'Al': 3, 'K': 1, 'O': -2,
    'Mg': 2, 'Ca': 2, 'Na': 1, 'Fe': 3,
    'P': 5, 'S': 6, 'Ti': 4, 'Mn': 2,
    'Zn': 2, 'Cu': 2, 'Ba': 2, 'Li': 1,
    'N': -3, 'Cl': -1, 'F': -1, 'H': 1,
}

# Maps id(slab) → net_charge for slabs that passed with allow_nonstoichiometric=True.
# Populated by _run_phases_on_scaffold; consumed by report_slab_detail and main.
_nonstoich_charges: dict = {}


def order_if_disordered(structure, max_structures=1):
    """
    Return an *ordered* copy of *structure*.

    Mineral CIFs very often carry fractional site occupancies — either genuine
    solid solution (Fe/Mg on one site) or, more importantly here, *vacancy
    ordering* (γ-Fe₂O₃ maghemite, where the chirality IS the vacancy pattern).
    A disordered ``Structure`` breaks everything downstream: ``site.specie``
    raises ``AttributeError`` on a mixed-occupancy site, and charge-balance
    guessing raises ``ValueError`` on fractional compositions.

    Returns the structure unchanged when it is already ordered.
    Raises ValueError with an actionable message when ordering is impossible.
    """
    if structure.is_ordered:
        return structure

    from pymatgen.transformations.standard_transformations import (
        OrderDisorderedStructureTransformation,
    )

    s = structure.copy()
    # Ewald ranking inside the transformation needs oxidation states; assign
    # them from OXI_MAP, which works on disordered sites too (unlike guessing).
    # Every element must get a value or Ewald ranking raises on the first
    # bare Element it meets, so fill gaps in OXI_MAP with a neutral default.
    _DEFAULT_OXI = {'H': 1, 'O': -2, 'C': 4, 'N': -3, 'S': 6, 'P': 5}
    try:
        s.add_oxidation_state_by_element({
            el.symbol: OXI_MAP.get(el.symbol, _DEFAULT_OXI.get(el.symbol, 0))
            for el in s.composition.elements
        })
    except Exception:
        pass

    try:
        trans = OrderDisorderedStructureTransformation(algo=0)
        out = trans.apply_transformation(s, return_ranked_list=max_structures)
        ordered = out[0]["structure"] if isinstance(out, list) else out
    except Exception as exc:
        raise ValueError(
            f"Structure has fractional site occupancies and could not be "
            f"ordered automatically ({type(exc).__name__}: {exc}). "
            f"Supply an explicitly ordered CIF — for vacancy-ordered phases "
            f"such as maghemite the ordering pattern determines the symmetry "
            f"(and the chirality), so it must be chosen deliberately, not guessed."
        ) from exc

    return ordered.remove_oxidation_states() or ordered


def ensure_oxidation_states(structure):
    """
    Return a *copy* of *structure* with oxidation states assigned.
    Falls back to a hard-coded map when pymatgen's guess fails on all sites
    or assigns obviously wrong values (e.g. Si³⁺ instead of Si⁴⁺).

    Disordered input is ordered first (see ``order_if_disordered``), because
    ``site.specie`` is undefined on a mixed-occupancy site.
    """
    s = structure.copy()

    # --- 0. Fractional occupancies must go before anything touches .specie ---
    if not s.is_ordered:
        s = order_if_disordered(s)

    # --- 1. Try pymatgen's guess, then its exhaustive mode, then OXI_MAP ---
    guessed = False
    for kwargs in ({}, {"max_sites": -1}):
        try:
            s.add_oxidation_state_by_guess(**kwargs)
            guessed = True
            break
        except Exception:
            continue

    if not guessed:
        # Both guesses failed (common for non-integer compositions such as
        # vacancy-bearing oxides). Assign directly from OXI_MAP instead of
        # letting the exception propagate.
        missing = sorted({el.symbol for el in s.composition.elements
                          if el.symbol not in OXI_MAP})
        if missing:
            raise ValueError(
                f"Could not assign oxidation states: pymatgen's guess failed "
                f"and these elements are absent from OXI_MAP: {missing}. "
                f"Add them to OXI_MAP in polycleaver2/slabgen.py."
            )
        s = Structure(
            s.lattice,
            [Species(site.specie.symbol if hasattr(site.specie, 'symbol')
                     else site.specie.element.symbol,
                     oxidation_state=OXI_MAP[
                         site.specie.symbol if hasattr(site.specie, 'symbol')
                         else site.specie.element.symbol])
             for site in s],
            s.frac_coords,
        )

    # --- 2. Check whether *any* site has a clearly wrong oxidation state. ---
    # A species is "wrong" if it appears in OXI_MAP but pymatgen's
    # assigned value doesn't match (e.g. Si³⁺ instead of Si⁴⁺).
    any_wrong = any(
        site.specie.oxi_state != OXI_MAP[site.specie.element.symbol]
        for site in s
        if site.specie.element.symbol in OXI_MAP
    )

    if any_wrong:
        fixed = []
        for site in s:
            sym = site.specie.element.symbol
            correct_oxi = OXI_MAP.get(sym)
            if correct_oxi is not None and site.specie.oxi_state != correct_oxi:
                fixed.append(Species(sym, oxidation_state=correct_oxi))
            else:
                fixed.append(site.specie)
        s = Structure(s.lattice, fixed, s.frac_coords)
    return s


def structure_to_slab(structure, ref_slab):
    """
    Convert a plain Structure back into a pymatgen Slab, preserving the
    Miller-index and oriented-unit-cell metadata from *ref_slab*.
    """
    return Slab(
        structure.lattice,
        structure.species,
        structure.frac_coords,
        ref_slab.miller_index,
        ref_slab.oriented_unit_cell,
        ref_slab.shift,
        ref_slab.scale_factor,
        site_properties=structure.site_properties,
    )


# =============================================================================
#  Configuration loader
# =============================================================================

def load_config(path="config.json"):
    """
    Load run parameters from a JSON file and return them as a dict.

    All keys are optional — missing ones fall back to sensible defaults.
    Keys whose names start with ``_`` are treated as inline documentation
    and ignored, so you can annotate your config.json freely.

    Parameters
    ----------
    path : str or Path
        Path to the JSON config file (default: ``config.json`` in the
        current working directory).

    Returns
    -------
    dict with keys:
        filename, target_thickness, vacuum_size, dipole_threshold,
        mode, single_hkl, max_index, custom_list

    Raises
    ------
    FileNotFoundError
        If the config file does not exist.
    ValueError
        If the JSON is malformed or a required parameter is invalid.
    """
    defaults = {
        "filename": None,
        "target_thickness": 15.0,
        "vacuum_size": 15.0,
        "dipole_threshold": 0.3,
        "mode": "single",
        "single_hkl": [0, 0, 1],
        "max_index": 2,
        "custom_list": [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        "thickness_range": None,
        "auto_select_stable": True,
        "allow_nonstoichiometric": False,
    }

    path = Path(path)
    try:
        raw = json.loads(path.read_text())
    except FileNotFoundError:
        raise FileNotFoundError(
            f"Config file not found: {path}\n"
            f"  Copy examples/configs/forsterite_110_range.json to config.json and set 'filename'."
        )
    except json.JSONDecodeError as exc:
        raise ValueError(f"Invalid JSON in {path}: {exc}") from exc

    # Overlay user values over defaults; skip documentation-comment keys
    cfg = {**defaults, **{k: v for k, v in raw.items() if not k.startswith("_")}}

    if not cfg["filename"]:
        raise ValueError(f"{path}: 'filename' must be set to your CIF file path.")
    if cfg["mode"] not in ("single", "max", "list"):
        raise ValueError(
            f"{path}: 'mode' must be 'single', 'max', or 'list'  (got {cfg['mode']!r})"
        )
    for key in ("target_thickness", "vacuum_size"):
        if cfg[key] is not None and cfg[key] <= 0:
            raise ValueError(f"{path}: '{key}' must be positive  (got {cfg[key]})")
    if cfg["target_thickness"] is None and cfg["thickness_range"] is None:
        raise ValueError(
            f"{path}: set either 'target_thickness' or 'thickness_range' (both are null)"
        )
    if cfg["dipole_threshold"] < 0:
        raise ValueError(
            f"{path}: 'dipole_threshold' must be ≥ 0  (got {cfg['dipole_threshold']})"
        )
    tr = cfg["thickness_range"]
    if tr is not None:
        if (not isinstance(tr, list) or len(tr) != 2
                or not all(isinstance(v, (int, float)) for v in tr)
                or tr[0] <= 0 or tr[1] <= tr[0]):
            raise ValueError(
                f"{path}: 'thickness_range' must be [min, max] with 0 < min < max "
                f"(got {tr!r})"
            )

    # Return Miller index fields as tuples so callers don't have to convert
    cfg["single_hkl"] = tuple(cfg["single_hkl"])
    cfg["custom_list"] = [tuple(x) for x in cfg["custom_list"]]
    return cfg


# =============================================================================
#  0.  SURFACE STABILITY UTILITIES
# =============================================================================

def _cn_of_site(site_sym, neighbors):
    """
    Count the coordination number of a site using element-pair-specific rules.

    Three cases:
    - Scaffold centre (Si, P, S, Al, …): count O neighbours with covalent
      bond distance < 2.3 Å.  Using the full 2.8 Å cutoff would pick up
      non-bonded O from adjacent polyhedra, inflating the apparent CN.
    - Oxygen: count all non-O, non-H neighbours within 2.8 Å (Si, cations,
      etc.).  This avoids counting intra-tetrahedral O–O contacts (~2.66 Å)
      which are close-packing artefacts, not bonds.
    - Cation (Mg, Ca, …): count O neighbours within 2.8 Å.

    Parameters
    ----------
    site_sym : str
        Element symbol of the central atom.
    neighbors : list
        Output of ``structure.get_neighbors`` / ``get_all_neighbors`` for
        this site; each neighbour exposes ``.specie.element.symbol`` and
        ``.nn_distance``.

    Returns
    -------
    int  coordination number
    """
    if site_sym == 'O':
        return sum(1 for n in neighbors
                   if n.specie.element.symbol not in ('O', 'H'))
    if site_sym in COVALENT_ELEMENTS:
        return sum(1 for n in neighbors
                   if n.specie.element.symbol == 'O'
                   and n.nn_distance < 2.3)
    # cation
    return sum(1 for n in neighbors
               if n.specie.element.symbol == 'O')


def _compute_bulk_cn_all(bulk):
    """
    Mean coordination number per element in the bulk structure, using
    element-pair-specific counting rules (see ``_cn_of_site``).

    Correctly gives CN ≈ 4 for Si (covalent bonds only), CN ≈ 4 for O
    (non-O non-H neighbours: 1 Si + n_cations), and the true ionic CN for
    cations (Mg → 6, Ca → 6–9, …).  H is excluded entirely.

    Returns
    -------
    dict  element → mean coordination number
    """
    stats = defaultdict(list)
    all_nn = bulk.get_all_neighbors(r=2.8)
    for i, site in enumerate(bulk):
        sym = site.specie.element.symbol
        if sym == 'H':
            continue
        stats[sym].append(_cn_of_site(sym, all_nn[i]))
    return {sym: float(np.mean(counts)) for sym, counts in stats.items()}


def score_surface_stability(slab, bulk_cn_all, surface_depth=3.5):
    """
    Score a slab's surface stability from the coordination saturation of
    ALL non-H surface atoms.

    Counting only cations misses the other side of every broken bond: when
    a Mg–O bond is severed at the surface, Mg loses an O neighbour AND that
    O atom loses a Mg neighbour.  Including scaffold centres (Si, P …) and
    anion O in the count captures both sides and gives a better proxy for
    the surface bond-breaking energy.

    A surface atom is any non-H atom within *surface_depth* Å of the
    topmost or bottommost atom in the slab.

    Returns
    -------
    (coord_saturation, n_undercoord, broken_per_area)
        coord_saturation : float  0–1  mean(CN_surf / CN_bulk) for surface atoms
                           — higher is more stable
        n_undercoord     : int    number of surface atoms with CN < bulk CN
                           — lower is more stable
        broken_per_area  : float  Σ(CN_bulk − CN_surf) for surface atoms / Å²
                           — lower is more stable
    """
    z = slab.cart_coords[:, 2]
    z_top, z_bot = float(z.max()), float(z.min())
    a_vec = np.array(slab.lattice.matrix[0])
    b_vec = np.array(slab.lattice.matrix[1])
    area = float(np.linalg.norm(np.cross(a_vec, b_vec)))

    all_nn = slab.get_all_neighbors(r=2.8)
    saturations = []
    n_undercoord = 0
    total_broken = 0

    for i, site in enumerate(slab):
        sym = site.specie.element.symbol
        if sym == 'H':
            continue
        zi = float(z[i])
        if not (zi <= z_bot + surface_depth or zi >= z_top - surface_depth):
            continue
        cn_surf = _cn_of_site(sym, all_nn[i])
        cn_bulk = bulk_cn_all.get(sym, cn_surf)
        if cn_bulk > 0:
            saturations.append(min(cn_surf / cn_bulk, 1.0))
            deficit = max(0, round(cn_bulk) - cn_surf)
            total_broken += deficit
            if deficit > 0:
                n_undercoord += 1

    coord_sat = float(np.mean(saturations)) if saturations else 0.0
    broken_per_area = total_broken / area if area > 0 else 0.0
    return coord_sat, n_undercoord, broken_per_area


def _deduplicate_slabs(slabs):
    """Remove slabs that share the same formula and thickness (within 0.1 Å)."""
    seen = set()
    unique = []
    for slab in slabs:
        key = (slab.formula, round(float(np.ptp(slab.cart_coords[:, 2])), 1))
        if key not in seen:
            seen.add(key)
            unique.append(slab)
    return unique


# =============================================================================
#  1.  CHEMICAL ANALYSER
# =============================================================================

class ChemicalAnalyzer:
    """
    Analyse a mineral Structure and classify every atom into one of:

    centres (M)
        The central atoms of the polyanion scaffold framework.
        These are covalent elements with positive oxidation state that have
        at least one oxygen neighbour (e.g. Si⁴⁺, P⁵⁺, S⁶⁺, Al³⁺ …).

    anions (O⁻)
        Oxygen atoms covalently bonded to one or more centres.  In
        isolated polyhedra each O is bonded to exactly one centre; in
        framework / sheet silicates they are shared (M–O–M bridges).

    cations (Catⁿ⁺)
        Ionic species (metals / alkali / alkaline-earth) that are *not*
        part of the covalent scaffold.  These are charge-balancing
        counter-ions (e.g. K⁺, Mg²⁺, NH₄⁺ …).

    molecules (H₂O, NH₃, OH⁻ …)
        Discrete molecular clusters identified via covalent bonding
        analysis.  A molecular cluster contains no centre atoms.

    The analyser also records:
      • covalent_map      centre → [centre, O₁, O₂, …]
      • centre_coordination   centre → number of bonded O
      • cn_by_species         element → expected coordination number
      • cation_hydration_target  element → modal water-O neighbours in bulk

    Usage
    -----
    an = ChemicalAnalyzer(bulk_structure)
    an.print_summary()
    """

    def __init__(self, structure):
        self.structure = structure
        self.centers = []                      # site indices of scaffold centres
        self.anions = []                       # site indices of O atoms in scaffold
        self.cations = []                      # site indices of counter-cations
        self.molecules = defaultdict(list)     # formula → list[list[int]]
        self.covalent_map = {}                 # c_idx → [c_idx, O_neighbour, …]
        self.center_coordination = {}          # c_idx → number of O neighbours
        self.cn_by_species = {}                # element → expected CN
        self.cation_hydration_target = {}      # cation elem → modal H₂O neighbours
        self.adj_list = None                   # covalent adjacency (site → set)

        if len(self.structure) > 0:
            self._analyze()

    # ------------------------------------------------------------------
    #  Core analysis
    # ------------------------------------------------------------------

    def _analyze(self):
        """
        Run the full chemical analysis pipeline on *self.structure*.

        1. Build a covalent adjacency list (neighbour radius 2.8 Å,
           X–H < 1.25 Å, other bonds < 2.3 Å).
        2. Find connected components (clusters) in the covalent graph.
        3. Cation = any non-covalent-element site.
        4. For each covalent cluster:
           a. If it contains at least one *centre* atom (a non-H/O
              covalent element with positive oxidation state and at
              least one O neighbour), extract centre–O units as
              scaffold; classify remaining cluster atoms as anions.
           b. Otherwise classify the whole cluster as a molecule.
        5. Compute per-species expected coordination numbers from
           the centres.
        6. Compute cation hydration statistics (modal water-O neighbours
           within 2.8 Å).
        """
        sites = self.structure.sites
        adj_list = defaultdict(set)
        self.adj_list = adj_list

        # --- 1. Covalent neighbour search (radius = 2.8 Å) ---
        try:
            all_nn = self.structure.get_all_neighbors(r=2.8)
        except Exception:
            all_nn = [[] for _ in sites]

        for i, neighbours in enumerate(all_nn):
            elem_i = sites[i].specie.element.symbol
            if elem_i not in COVALENT_ELEMENTS:
                continue
            for n in neighbours:
                elem_j = sites[n.index].specie.element.symbol
                if elem_j not in COVALENT_ELEMENTS:
                    continue
                dist = n.nn_distance
                # X–H bonds are shorter than other bonds
                threshold = 1.25 if (elem_i == 'H' or elem_j == 'H') else 2.3
                if dist < threshold:
                    adj_list[i].add(n.index)
                    adj_list[n.index].add(i)

        # --- 2. Connected-component clustering ---
        visited = set()
        clusters = []
        for i in range(len(self.structure)):
            elem = sites[i].specie.element.symbol
            if elem not in COVALENT_ELEMENTS:
                self.cations.append(i)
                visited.add(i)
                continue
            if i not in visited:
                component = []
                stack = [i]
                visited.add(i)
                while stack:
                    cur = stack.pop()
                    component.append(cur)
                    for nb in adj_list[cur]:
                        if nb not in visited:
                            visited.add(nb)
                            stack.append(nb)
                clusters.append(component)

        # --- 3. Determine which covalent species count as centres ---
        # Centres = non-H/O covalent elements with positive oxidation state.
        try:
            unique_species = sorted(
                list({s.specie for s in self.structure}),
                key=lambda x: x.oxi_state,
            )
            centre_symbols = {
                sp.element.symbol for sp in unique_species
                if sp.element.symbol not in ('H', 'O')
                   and sp.oxi_state > 0
                   and sp.element.symbol in COVALENT_ELEMENTS
            }
        except Exception:
            centre_symbols = set()

        # --- 4. Classify each cluster as scaffold or molecule ---
        all_molecular_indices = set()
        for cluster in clusters:
            c_idxs = [
                x for x in cluster
                if sites[x].specie.element.symbol in centre_symbols
                   and any(sites[nb].specie.element.symbol == 'O'
                           for nb in adj_list[x])
            ]
            if c_idxs:
                # ── This cluster is a scaffold ──
                c_set = set(c_idxs)
                for x in cluster:
                    if x in c_set:
                        local_o = [nb for nb in adj_list[x]
                                   if sites[nb].specie.element.symbol == 'O']
                        # Structural hydroxyl H bonded to this centre's O
                        # (kaolinite, micas, amphiboles …) travels with the
                        # unit, so stripping/thinning never orphans it.
                        local_h = [
                            nb2 for o in local_o for nb2 in adj_list[o]
                            if sites[nb2].specie.element.symbol == 'H'
                        ]
                        self.centers.append(x)
                        self.covalent_map[x] = [x] + local_o + local_h
                        self.center_coordination[x] = len(local_o)
                    else:
                        self.anions.append(x)
            else:
                # ── This cluster is a discrete molecule ──
                elements = [sites[x].specie.element.symbol for x in cluster]
                counts = Counter(elements)
                formula = "".join(
                    sorted(f"{k}{v if v > 1 else ''}" for k, v in counts.items())
                )
                self.molecules[formula].append(cluster)
                all_molecular_indices.update(cluster)

        # --- 5. Per-species coordination numbers ---
        for c_idx in self.centers:
            sym = sites[c_idx].specie.element.symbol
            if sym not in self.cn_by_species:
                self.cn_by_species[sym] = self.center_coordination.get(c_idx, 4)

        # --- 6. Cation hydration statistics ---
        molecular_oxygens = {
            idx for idx in all_molecular_indices
            if sites[idx].specie.element.symbol == 'O'
        }
        if self.cations and molecular_oxygens:
            # Number of molecular O (H2O, OH-, …) within 2.8 Å of each cation;
            # the modal value per element is the bulk hydration number
            # (e.g. 6 for the hexaaqua Mg2+ in struvite).
            stats = defaultdict(list)
            for cat_idx in self.cations:
                cat_sym = sites[cat_idx].specie.element.symbol
                n_mol_o = sum(
                    1 for n in self.structure.get_neighbors(sites[cat_idx], 2.8)
                    if n.index in molecular_oxygens
                )
                stats[cat_sym].append(n_mol_o)
            for sym, counts in stats.items():
                self.cation_hydration_target[sym] = (
                    Counter(counts).most_common(1)[0][0] if counts else 0
                )
        for cat_idx in self.cations:
            sym = sites[cat_idx].specie.element.symbol
            if sym not in self.cation_hydration_target:
                self.cation_hydration_target[sym] = 0

    # ------------------------------------------------------------------
    #  Public helpers
    # ------------------------------------------------------------------

    def formula_summary(self):
        """Return a human-readable breakdown string."""
        parts = []
        for c in sorted(set(self.structure[i].specie.element.symbol
                            for i in self.centers)):
            n = sum(1 for i in self.centers
                    if self.structure[i].specie.element.symbol == c)
            parts.append(f"{c}{n}")
        return "Scaffold: " + " + ".join(parts) if parts else "No scaffold"

    def is_polymerized(self):
        """
        True if the scaffold centres form a *single connected* network
        via shared oxygen bridges (i.e. a framework or sheet silicate).
        """
        if not self.centers or not self.adj_list:
            return False
        # Build centre–centre adjacency through O bridges
        centre_set = set(self.centers)
        cx_adj = defaultdict(set)
        for c1 in self.centers:
            for n_idx in self.adj_list[c1]:
                if self.structure[n_idx].specie.element.symbol != 'O':
                    continue
                for n2_idx in self.adj_list[n_idx]:
                    if n2_idx in centre_set and n2_idx != c1:
                        cx_adj[c1].add(n2_idx)
        if not cx_adj:
            return False
        # BFS from first centre
        start = self.centers[0]
        seen = {start}
        stack = [start]
        while stack:
            cur = stack.pop()
            for nb in cx_adj[cur]:
                if nb not in seen:
                    seen.add(nb)
                    stack.append(nb)
        return len(seen) == len(self.centers)

    def has_shared_oxygens(self):
        """
        True if any oxygen in the structure is bonded to *two or more*
        scaffold centres (i.e. M–O–M bridges exist).  This is the
        defining characteristic of framework and sheet silicates.
        Unlike *is_polymerized* it does not require all centres to be
        in a single connected component — important for layered /
        phyllosilicate bulks that contain multiple disconnected T-O-T
        sheets but are still chemically frameworks.
        """
        if not self.centers or not self.adj_list:
            return False
        centre_idx = set(self.centers)
        o_adjacent_centres = defaultdict(set)
        for c in self.centers:
            for nb in self.adj_list[c]:
                if self.structure[nb].specie.element.symbol == 'O':
                    o_adjacent_centres[nb].add(c)
        return any(len(cs) > 1 for cs in o_adjacent_centres.values())

    def print_summary(self):
        """Short one-line summary for the user."""
        if self.centers:
            print(f"  > Scaffold         : {len(self.centers)} centres detected")
        else:
            print(f"  > Scaffold         : None detected")
        print(f"  > Anions (O)       : {len(self.anions)}")
        print(f"  > Cations          : {len(self.cations)} detected" if self.cations
              else "  > Cations          : None")
        if self.molecules:
            for form, units in self.molecules.items():
                print(f"  > Molecule {form:<8s}: {len(units)} units")

    def print_detailed(self):
        """Comprehensive bulk-analysis report."""
        sites = self.structure.sites

        print("  ┌─ Bulk composition ──────────────────────────────────────")
        print(f"  │ Formula              : {self.structure.formula}")
        print(f"  │ Space group (guess)  : {self.structure.get_space_group_info()[0]}")
        lat = self.structure.lattice
        print(f"  │ Lattice              : a={lat.a:.4f}  b={lat.b:.4f}  c={lat.c:.4f}  "
              f"α={lat.alpha:.2f} β={lat.beta:.2f} γ={lat.gamma:.2f}")
        print(f"  │ Volume               : {lat.volume:.2f} Å³")
        print(f"  │ Total atoms          : {len(self.structure)}")
        print("  ├─ Poly-anion scaffold ───────────────────────────────────")

        if self.centers:
            by_sym = Counter(sites[i].specie.element.symbol for i in self.centers)
            for sym, cnt in sorted(by_sym.items()):
                c_idx = next(i for i in self.centers
                             if sites[i].specie.element.symbol == sym)
                oxi = sites[c_idx].specie.oxi_state
                cn_vals = [self.center_coordination.get(i) for i in self.centers
                           if sites[i].specie.element.symbol == sym]
                cn_str = f", CN = {cn_vals[0]}" if cn_vals and cn_vals[0] is not None else ""
                print(f"  │   Centre {sym:<4s}    : {cnt:3d} ×  "
                      f"(oxidation = {oxi:+.0f}{cn_str})")
            shared = self.has_shared_oxygens()
            print(f"  │   Connectivity      : "
                  f"{'Framework/sheet (M–O–M bridges)' if shared else 'Isolated polyhedra'}")
        else:
            print("  │   (none detected)")

        if self.anions:
            by_sym = Counter(sites[i].specie.element.symbol for i in self.anions)
            for sym, cnt in sorted(by_sym.items()):
                print(f"  │   Anion {sym:<5s}    : {cnt:3d} ×")
        print("  ├─ Counter-cations ───────────────────────────────────────")
        if self.cations:
            by_sym = Counter(sites[i].specie.element.symbol for i in self.cations)
            for sym, cnt in sorted(by_sym.items()):
                c_idx = next(i for i in self.cations
                             if sites[i].specie.element.symbol == sym)
                oxi = sites[c_idx].specie.oxi_state if hasattr(
                    sites[c_idx].specie, 'oxi_state') else 0
                hydr = self.cation_hydration_target.get(sym, '—')
                print(f"  │   {sym:<6s}    : {cnt:3d} ×  (oxidation = {oxi:+.0f}, "
                      f"Hydration target = {hydr} H₂O)")
        else:
            print("  │   (none detected)")

        print("  ├─ Molecular species ─────────────────────────────────────")
        if self.molecules:
            for formula, units in sorted(self.molecules.items()):
                n_atoms = sum(len(u) for u in units)
                composition = Counter(sites[i].specie.element.symbol
                                      for unit in units for i in unit)
                comp_str = " ".join(f"{s}{n}" for s, n in sorted(composition.items()))
                print(f"  │   {formula:<8s}  : {len(units):3d} units  ({n_atoms} atoms: {comp_str})")
        else:
            print("  │   (none detected)")

        n_c = len(self.centers) or 1
        cat_ratio = len(self.cations) / n_c
        print("  ├─ Stoichiometry ─────────────────────────────────────────")
        print(f"  │   Cations per centre : {cat_ratio:.4f}")
        for form, units in self.molecules.items():
            print(f"  │   {form} per centre    : {len(units) / n_c:.4f}")
        print("  └─────────────────────────────────────────────────────────")


# =============================================================================
#  2.  SLAB PROCESSOR  (stripping & thinning)
# =============================================================================

class SlabProcessor:
    """
    Manages stripping a slab down to its covalent scaffold and then thinning
    it to the target thickness while controlling the z-dipole moment.

    Parameters
    ----------
    slab : Slab
        The full (un-stripped) slab, already orthogonalised along c.
    bulk_analyzer : ChemicalAnalyzer
        Analysis of the bulk structure, used to obtain per-species
        coordination numbers and stoichiometric ratios.
    """

    def __init__(self, slab, bulk_analyzer):
        self.slab = slab
        self.bulk_analyzer = bulk_analyzer
        self.local_analysis = ChemicalAnalyzer(self.slab)

        n_bulk_centers = len(bulk_analyzer.centers)
        if n_bulk_centers == 0:
            self.ratio_cation = 0.0
            self.ratio_mols = {}
        else:
            self.ratio_cation = len(bulk_analyzer.cations) / n_bulk_centers
            self.ratio_mols = {
                formula: len(units) / n_bulk_centers
                for formula, units in bulk_analyzer.molecules.items()
            }

    # ------------------------------------------------------------------
    #  Layer detection (framework silicates)
    # ------------------------------------------------------------------

    def _is_framework(self):
        """
        Detect whether this slab exhibits framework-like connectivity — i.e.
        any oxygen atom is shared between *two or more* centres (an M–O–M
        bridge).  Used during thinning: framework slabs must remove entire
        crystallographic z-layers at once, while isolated-polyhedra slabs
        can remove single units.
        """
        if not self.local_analysis:
            return False
        o_count = defaultdict(int)
        for unit in self.local_analysis.covalent_map.values():
            for idx in unit:
                if self.slab[idx].specie.element.symbol == 'O':
                    o_count[idx] += 1
        return any(c > 1 for c in o_count.values())

    def _get_center_layers(self):
        """
        Group centres into crystallographic z-layers.

        The algorithm sorts centres by their Cartesian z-coordinate,
        computes all consecutive z-gaps, then derives a threshold:
            threshold = max(1.5, median(gaps) × 2.5)
        Splits wherever gap > threshold.  This cleanly separates the
        ～10 Å inter-sheet gaps in phyllosilicates from the ～1–3 Å
        intra-sheet spread.

        Returns
        -------
        list of [avg_z, set_of_centre_indices], sorted bottom → top.
        """
        cart = self.slab.cart_coords
        sorted_centers = sorted(
            [(float(cart[c][2]), c) for c in self.local_analysis.covalent_map],
            key=lambda x: x[0],
        )
        if not sorted_centers:
            return []

        gaps = [sorted_centers[i + 1][0] - sorted_centers[i][0]
                for i in range(len(sorted_centers) - 1)]
        threshold = max(1.5, float(np.median(gaps)) * 2.5) if gaps else 1.5

        layers = []
        layer_zs = [sorted_centers[0][0]]
        layer_idxs = {sorted_centers[0][1]}
        for i in range(1, len(sorted_centers)):
            z, c = sorted_centers[i]
            if (z - sorted_centers[i - 1][0]) > threshold:
                layers.append([float(np.mean(layer_zs)), layer_idxs])
                layer_zs = [z]
                layer_idxs = {c}
            else:
                layer_zs.append(z)
                layer_idxs.add(c)
        layers.append([float(np.mean(layer_zs)), layer_idxs])
        return layers

    # ------------------------------------------------------------------
    #  Stripping
    # ------------------------------------------------------------------

    def remove_indices(self, indices_to_remove):
        """Remove a set of site indices from the slab and re-analyse."""
        if not indices_to_remove:
            return
        self.slab.remove_sites(list(indices_to_remove))
        if self.slab:
            self.local_analysis = ChemicalAnalyzer(self.slab)
        else:
            self.local_analysis = None

    def step_5_strip_to_scaffold(self):
        """
        Strip the slab to its bare covalent scaffold.

        1. Remove all molecular units and counter-cations.
        2. If the structure is *not* polymerised (i.e. has isolated
           polyhedra like orthosilicates, phosphates …), also remove
           surface-incomplete polyhedra whose O-coordination number
           is lower than the bulk value.  This ensures only
           stoichiometrically complete scaffold units remain.
        3. Cascade-remove any oxygen atoms left orphaned by the
           removal of their parent centres.

        For polymerised structures (framework / sheet silicates) the
        under-coordinated surface tetrahedra are *kept* — they are
        inherently broken at the cleavage plane and must be preserved.
        """
        to_remove = set()
        # 1. Molecules + cations
        for mol_list in self.local_analysis.molecules.values():
            for unit in mol_list:
                to_remove.update(unit)
        to_remove.update(self.local_analysis.cations)

        # 2. Incomplete polyhedra (only for non-framework scaffolds)
        if self.bulk_analyzer.centers and not self.bulk_analyzer.has_shared_oxygens():
            incomplete = set()
            for c_idx in self.local_analysis.covalent_map:
                sym = self.slab[c_idx].specie.element.symbol
                target_cn = self.bulk_analyzer.cn_by_species.get(sym, 4)
                current_cn = self.local_analysis.center_coordination.get(c_idx, 0)
                if current_cn < target_cn:
                    incomplete.add(c_idx)
            to_remove.update(incomplete)

            # Cascade: O atoms that are no longer bonded to any remaining centre
            remaining_centers = set(self.local_analysis.covalent_map) - incomplete
            bonded_O = set()
            for c_idx in remaining_centers:
                bonded_O.update(self.local_analysis.covalent_map[c_idx])
            for a_idx in self.local_analysis.anions:
                if a_idx not in bonded_O:
                    to_remove.add(a_idx)

        # 3. Any atom that slipped through classification
        classified = set(self.local_analysis.cations)
        classified.update(self.local_analysis.anions)
        for u in self.local_analysis.covalent_map.values():
            classified.update(u)
        for mlist in self.local_analysis.molecules.values():
            for u in mlist:
                classified.update(u)
        for i in range(len(self.slab)):
            if i not in classified:
                to_remove.add(i)

        self.remove_indices(to_remove)

    # ------------------------------------------------------------------
    #  Thinning
    # ------------------------------------------------------------------

    def step_5b_thin_to_target(self, target_thickness, dipole_tol):
        """
        Thin the scaffold to find the minimum-period non-polar termination.

        Removes outermost scaffold layers iteratively (framework: whole
        layers; isolated polyhedra: one unit per side) until no further
        removal is possible.  At every thickness where |dipole_z| ≤
        *dipole_tol*, the slab state is recorded.  After exhaustion:

        - If any valid state has thickness ≥ *target_thickness*, the
          **thinnest** such state is selected (smallest slab that still
          meets the user's minimum thickness requirement).
        - Otherwise the **thickest** valid state is selected — the
          irreducible minimum period for this termination.  A warning is
          printed by the caller.

        Returns True if at least one valid state was found, False if the
        scaffold has no configuration with |dipole_z| ≤ *dipole_tol*.
        """
        if not self.local_analysis or not self.local_analysis.covalent_map:
            return False

        valid_states = []  # list of (thickness, slab_copy)
        # For framework silicates the scaffold is stoichiometrically self-contained,
        # so its own dipole must be small.  For isolated-polyhedra minerals (e.g.
        # orthosilicates, phosphates) the scaffold alone is always charged, and the
        # counter-cations cancel the dipole — so we collect every thickness state
        # and let Phase 1 apply the dipole filter.
        is_framework = self._is_framework()

        while True:
            cart_coords = self.slab.cart_coords
            current_thickness = np.ptp(cart_coords[:, 2])
            current_dipole = abs(self.slab.dipole[2])

            if is_framework:
                if current_dipole <= dipole_tol:
                    valid_states.append((current_thickness, self.slab.copy()))
            else:
                valid_states.append((current_thickness, self.slab.copy()))

            if not self.local_analysis.centers:
                break

            if is_framework:
                layers = self._get_center_layers()
                if len(layers) < 2:
                    break
                centers_to_remove = layers[0][1] | layers[-1][1]
            else:
                # Isolated polyhedra: outermost units from fractional‑z extremes
                centre_z = sorted(
                    [(self.slab[c].frac_coords[2], c)
                     for c in self.local_analysis.covalent_map],
                    key=lambda x: x[0],
                )
                if len(centre_z) < 2:
                    break
                centers_to_remove = {centre_z[0][1], centre_z[-1][1]}

            # Cascade-remove orphaned O atoms
            remaining = set(self.local_analysis.covalent_map) - centers_to_remove
            if not remaining:
                break
            bonded_to_remaining = set()
            for c_idx in remaining:
                bonded_to_remaining.update(self.local_analysis.covalent_map[c_idx])

            to_remove = set(centers_to_remove)
            for a_idx in self.local_analysis.anions:
                if a_idx not in bonded_to_remaining:
                    to_remove.add(a_idx)

            self.remove_indices(to_remove)
            if not self.local_analysis or not self.local_analysis.centers:
                break

        if not valid_states:
            return False

        # Prefer the thinnest valid slab that is still ≥ target_thickness.
        # Fall back to the thickest valid state (irreducible minimum period)
        # when no state meets the target.
        thick_enough = [(t, s) for t, s in valid_states if t >= target_thickness]
        _, best_slab = (min(thick_enough, key=lambda x: x[0]) if thick_enough
                        else max(valid_states, key=lambda x: x[0]))

        self.slab = best_slab
        self.local_analysis = ChemicalAnalyzer(self.slab)
        return True

    def step_5b_get_all_valid_states(self, min_thick, max_thick, dipole_tol):
        """
        Return ALL scaffold states in [min_thick, max_thick] with |dipole_z| ≤ dipole_tol.

        Thins the scaffold destructively (same logic as step_5b_thin_to_target)
        but collects every qualifying state instead of selecting one.  Stops
        early once thickness drops below min_thick − 2 Å.

        Returns
        -------
        list of (thickness, slab_copy) sorted thickest → thinnest.
        """
        if not self.local_analysis or not self.local_analysis.covalent_map:
            return []

        valid_states = []
        is_framework = self._is_framework()

        while True:
            cart_coords = self.slab.cart_coords
            current_thickness = float(np.ptp(cart_coords[:, 2]))
            current_dipole = abs(self.slab.dipole[2])

            in_range = min_thick <= current_thickness <= max_thick
            if in_range:
                if is_framework:
                    if current_dipole <= dipole_tol:
                        valid_states.append((current_thickness, self.slab.copy()))
                else:
                    valid_states.append((current_thickness, self.slab.copy()))

            if current_thickness < min_thick - 2.0:
                break

            if not self.local_analysis.centers:
                break

            if is_framework:
                layers = self._get_center_layers()
                if len(layers) < 2:
                    break
                centers_to_remove = layers[0][1] | layers[-1][1]
            else:
                centre_z = sorted(
                    [(self.slab[c].frac_coords[2], c)
                     for c in self.local_analysis.covalent_map],
                    key=lambda x: x[0],
                )
                if len(centre_z) < 2:
                    break
                centers_to_remove = {centre_z[0][1], centre_z[-1][1]}

            remaining = set(self.local_analysis.covalent_map) - centers_to_remove
            if not remaining:
                break
            bonded_to_remaining = set()
            for c_idx in remaining:
                bonded_to_remaining.update(self.local_analysis.covalent_map[c_idx])

            to_remove = set(centers_to_remove)
            for a_idx in self.local_analysis.anions:
                if a_idx not in bonded_to_remaining:
                    to_remove.add(a_idx)

            self.remove_indices(to_remove)
            if not self.local_analysis or not self.local_analysis.centers:
                break

        return sorted(valid_states, key=lambda x: x[0], reverse=True)


# =============================================================================
#  3.  MAIN GENERATION PIPELINE
# =============================================================================

def _make_slab_from_sites(sites, lattice, ref_slab):
    """Build a Slab from a flat list of PeriodicSite objects."""
    return structure_to_slab(
        Structure(lattice, [s.specie for s in sites], [s.frac_coords for s in sites]),
        ref_slab,
    )


_PROGRESS_INTERVAL = 5_000   # print a progress line every N combos when search is large
_STAGNATION_PATIENCE = 8_000  # abort after this many combos with no meaningful improvement
_STAGNATION_MIN_IMPROVEMENT = 1e-3  # improvement threshold to reset the stagnation counter


def _combinatorial_search(pool, req_n, base_sites, lattice, ref_slab, dipole_tol,
                           total_combos=None):
    """
    Find the combination of *req_n* units from *pool* that minimises |dipole_z|.

    Parameters
    ----------
    pool : list[list[PeriodicSite]]
        Candidate groups. Each element is a list of sites forming one unit.
        Wrap single-atom cations as ``[[s] for s in candidates]``.
    req_n : int
        Number of units to select.
    base_sites : list[PeriodicSite]
        Fixed sites already placed in the slab (scaffold + earlier species).
    lattice, ref_slab
        Passed to ``_make_slab_from_sites``.
    dipole_tol : float
        Early-exit threshold for |dipole_z| in eÅ.
    total_combos : int or None
        Pre-computed C(len(pool), req_n) for progress reporting.

    Returns
    -------
    (best_combo, best_dipole, best_slab)
    """
    best_combo = None
    best_dipole = float('inf')
    best_slab = None
    n_total = total_combos or math.comb(len(pool), req_n)
    report = n_total >= _PROGRESS_INTERVAL * 2   # only print for large searches
    last_improvement_step = 0

    for step, combo in enumerate(itertools.combinations(pool, req_n)):
        trial_sites = base_sites + [site for unit in combo for site in unit]
        trial_slab = _make_slab_from_sites(trial_sites, lattice, ref_slab)
        d = abs(trial_slab.dipole[2])
        if d < best_dipole - _STAGNATION_MIN_IMPROVEMENT:
            last_improvement_step = step
        if d < best_dipole:
            best_dipole = d
            best_combo = combo
            best_slab = trial_slab
        if d <= dipole_tol:
            break
        if report and (step + 1) % _PROGRESS_INTERVAL == 0:
            print(f"│     ... {step + 1:,}/{n_total:,} combos  best so far = {best_dipole:.4f} eÅ",
                  flush=True)
        # Abort early if stuck above tolerance with no recent improvement.
        # Factor 1.2 keeps searching when within 20% of target; aborts otherwise.
        if (step - last_improvement_step >= _STAGNATION_PATIENCE
                and best_dipole > 1.2 * dipole_tol):
            if report:
                print(f"│     ... stagnated at {best_dipole:.4f} eÅ after {step + 1:,} combos"
                      f" — aborting early", flush=True)
            break

    return best_combo, best_dipole, best_slab


def _run_phases_on_scaffold(scaffold, proc, cat_cands, mol_cands_map,
                            proc_slab, dipole_tol,
                            allow_nonstoichiometric=False):
    """
    Run Phase 1 (cation placement) and Phase 2 (molecule placement) on a
    given scaffold slab.

    Parameters
    ----------
    scaffold : Slab
        The stripped + thinned scaffold (proc.slab must equal this).
    proc : SlabProcessor
        Processor whose local_analysis reflects *scaffold*.
    cat_cands : list[PeriodicSite]
        All cation candidate sites captured from the full raw slab.
        Sorted in-place by proximity to the scaffold CoM each call.
    mol_cands_map : dict[str, list[list[PeriodicSite]]]
        Molecular candidate groups from the full raw slab.
    proc_slab : Slab
        Full orthogonalised slab — used as the Slab metadata reference.
    dipole_tol : float
        Maximum acceptable |dipole_z| in eÅ.
    allow_nonstoichiometric : bool
        If True, accept slabs with net charge ≠ 0 (non-stoichiometric).
        The net charge is stored in _nonstoich_charges for later annotation.

    Returns
    -------
    Slab if both phases succeed and |dipole_z| ≤ dipole_tol, else None.
    """
    # ------------------------------------------------------------------
    #  Phase 1 — Add cations
    # ------------------------------------------------------------------
    ase_scaffold = AseAtomsAdaptor.get_atoms(scaffold)
    com_frac_z = scaffold.lattice.get_fractional_coords(
        ase_scaffold.get_center_of_mass())[2]

    n_centers = len(proc.local_analysis.centers)
    req_cats = int(round(n_centers * proc.ratio_cation))

    cat_cands.sort(key=lambda s: abs(s.frac_coords[2] - com_frac_z))

    if len(cat_cands) < req_cats:
        print(f"│     ✕ Not enough cations: need {req_cats}, have {len(cat_cands)}")
        return None

    pool_cat = cat_cands[:min(len(cat_cands), req_cats + 4)]
    n_combinations_checked = math.comb(len(pool_cat), req_cats)

    print(f"│     Phase-1 cations: req={req_cats}, "
          f"pool={len(pool_cat)} → {n_combinations_checked:,} combos …", flush=True)

    _, lowest_dipole, best_base = _combinatorial_search(
        [[s] for s in pool_cat], req_cats,
        [s for s in scaffold], scaffold.lattice, proc_slab, dipole_tol,
        total_combos=n_combinations_checked,
    )

    print(f"│     Phase-1 result : best dipole = {lowest_dipole:.4f} eÅ")

    if lowest_dipole > dipole_tol:
        print(f"│     ✕ Cation search failed to reach dipole ≤ {dipole_tol}")
        return None
    print(f"│     ✓ Cation search accepted (dipole = {lowest_dipole:.4f} eÅ)")

    # ------------------------------------------------------------------
    #  Phase 2 — Add molecules
    # ------------------------------------------------------------------
    ase_base = AseAtomsAdaptor.get_atoms(best_base)
    ase_cell = ase_base.get_cell()
    ase_pbc = ase_base.get_pbc()
    mol_center_z = best_base.lattice.get_fractional_coords(
        ase_base.get_center_of_mass())[2]

    ase_cation_indices = [
        idx for idx, atom in enumerate(ase_base)
        if atom.symbol not in COVALENT_ELEMENTS
    ]
    cation_positions = ase_base.get_positions()[ase_cation_indices]

    def _water_key(unit):
        o_site = next((s for s in unit if s.specie.symbol == 'O'), unit[0])
        _, dists = get_distances(np.array([o_site.coords]), cation_positions,
                                 cell=ase_cell, pbc=ase_pbc)
        return np.min(dists)

    final_sites_list = [s for s in best_base]

    for form, ratio in proc.ratio_mols.items():
        req_n = int(round(n_centers * ratio))
        if req_n == 0:
            continue

        cands = mol_cands_map.get(form, [])
        if len(cands) < req_n:
            print(f"│     ✕ Not enough {form}: need {req_n}, have {len(cands)}")
            return None

        if 'H' in form and 'O' in form and len(ase_cation_indices) > 0:
            cands.sort(key=_water_key)
        else:
            cands.sort(key=lambda unit: abs(unit[0].frac_coords[2] - mol_center_z))

        mol_pool = cands[:min(len(cands), req_n + 4)]
        n_mol_combos = math.comb(len(mol_pool), req_n)

        best_mol_combo, best_mol_dipole, _ = _combinatorial_search(
            mol_pool, req_n,
            list(final_sites_list), scaffold.lattice, proc_slab, dipole_tol,
            total_combos=n_mol_combos,
        )

        final_sites_list.extend(site for unit in best_mol_combo for site in unit)
        print(f"│     Phase-2 {form:<6s}: need {req_n}, pool={len(mol_pool)} "
              f"→ {n_mol_combos} combos, best dipole = {best_mol_dipole:.4f} eÅ")

    # ------------------------------------------------------------------
    #  Build final slab and check dipole
    # ------------------------------------------------------------------
    final_slab = _make_slab_from_sites(final_sites_list, scaffold.lattice, proc_slab)
    final_dipole = abs(final_slab.dipole[2])

    if final_dipole > dipole_tol:
        print(f"│     ✕ Final dipole too large: {final_dipole:.4f} > {dipole_tol}")
        return None

    net_charge = sum(
        site.specie.oxi_state for site in final_slab
        if hasattr(site.specie, 'oxi_state')
    )
    if abs(net_charge) > 0.01:
        if not allow_nonstoichiometric:
            print(f"│     ✕ Non-stoichiometric slab (net charge = {net_charge:+.2f} e) — discarded")
            return None
        print(f"│     ! Non-stoichiometric (net charge = {net_charge:+.2f} e) — kept (allow_nonstoichiometric=true)")
        _nonstoich_charges[id(final_slab)] = net_charge

    print(f"│     ✓ SLAB ACCEPTED.  Final dipole = {final_dipole:.4f} eÅ")
    return final_slab


def generate_slabs_for_hkl(bulk, bulk_an, hkl,
                           target_thickness=23, vacuum=15, dipole_tol=2.0,
                           thickness_range=None, allow_nonstoichiometric=False):
    """
    Generate non-polar surface slabs for a single Miller index *hkl*.

    Pipeline
    --------
    1. Generate raw slabs with SlabGenerator (thickness = target + 10 Å buffer).
    2. For each cut:
       a. Convert to orthogonal-c representation.
       b. Capture candidate cations and molecules from the raw slab.
       c. Strip to bare scaffold.
       d. Thin to target thickness (single-target mode) OR collect all valid
          scaffold states in [min, max] (range mode).
       e. Phase 1: re-add cations — combinatorial search over the pool of
          candidates closest to the scaffold centre-of-mass.
       f. Phase 2: re-add molecules — sort by distance to cations (water) or
          distance to slab centre (other).
       g. Accept only if the final |dipole_z| ≤ *dipole_tol*.

    Parameters
    ----------
    bulk : Structure
        Bulk structure (with oxidation states).
    bulk_an : ChemicalAnalyzer
        Pre-computed bulk analysis.
    hkl : tuple[int, int, int]
        Miller index of the surface.
    target_thickness : float
        Desired slab thickness in Å (lower-bound preference, single-target mode).
    vacuum : float
        Vacuum thickness in Å.
    dipole_tol : float
        Maximum acceptable absolute z-dipole (eÅ).
    thickness_range : list[float] or None
        If [min, max], scan ALL valid scaffold thicknesses in that window and
        run Phase 1+2 for each, potentially yielding multiple terminations.
        If None, use single-target mode with *target_thickness*.

    Returns
    -------
    list[Slab] — accepted slab structures.
    """
    buffer_thickness = 10.0
    if thickness_range is not None:
        gen_thickness = thickness_range[1] + buffer_thickness
        range_str = f"[{thickness_range[0]:.1f}, {thickness_range[1]:.1f}] Å"
    else:
        gen_thickness = target_thickness + buffer_thickness
        range_str = f"target {target_thickness:.1f} Å"

    print(f"\n┌─ Processing Miller index {hkl} ".ljust(68, '─') + "┐")
    print(f"│  Generation thickness: {gen_thickness:.1f} Å  "
          f"({range_str} + buffer {buffer_thickness:.0f})")

    try:
        slabgen = SlabGenerator(bulk, hkl, gen_thickness, vacuum, center_slab=True)
        raw_slabs = slabgen.get_slabs()
    except Exception as e:
        print(f"│  ERROR: {e}")
        print("└" + "─" * 68 + "┘")
        return []

    print(f"│  Unique cuts found  : {len(raw_slabs)}")
    accepted = []

    for i, raw_slab in enumerate(raw_slabs):
        print(f"│  ─── Cut {i} ───────────────────────────────────────────────────")

        # ------------------------------------------------------------------
        #  Prepare: orthogonalise + capture candidates
        # ------------------------------------------------------------------
        proc_slab = raw_slab.get_orthogonal_c_slab()
        raw_an = ChemicalAnalyzer(proc_slab)

        cat_cands = [proc_slab[x] for x in raw_an.cations]
        mol_cands_map = {
            form: [[proc_slab[x] for x in unit] for unit in units]
            for form, units in raw_an.molecules.items()
        }
        print(f"│     Cation pool    : {len(cat_cands)}  |  Molecular pool: "
              f"{sum(len(v) for v in mol_cands_map.values())} units "
              f"({dict((k, len(v)) for k, v in mol_cands_map.items())})")

        # ------------------------------------------------------------------
        #  A. Strip to scaffold
        # ------------------------------------------------------------------
        proc = SlabProcessor(proc_slab.copy(), bulk_an)
        proc.step_5_strip_to_scaffold()
        if not proc.local_analysis or not proc.local_analysis.centers:
            print(f"│     ✕ Scaffold invalid or empty after stripping")
            continue

        print(f"│     Strip OK       : {len(proc.local_analysis.centers)} centres, "
              f"{len(proc.local_analysis.anions)} O")

        # ------------------------------------------------------------------
        #  B. Thin + Phase 1+2
        # ------------------------------------------------------------------
        if thickness_range is not None:
            # Range mode — find all valid scaffold states, run phases for each
            min_thick, max_thick = thickness_range
            valid_scaffolds = proc.step_5b_get_all_valid_states(
                min_thick, max_thick, dipole_tol)

            if not valid_scaffolds:
                print(f"│     ✕ No valid scaffold in [{min_thick:.1f}, {max_thick:.1f}] Å")
                continue

            thicknesses_str = ", ".join(f"{t:.2f}" for t, _ in valid_scaffolds)
            print(f"│     Found {len(valid_scaffolds)} scaffold state(s): {thicknesses_str} Å")

            for sc_thick, sc_slab in valid_scaffolds:
                proc.slab = sc_slab
                proc.local_analysis = ChemicalAnalyzer(sc_slab)
                scart = sc_slab.cart_coords
                print(f"│     ── Scaffold {sc_thick:.2f} Å, "
                      f"dipole = {abs(sc_slab.dipole[2]):.4f} eÅ")
                print(f"│        {len(proc.local_analysis.centers)} centres, "
                      f"{len(proc.local_analysis.anions)} O")
                result = _run_phases_on_scaffold(
                    sc_slab, proc, list(cat_cands),
                    {k: [list(u) for u in v] for k, v in mol_cands_map.items()},
                    proc_slab, dipole_tol,
                    allow_nonstoichiometric=allow_nonstoichiometric,
                )
                if result is not None:
                    accepted.append(result)
        else:
            # Single-target mode — original behaviour
            if not proc.step_5b_thin_to_target(target_thickness, dipole_tol):
                print(f"│     ✕ Could not find valid scaffold "
                      f"with |dipole| ≤ {dipole_tol:.2f}")
                continue

            scaffold = proc.slab
            scart = scaffold.cart_coords
            print(f"│     Thin OK        : thickness = {np.ptp(scart[:, 2]):.2f} Å, "
                  f"dipole = {abs(scaffold.dipole[2]):.4f} eÅ")

            result = _run_phases_on_scaffold(
                scaffold, proc, cat_cands, mol_cands_map, proc_slab, dipole_tol,
                allow_nonstoichiometric=allow_nonstoichiometric)
            if result is not None:
                accepted.append(result)

    print(f"│  Accepted slabs    : {len(accepted)}")
    print("└" + "─" * 68 + "┘")
    return accepted


# =============================================================================
#  4.  REPORTING
# =============================================================================

def _get_cation_hydration(slab):
    """
    For each unique cation species at the top and bottom surfaces,
    count the number of water-O neighbours within 2.8 Å.

    A "water O" is an oxygen that has exactly two H neighbours within
    1.2 Å (i.e. it belongs to H₂O, not OH⁻).  Oxygens with 0 or 1 H
    neighbour are ignored (framework oxygen or hydroxide).
    """
    cats = [s for s in slab if s.specie.element.symbol not in COVALENT_ELEMENTS]
    if not cats:
        return []
    cats.sort(key=lambda s: s.frac_coords[2])
    outermost = [cats[0], cats[-1]]

    # Pre-build a dict of O-site index → number of H within 1.2 Å, so we
    # avoid calling get_neighbors inside the outer loop.
    h_count = {
        i: sum(1 for n in nns if n.specie.symbol == 'H')
        for i, nns in enumerate(slab.get_all_neighbors(r=1.2))
        if slab[i].specie.symbol == 'O'
    }

    rows = []
    for c in outermost:
        waters = sum(
            1 for n in slab.get_neighbors(c, 2.8)
            if n.specie.symbol == 'O' and h_count.get(n.index, 0) == 2
        )
        rows.append((c.specie, c.frac_coords[2], waters))
    return rows


def print_multiline(label, value, width=58):
    """Print key: value with continuation lines if value is long."""
    line = f"  │  {label}: {value}"
    if len(line) <= width + 8:
        print(line)
        return
    # split value across lines
    print(f"  │  {label}:")
    for chunk in str(value).split(';'):
        chunk = chunk.strip()
        if chunk:
            print(f"  │    {chunk}")


def report_slab_detail(slab, index, bulk_len, target_thickness, dipole_tol,
                       stability=None, stability_rank=None, thickness_range=None,
                       net_charge=None):
    """
    Print a richly formatted report for one accepted slab.

    Parameters
    ----------
    stability : (coord_sat, n_undercoord, broken_per_area) or None
    stability_rank : int (1 = most stable) or None
    """
    an = ChemicalAnalyzer(slab)
    an_centers = len(an.centers)
    mol_counts = {k: len(v) for k, v in an.molecules.items()}
    hydration = _get_cation_hydration(slab)

    miller = slab.miller_index
    hkl_str = "".join(str(x) for x in miller)
    thickness = float(np.ptp(slab.cart_coords[:, 2]))
    dipole_z = float(slab.dipole[2])
    is_sym = slab.is_symmetric()
    mult = len(slab) / bulk_len if bulk_len else 0

    a_vec = np.array(slab.lattice.matrix[0])
    b_vec = np.array(slab.lattice.matrix[1])
    area = float(np.linalg.norm(np.cross(a_vec, b_vec)))

    is_ns = net_charge is not None and abs(net_charge) > 0.01
    stem = f"slab_{hkl_str}_{index}_NS" if is_ns else f"slab_{hkl_str}_{index}"

    print(f"\n  ┌─ SLAB #{index}  ({hkl_str}) ".ljust(62, '─') + "┐")
    print(f"  │  Formula              : {slab.formula}")
    print(f"  │  Miller index         : ({', '.join(str(x) for x in miller)})")
    print(f"  │  File stem            : {stem}")
    print(f"  ├─ Geometry ────────────────────────────────────────────────────")
    if thickness_range is not None:
        target_str = f"range = [{thickness_range[0]:.1f}, {thickness_range[1]:.1f}]"
    elif target_thickness is not None:
        target_str = f"target = {target_thickness:.1f}"
    else:
        target_str = "target = —"
    print(f"  │  Thickness            : {thickness:.3f}  Å  ({target_str})")
    print(f"  │  Surface area         : {area:.2f}  Å²")
    print(f"  │  Multiplicity         : {mult:.2f} × bulk")
    print(f"  │  Symmetry             : {'Symmetric' if is_sym else 'Asymmetric'}")
    print(f"  ├─ Polarity ────────────────────────────────────────────────────")
    print(f"  │  z-dipole moment      : {dipole_z:+.4f}  eÅ  "
          f"(tolerance = {dipole_tol})")
    print(f"  │  Dipole filter        : {'✓ PASS' if abs(dipole_z) <= dipole_tol else '✕ FAIL'}")
    print(f"  ├─ Composition ─────────────────────────────────────────────────")
    print(f"  │  Total atoms          : {len(slab)}")
    print(f"  │  Scaffold centres     : {an_centers}")
    print(f"  │  Non-scaffold         : {len(slab) - an_centers}")

    if mol_counts:
        print(f"  │  Molecular species    :")
        for form, cnt in sorted(mol_counts.items()):
            print(f"  │    {form:<12s} : {cnt} unit{'s' if cnt != 1 else ''}")

    if hydration:
        print(f"  │  Surface cation hydration:")
        for spec, z, n in hydration:
            print(f"  │    {str(spec):<12s} at z = {z:.3f}  →  {n} H₂O neighbour{'s' if n != 1 else ''}")

    if stability is not None:
        coord_sat, n_undercoord, broken_per_area = stability
        rank_str = (f"  #{stability_rank}"
                    + (" ← most stable" if stability_rank == 1 else "")
                    if stability_rank is not None else "")
        print(f"  ├─ Surface stability{rank_str}  (lower broken bonds/Å² = more stable)")
        print(f"  │  Broken bonds / Å²   : {broken_per_area:.4f}  ← primary ranking metric")
        print(f"  │  Coord. saturation    : {coord_sat:.3f}  "
              f"(mean surface CN / bulk CN, all non-H species)")
        print(f"  │  Undercoord. atoms    : {n_undercoord}")

    if is_ns:
        print(f"  ├─ *** NON-STOICHIOMETRIC ***  ──────────────────────────────────")
        print(f"  │  Net charge           : {net_charge:+.2f} e  (slab is not charge-neutral)")
        print(f"  │  This surface is intended for reactivity studies.")
        print(f"  │  Add H atoms manually to neutralise and passivate dangling bonds.")

    print(f"  └" + "─" * 60 + "┘")
    return f"{stem}.cif"


# =============================================================================
#  5.  DIAGNOSTIC TABLE  (when no surfaces pass the filter)
# =============================================================================

def print_no_surfaces_diagnostic():
    """Suggest parameter tuning when no slabs are accepted."""
    print()
    print("  ┌" + "─" * 60 + "┐")
    print("  │  No surfaces passed the dipole filter.  Possible fixes:   │")
    print("  │  1. Increase dipole_threshold (currently set in config).  │")
    print("  │  2. Increase target_thickness (thicker slabs are better    │")
    print("  │     at balancing dipole).                                  │")
    print("  │  3. Try different Miller indices — some cuts are inherently │")
    print("  │     polar.  For sheet silicates, (001) is usually the      │")
    print("  │     best choice.                                           │")
    print("  └" + "─" * 60 + "┘")
