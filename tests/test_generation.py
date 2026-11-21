"""
Slow, end-to-end tests: generate slabs for the bundled examples and check the
results against the documented reference outcomes.

    pytest tests/ -m slow          # run only these
    pytest tests/ -m "not slow"    # skip them
"""

from pathlib import Path

import numpy as np
import pytest
from pymatgen.core import Structure

from polycleaver2 import (
    ChemicalAnalyzer,
    ensure_oxidation_states,
    generate_slabs_for_hkl,
    score_surface_stability,
    _compute_bulk_cn_all,
    _deduplicate_slabs,
)

STRUCTURES = Path(__file__).resolve().parent.parent / "examples" / "structures"

pytestmark = pytest.mark.slow


def _load(name):
    bulk = ensure_oxidation_states(Structure.from_file(STRUCTURES / name))
    return bulk, ChemicalAnalyzer(bulk)


def _net_charge(slab):
    return sum(s.specie.oxi_state for s in slab)


def test_forsterite_110_range_mode_symmetric_termination():
    bulk, an = _load("forsterite.cif")
    slabs = generate_slabs_for_hkl(bulk, an, (1, 1, 0), vacuum=15.0, dipole_tol=0.3,
                                   thickness_range=[12, 20])
    assert slabs, "no (110) forsterite termination found"
    for s in slabs:
        assert abs(s.dipole[2]) <= 0.3
        assert abs(_net_charge(s)) < 1e-6
        assert s.composition.reduced_formula == "Mg2SiO4"

    cn = _compute_bulk_cn_all(bulk)
    scored = sorted(zip([score_surface_stability(s, cn) for s in slabs], slabs),
                    key=lambda x: (x[0][2], -int(x[1].is_symmetric()), -x[0][0]))
    unique = _deduplicate_slabs([s for _, s in scored])
    # A symmetric Mg32Si16O64 slab near 15 Å is the DFT-validated termination.
    assert any(s.is_symmetric()
               and s.composition.get_el_amt_dict()["Si"] == 16
               and abs(np.ptp(s.cart_coords[:, 2]) - 15.0) < 0.5
               for s in unique)


def test_hydroxyapatite_010_reference_slab():
    bulk, an = _load("hydroxyapatite.cif")
    slabs = generate_slabs_for_hkl(bulk, an, (0, 1, 0), target_thickness=17.0,
                                   vacuum=15.0, dipole_tol=0.3)
    assert slabs
    # Ca20(PO4)12(OH)4, 88 atoms, symmetric — matches POSCAR_hydroxyapatite_010_DFT
    ok = [s for s in slabs
          if len(s) == 88 and s.is_symmetric() and abs(s.dipole[2]) < 1e-3]
    assert ok, [ (len(s), s.is_symmetric(), float(s.dipole[2])) for s in slabs ]
    for s in slabs:
        assert abs(_net_charge(s)) < 1e-6


def test_forsterite_010_isolated_polyhedra_path():
    """The bare SiO4 scaffold is polar on (010); Mg placement must cancel it."""
    bulk, an = _load("forsterite.cif")
    slabs = generate_slabs_for_hkl(bulk, an, (0, 1, 0), target_thickness=10.0,
                                   vacuum=15.0, dipole_tol=0.3)
    assert slabs
    assert all(abs(s.dipole[2]) <= 0.3 for s in slabs)
