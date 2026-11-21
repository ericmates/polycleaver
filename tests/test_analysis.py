"""
Fast tests: bulk chemical analysis on the bundled example structures.

These exercise the classification logic (scaffold centres, anions, cations,
molecules, framework detection, hydration targets) without generating slabs.
Run with:  pytest tests/
"""

from pathlib import Path

import pytest
from pymatgen.core import Structure

from polycleaver2 import ChemicalAnalyzer, ensure_oxidation_states, _compute_bulk_cn_all

STRUCTURES = Path(__file__).resolve().parent.parent / "examples" / "structures"


def load(name):
    return ensure_oxidation_states(Structure.from_file(STRUCTURES / name))


def test_forsterite_isolated_tetrahedra():
    bulk = load("forsterite.cif")
    an = ChemicalAnalyzer(bulk)
    syms = lambda idxs: sorted({bulk[i].specie.element.symbol for i in idxs})
    assert syms(an.centers) == ["Si"]
    assert len(an.centers) == 4
    assert len(an.anions) == 16
    assert syms(an.cations) == ["Mg"]
    assert len(an.cations) == 8
    assert not an.molecules
    assert not an.has_shared_oxygens()          # isolated SiO4
    assert an.cn_by_species["Si"] == 4
    cn = _compute_bulk_cn_all(bulk)
    assert cn["Si"] == pytest.approx(4.0)
    assert cn["Mg"] == pytest.approx(6.0)
    assert cn["O"] == pytest.approx(4.0)


def test_hydroxyapatite_hydroxide_is_molecular():
    bulk = load("hydroxyapatite.cif")
    an = ChemicalAnalyzer(bulk)
    assert {bulk[i].specie.element.symbol for i in an.centers} == {"P"}
    assert len(an.centers) == 6
    assert {bulk[i].specie.element.symbol for i in an.cations} == {"Ca"}
    assert len(an.cations) == 10
    assert set(an.molecules) == {"HO"}
    assert len(an.molecules["HO"]) == 2
    assert not an.has_shared_oxygens()


def test_struvite_hydration_and_ammonium():
    bulk = load("struvite.cif")                 # MgNH4PO4·6H2O, Z = 2
    an = ChemicalAnalyzer(bulk)
    assert {bulk[i].specie.element.symbol for i in an.centers} == {"P"}
    assert len(an.centers) == 2
    assert {bulk[i].specie.element.symbol for i in an.cations} == {"Mg"}
    assert set(an.molecules) == {"H2O", "H4N"}
    assert len(an.molecules["H2O"]) == 12
    assert len(an.molecules["H4N"]) == 2
    assert an.cation_hydration_target["Mg"] == 6   # hexaaqua Mg2+


def test_kaolinite_hydroxyls_travel_with_scaffold():
    """Structural OH: H must be attached to the Al/Si unit, never orphaned."""
    bulk = load("kaolinite.cif")                # Al2Si2O5(OH)4, COD 9009234
    an = ChemicalAnalyzer(bulk)
    assert {bulk[i].specie.element.symbol for i in an.centers} == {"Al", "Si"}
    assert not an.cations and not an.molecules
    assert an.has_shared_oxygens()              # 1:1 T–O layer
    h_in_units = {i for unit in an.covalent_map.values() for i in unit
                  if bulk[i].specie.element.symbol == "H"}
    assert len(h_in_units) == 8                 # all 8 H belong to some unit
    assert abs(sum(s.specie.oxi_state for s in bulk)) < 1e-9


def test_diopside_chain_is_framework():
    bulk = load("diopside.cif")
    an = ChemicalAnalyzer(bulk)
    assert {bulk[i].specie.element.symbol for i in an.centers} == {"Si"}
    assert {bulk[i].specie.element.symbol for i in an.cations} == {"Ca", "Mg"}
    assert an.has_shared_oxygens()              # single-chain SiO3


def test_oxidation_states_assigned():
    bulk = load("forsterite.cif")
    oxi = {s.specie.element.symbol: s.specie.oxi_state for s in bulk}
    assert oxi == {"Mg": 2, "Si": 4, "O": -2}
    assert abs(sum(s.specie.oxi_state for s in bulk)) < 1e-9
