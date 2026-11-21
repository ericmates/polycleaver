"""
Machine-readable records of a polycleaver run.

The command-line interface writes these records to ``polycleaver_summary.json``
next to the CIF files (and prints them with ``--json``); the MCP server returns
them as tool results.  Every quantity is a plain JSON type, so that a script,
a workflow manager or an LLM agent can select slabs without parsing the
human-readable report.
"""

import platform
from collections import Counter
from importlib import metadata

import numpy as np

from .slabgen import ChemicalAnalyzer, _get_cation_hydration

#: Version of the layout of the records below.  It changes only when a key is
#: removed or changes meaning; new keys may be added in minor releases.
SUMMARY_SCHEMA_VERSION = 1


def environment_versions():
    """Versions of Python and of the libraries that determine the slab cuts."""
    versions = {"python": platform.python_version()}
    for pkg in ("polycleaver", "pymatgen", "ase", "numpy", "spglib"):
        try:
            versions[pkg] = metadata.version(pkg)
        except metadata.PackageNotFoundError:
            versions[pkg] = None
    return versions


def _counts(structure, indices):
    return dict(sorted(Counter(structure[i].specie.element.symbol for i in indices).items()))


def bulk_record(bulk, analyzer, bulk_cn=None):
    """Return a JSON-serialisable description of the bulk chemical analysis."""
    try:
        space_group = bulk.get_space_group_info()[0]
    except Exception:                                 # pragma: no cover - spglib failure
        space_group = None
    net_charge = float(sum(site.specie.oxi_state for site in bulk))
    centre_oxi = {}
    for i in analyzer.centers:
        centre_oxi.setdefault(bulk[i].specie.element.symbol, float(bulk[i].specie.oxi_state))
    cation_oxi = {}
    for i in analyzer.cations:
        cation_oxi.setdefault(bulk[i].specie.element.symbol, float(bulk[i].specie.oxi_state))
    if not analyzer.centers:
        connectivity = "none"
    elif analyzer.has_shared_oxygens():
        connectivity = "framework_or_sheet"
    else:
        connectivity = "isolated_polyhedra"
    return {
        "formula": bulk.formula,
        "reduced_formula": bulk.composition.reduced_formula,
        "n_atoms": len(bulk),
        "space_group": space_group,
        "lattice": {
            "a": bulk.lattice.a, "b": bulk.lattice.b, "c": bulk.lattice.c,
            "alpha": bulk.lattice.alpha, "beta": bulk.lattice.beta,
            "gamma": bulk.lattice.gamma, "volume": bulk.lattice.volume,
        },
        "net_charge_e": net_charge,
        "charge_neutral": abs(net_charge) <= 0.01,
        "scaffold_connectivity": connectivity,
        "scaffold_centres": _counts(bulk, analyzer.centers),
        "scaffold_centre_oxidation_states": centre_oxi,
        "scaffold_anions": _counts(bulk, analyzer.anions),
        "cations": _counts(bulk, analyzer.cations),
        "cation_oxidation_states": cation_oxi,
        "molecules": {k: len(v) for k, v in sorted(analyzer.molecules.items())},
        "cation_hydration_target": dict(sorted(analyzer.cation_hydration_target.items())),
        "bulk_coordination_numbers": (
            {k: float(v) for k, v in sorted(bulk_cn.items())} if bulk_cn else None
        ),
    }


def slab_record(slab, filename, bulk_len, stability=None, rank=None,
                net_charge=None, dipole_tol=None):
    """Return a JSON-serialisable description of one accepted slab."""
    an = ChemicalAnalyzer(slab)
    thickness = float(np.ptp(slab.cart_coords[:, 2]))
    a_vec = np.array(slab.lattice.matrix[0])
    b_vec = np.array(slab.lattice.matrix[1])
    area = float(np.linalg.norm(np.cross(a_vec, b_vec)))
    dipole_z = float(slab.dipole[2])
    if net_charge is None:
        net_charge = float(sum(site.specie.oxi_state for site in slab))
    rec = {
        "rank": rank,
        "file": filename,
        "miller_index": [int(x) for x in slab.miller_index],
        "formula": slab.formula,
        "reduced_formula": slab.composition.reduced_formula,
        "n_atoms": len(slab),
        "multiplicity_vs_bulk": len(slab) / bulk_len if bulk_len else None,
        "thickness_A": thickness,
        "surface_area_A2": area,
        "symmetric": bool(slab.is_symmetric()),
        "dipole_z_eA": dipole_z,
        "dipole_within_tolerance": (abs(dipole_z) <= dipole_tol) if dipole_tol is not None else None,
        "net_charge_e": float(net_charge),
        "stoichiometric": abs(net_charge) <= 0.01,
        "scaffold_centres": _counts(slab, an.centers),
        "molecules": {k: len(v) for k, v in sorted(an.molecules.items())},
        "surface_cation_hydration": [
            {"species": str(spec), "z_frac": float(z), "n_H2O": int(n)}
            for spec, z, n in _get_cation_hydration(slab)
        ],
        "lattice": {
            "a": slab.lattice.a, "b": slab.lattice.b, "c": slab.lattice.c,
            "alpha": slab.lattice.alpha, "beta": slab.lattice.beta,
            "gamma": slab.lattice.gamma,
        },
    }
    if stability is not None:
        coord_sat, n_undercoord, broken_per_area = stability
        rec.update({
            "broken_bonds_per_A2": float(broken_per_area),
            "coordination_saturation": float(coord_sat),
            "n_undercoordinated": int(n_undercoord),
        })
    return rec


NO_SLAB_HINTS = [
    "Raise dipole_threshold moderately (e.g. 0.5-1.0 eÅ) and report the value used.",
    "Increase target_thickness, or scan a wider thickness_range; some cuts need "
    "at least two crystallographic periods along the surface normal.",
    "Try other Miller indices; some cuts are polar by construction.",
    "Check that the bulk cell is ordered and charge-neutral (bulk.charge_neutral).",
]
