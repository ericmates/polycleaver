"""
polycleaver: non-polar, stoichiometric surface slabs of complex minerals.

polycleaver classifies the atoms of a bulk structure into a polyanion scaffold,
counter-cations and discrete molecules (H2O, OH-, NH3, NH4+ ...), cuts it with
pymatgen's SlabGenerator, rebuilds each slab by placing cations and molecules
combinatorially to cancel the z-dipole, and ranks the resulting terminations by
broken bonds per unit area.  The 0.x series is described in Mates-Torres &
Rimola, J. Appl. Cryst. 57, 503-508 (2024); version 1.0 is a rewrite.

Typical programmatic use::

    import polycleaver

    summary = polycleaver.run({"filename": "forsterite.cif",
                               "single_hkl": [1, 1, 0],
                               "thickness_range": [12, 20]}, outdir="out")
    best = summary["slabs"][0]          # rank 1: fewest broken bonds per area
    print(best["path"], best["thickness_A"], best["dipole_z_eA"])

or, one step at a time::

    from pymatgen.core import Structure
    from polycleaver import ChemicalAnalyzer, ensure_oxidation_states, generate_slabs_for_hkl

    bulk = ensure_oxidation_states(Structure.from_file("forsterite.cif"))
    bulk_an = ChemicalAnalyzer(bulk)
    slabs = generate_slabs_for_hkl(bulk, bulk_an, (1, 1, 0),
                                   target_thickness=15.0, vacuum=15.0, dipole_tol=0.3)
    for i, s in enumerate(slabs):
        s.to(filename=f"slab_110_{i}.cif")
"""

__version__ = "1.0.0"

from .slabgen import (                       # noqa: E402,F401
    COVALENT_ELEMENTS,
    OXI_MAP,
    ChemicalAnalyzer,
    SlabProcessor,
    ensure_oxidation_states,
    order_if_disordered,
    structure_to_slab,
    load_config,
    generate_slabs_for_hkl,
    score_surface_stability,
    report_slab_detail,
    print_no_surfaces_diagnostic,
    _compute_bulk_cn_all,
    _deduplicate_slabs,
    _run_phases_on_scaffold,
    _combinatorial_search,
    _nonstoich_charges,
)
from .cli import run, analyse                # noqa: E402,F401
from .summary import bulk_record, slab_record  # noqa: E402,F401

__all__ = [
    "__version__",
    "COVALENT_ELEMENTS",
    "OXI_MAP",
    "ChemicalAnalyzer",
    "SlabProcessor",
    "ensure_oxidation_states",
    "order_if_disordered",
    "structure_to_slab",
    "load_config",
    "generate_slabs_for_hkl",
    "score_surface_stability",
    "report_slab_detail",
    "print_no_surfaces_diagnostic",
    "run",
    "analyse",
    "bulk_record",
    "slab_record",
]
