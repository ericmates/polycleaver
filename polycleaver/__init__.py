"""
polycleaver — non-polar, stoichiometric surface slabs of complex minerals.

Successor to PolyCleaver (Mates-Torres & Rimola, J. Appl. Cryst. 57, 503, 2024).
Handles isolated-polyhedra, chain, sheet and framework polyanion scaffolds,
counter-cations, and discrete molecular species (H2O, OH-, NH3, NH4+ …).

Typical programmatic use::

    from pymatgen.core import Structure
    from polycleaver import ChemicalAnalyzer, ensure_oxidation_states, generate_slabs_for_hkl

    bulk = ensure_oxidation_states(Structure.from_file("forsterite.cif"))
    bulk_an = ChemicalAnalyzer(bulk)
    slabs = generate_slabs_for_hkl(bulk, bulk_an, (1, 1, 0),
                                   target_thickness=15.0, vacuum=15.0, dipole_tol=0.3)
    for i, s in enumerate(slabs):
        s.to(filename=f"slab_110_{i}.cif")
"""

__version__ = "2.0.0a1"

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
]
