#!/usr/bin/env python3
"""
Report the z-dipole (polarity) of one or more slab files.

    python scripts/analyze_polarity.py slab_110_0.cif [more.cif …]

Oxidation states are assigned with the same logic as the generator
(pymatgen guess → OXI_MAP fallback), so the number printed here is the same
quantity the generator filters on (|dipole_z| in eÅ).
"""

import sys
import numpy as np
from pymatgen.core import Structure
from pymatgen.core.surface import Slab

from polycleaver import ensure_oxidation_states


def slab_dipole(path):
    structure = ensure_oxidation_states(Structure.from_file(path))
    slab = Slab(
        lattice=structure.lattice,
        species=structure.species,
        coords=structure.frac_coords,
        miller_index=(0, 0, 1),
        oriented_unit_cell=structure,
        shift=0,
        scale_factor=np.eye(3),
    )
    return slab.dipole, float(np.ptp(slab.cart_coords[:, 2]))


def main(paths):
    if not paths:
        print(__doc__)
        sys.exit(1)
    for path in paths:
        dip, thick = slab_dipole(path)
        pz = float(dip[2])
        verdict = "non-polar (Tasker type 1/2)" if abs(pz) < 1e-2 else "polar (Tasker type 3)"
        print(f"{path}")
        print(f"  dipole vector (eÅ) : {dip[0]:+.4f} {dip[1]:+.4f} {dip[2]:+.4f}")
        print(f"  |dipole_z|         : {abs(pz):.4f} eÅ  →  {verdict}")
        print(f"  thickness          : {thick:.2f} Å")


if __name__ == "__main__":
    main(sys.argv[1:])
