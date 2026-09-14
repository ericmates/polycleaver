import math
import numpy as np

CPK_COLORS = {
    'H':  '#FFFFFF', 'C':  '#909090', 'N':  '#3050F8', 'O':  '#FF0D0D',
    'F':  '#90E050', 'P':  '#FF8000', 'S':  '#FFFF30', 'Cl': '#1FF01F',
    'Ca': '#3DFF00', 'Na': '#AB5CF2', 'K':  '#8F40D4', 'Mg': '#8AFF00',
    'Al': '#BFA6A6', 'Si': '#F0C8A0', 'Fe': '#E06633', 'Mn': '#9C7AC7',
    'Zn': '#7D80B0', 'Cu': '#C88033', 'B':  '#FFB5B5', 'Br': '#A62929',
    'I':  '#940094', 'As': '#BD80E3', 'Se': '#FFA100',
}

VDW_RADII = {
    'H':  1.20, 'C':  1.70, 'N':  1.55, 'O':  1.52, 'F':  1.47,
    'P':  1.80, 'S':  1.80, 'Cl': 1.75, 'Ca': 2.31, 'Na': 2.27,
    'K':  2.75, 'Mg': 1.73, 'Al': 1.84, 'Si': 2.10, 'Fe': 2.04,
    'Mn': 2.05, 'Zn': 2.10, 'Cu': 1.96, 'B':  1.92, 'Br': 1.85,
    'I':  1.98, 'As': 1.85, 'Se': 1.90,
}

# Per-element display radius overrides (Å, already scaled for aesthetics)
DISPLAY_RADII = {
    'H': 0.14,
    'O': 0.34,
    'N': 0.36,
}
ATOM_SCALE = 0.4


def _display_radius(sym):
    if sym in DISPLAY_RADII:
        return DISPLAY_RADII[sym]
    return VDW_RADII.get(sym, 2.0) * ATOM_SCALE


def _safe(v):
    if isinstance(v, np.ndarray):
        return v.tolist()
    if isinstance(v, np.generic):   # catches all numpy scalars (integer, floating, longlong, …)
        pv = v.item()
        if isinstance(pv, float) and (math.isinf(pv) or math.isnan(pv)):
            return None
        return pv
    if isinstance(v, float) and (math.isinf(v) or math.isnan(v)):
        return None
    return v


def _chemistry(structure, analyzer, display_cart=None, include_polyhedra=False):
    """
    Extract from ChemicalAnalyzer:
      - water_atom_indices  : site indices belonging to H2O molecules
      - water_bonds         : O-H bond pairs within water  {a, b, elem_a, elem_b}
      - chemical_bonds      : scaffold center-O bonds + non-water molecule bonds
      - polyhedra           : coordination polyhedra with PBC-corrected vertices

    display_cart  : optional (N,3) array of cartesian coords to use for atom
                    positions (wrapped coords for slab display). Original
                    fractional coords are always used for PBC/distance logic.
    include_polyhedra : only True for slabs (bulk shows empty list)
    """
    if analyzer is None:
        return [], [], [], []

    disp = display_cart if display_cart is not None else np.array([s.coords for s in structure])
    water_atom_indices, water_bonds, chemical_bonds, polyhedra = [], [], [], []

    # ---- Molecular clusters (water and other) ----
    for formula, mol_list in analyzer.molecules.items():
        is_water = 'H' in formula and 'O' in formula
        for cluster in mol_list:
            if is_water:
                water_atom_indices.extend(cluster)
            for ki in range(len(cluster)):
                i = cluster[ki]
                ei = structure[i].specie.element.symbol
                for kj in range(ki + 1, len(cluster)):
                    j = cluster[kj]
                    ej = structure[j].specie.element.symbol
                    # Use minimum-image distance — correct across PBC
                    delta_frac = structure[j].frac_coords - structure[i].frac_coords
                    delta_frac -= np.round(delta_frac)
                    dc = structure.lattice.get_cartesian_coords(delta_frac)
                    dist = float(np.linalg.norm(dc))
                    is_xh = 'H' in (ei, ej)
                    thresh = 1.25 if is_xh else 2.3
                    if dist < thresh:
                        bd = {'a': int(i), 'b': int(j), 'elem_a': ei, 'elem_b': ej}
                        if is_water and is_xh:
                            # Store PBC-corrected positions so water bonds never
                            # span across the cell (hydrogen placed relative to O)
                            o_loc = i if ei == 'O' else j
                            h_loc = j if ei == 'O' else i
                            o_pos = disp[o_loc]
                            df = structure[h_loc].frac_coords - structure[o_loc].frac_coords
                            df -= np.round(df)
                            h_pos = o_pos + structure.lattice.get_cartesian_coords(df)
                            bd['pos_a'] = [_safe(c) for c in (o_pos if i == o_loc else h_pos)]
                            bd['pos_b'] = [_safe(c) for c in (h_pos if i == o_loc else o_pos)]
                        (water_bonds if is_water else chemical_bonds).append(bd)

    # ---- Scaffold bonds + coordination polyhedra ----
    for c_idx, cluster in analyzer.covalent_map.items():
        center_sym = structure[c_idx].specie.element.symbol
        oxy_idx = [x for x in cluster if structure[x].specie.element.symbol == 'O']

        for ox in oxy_idx:
            chemical_bonds.append({'a': int(c_idx), 'b': int(ox),
                                    'elem_a': center_sym, 'elem_b': 'O'})

        if include_polyhedra and len(oxy_idx) >= 3:
            # Use minimum-image convention to bring oxygens near the (wrapped) center
            center_frac = structure[c_idx].frac_coords
            center_cart = disp[c_idx]
            vertices = []
            for ox in oxy_idx:
                df = structure[ox].frac_coords - center_frac
                df -= np.round(df)            # minimum image in fractional space
                dc = structure.lattice.get_cartesian_coords(df)
                vertices.append([_safe(c) for c in (center_cart + dc)])
            polyhedra.append({
                'color': CPK_COLORS.get(center_sym, '#FF69B4'),
                'center': [_safe(c) for c in center_cart],
                'vertices': vertices,
            })

    return water_atom_indices, water_bonds, chemical_bonds, polyhedra


def _lattice_dict(lat):
    return {
        'matrix': [row.tolist() for row in lat.matrix],
        'a': _safe(lat.a), 'b': _safe(lat.b), 'c': _safe(lat.c),
        'alpha': _safe(lat.alpha), 'beta': _safe(lat.beta), 'gamma': _safe(lat.gamma),
    }


def serialize_structure(structure, analyzer=None):
    """Serialise a bulk Structure. Bonds = scaffold + molecules only (no all-pairs)."""
    atoms = []
    for site in structure:
        sym = site.specie.element.symbol
        x, y, z = [_safe(c) for c in site.coords]
        atoms.append({
            'element': sym,
            'x': x, 'y': y, 'z': z,
            'radius': _display_radius(sym),
            'color': CPK_COLORS.get(sym, '#FF69B4'),
        })

    water_indices, water_bonds, chem_bonds, _ = _chemistry(
        structure, analyzer, include_polyhedra=False)

    return {
        'atoms': atoms,
        'lattice': _lattice_dict(structure.lattice),
        'bonds': chem_bonds,
        'water_atom_indices': water_indices,
        'water_bonds': water_bonds,
        'polyhedra': [],     # never shown for bulk
    }


def serialize_slab(slab, bulk_len=None):
    """Serialise a Slab with PBC-wrapped display coords and coordination polyhedra."""
    from polycleaver2 import ChemicalAnalyzer

    an = ChemicalAnalyzer(slab)

    # Wrap x, y fractional coords to [0, 1) — fixes atoms displaced to
    # neighbouring cells in the slab plane. Z is NOT wrapped (vacuum direction).
    orig_frac = slab.frac_coords          # keep original for PBC/min-image logic
    disp_frac = orig_frac.copy()
    disp_frac[:, 0] %= 1.0
    disp_frac[:, 1] %= 1.0
    disp_cart = disp_frac @ slab.lattice.matrix   # (N, 3) display positions

    atoms = []
    for i, site in enumerate(slab):
        sym = site.specie.element.symbol
        x, y, z = [_safe(c) for c in disp_cart[i]]
        atoms.append({
            'element': sym,
            'x': x, 'y': y, 'z': z,
            'radius': _display_radius(sym),
            'color': CPK_COLORS.get(sym, '#FF69B4'),
        })

    water_indices, water_bonds, chem_bonds, polyhedra = _chemistry(
        slab, an, display_cart=disp_cart, include_polyhedra=True)

    # Slab geometry (use original cart_coords for thickness — wrapping doesn't change z)
    thickness   = float(np.ptp(slab.cart_coords[:, 2]))
    lat         = slab.lattice
    a_vec, b_vec = np.array(lat.matrix[0]), np.array(lat.matrix[1])
    surface_area = float(np.linalg.norm(np.cross(a_vec, b_vec)))
    dipole_z    = _safe(abs(slab.dipole[2]))
    is_sym      = bool(slab.is_symmetric())
    mult        = _safe(len(slab) / bulk_len) if bulk_len else None

    return {
        'atoms': atoms,
        'lattice': _lattice_dict(lat),
        'bonds': chem_bonds,
        'water_atom_indices': water_indices,
        'water_bonds': water_bonds,
        'polyhedra': polyhedra,
        'meta': {
            'miller_index':  list(slab.miller_index),
            'formula':       slab.formula,
            'dipole':        dipole_z,
            'is_symmetric':  is_sym,
            'thickness':     _safe(thickness),
            'surface_area':  _safe(surface_area),
            'multiplicity':  mult,
            'molecules':     {k: len(v) for k, v in an.molecules.items()},
            'num_atoms':     len(slab),
            'n_centers':     len(an.centers),
            'n_cations':     len(an.cations),
        },
    }


def serialize_analysis(bulk_an):
    return {
        'n_centers': len(bulk_an.centers),
        'n_cations': len(bulk_an.cations),
        'n_anions':  len(bulk_an.anions),
        'molecules': {k: len(v) for k, v in bulk_an.molecules.items()},
        'cation_hydration_target': dict(bulk_an.cation_hydration_target),
    }
