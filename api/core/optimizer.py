"""
MACE-MP structure optimization and surface energy / Wulff construction.

All heavy imports (mace, ase) are deferred to function bodies so the rest of
the API starts up without them installed.
"""

import traceback
import numpy as np

from pymatgen.core import Structure
from pymatgen.io.ase import AseAtomsAdaptor

from polycleaver2 import ensure_oxidation_states


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _get_calculator():
    from mace.calculators import mace_mp
    return mace_mp(model="medium", dispersion=False,
                   default_dtype="float32", device="cpu")


def _phase_label(phase, total):
    if phase == 'bulk':
        return 'Optimizing bulk…'
    if phase == 'done':
        return 'Done'
    try:
        idx = int(phase.split('_')[1])
        return f'Optimizing slab {idx + 1}/{total}…'
    except Exception:
        return phase


# ---------------------------------------------------------------------------
# Public optimizers
# ---------------------------------------------------------------------------

def optimize_bulk(structure, fmax=0.05, steps=300):
    """Relax bulk cell + atoms with ExpCellFilter. Returns (optimized_structure, energy_eV)."""
    from ase.optimize import LBFGS
    from ase.constraints import ExpCellFilter

    calc = _get_calculator()
    adaptor = AseAtomsAdaptor()
    atoms = adaptor.get_atoms(structure)
    atoms.calc = calc

    ecf = ExpCellFilter(atoms)
    opt = LBFGS(ecf, logfile=None)
    opt.run(fmax=fmax, steps=steps)

    energy = float(atoms.get_potential_energy())
    optimized = adaptor.get_structure(atoms)
    return optimized, energy


def optimize_slab(slab, fmax=0.05, steps=300):
    """Relax slab atoms at fixed cell. Returns energy_eV."""
    from ase.optimize import LBFGS

    calc = _get_calculator()
    adaptor = AseAtomsAdaptor()
    atoms = adaptor.get_atoms(slab)
    atoms.calc = calc

    opt = LBFGS(atoms, logfile=None)
    opt.run(fmax=fmax, steps=steps)

    return float(atoms.get_potential_energy())


def surface_energy(e_slab, n_slab, e_bulk, n_bulk, area_ang2):
    """
    E_surf = (E_slab - (N_slab/N_bulk)*E_bulk) / (2*A)
    Returns (eV/Å², J/m²).
    """
    ev_per_ang2 = (e_slab - (n_slab / n_bulk) * e_bulk) / (2.0 * area_ang2)
    j_per_m2 = ev_per_ang2 * 16.0218
    return float(ev_per_ang2), float(j_per_m2)


# ---------------------------------------------------------------------------
# Wulff shape serialisation
# ---------------------------------------------------------------------------

def serialize_wulff(wulff_shape):
    """Convert WulffShape → JSON-serialisable dict for the frontend."""
    faces = []
    for facet in wulff_shape.facets:
        verts = facet.vertices
        if len(verts) < 3:
            continue
        faces.append({
            'hkl':      list(facet.miller_index),
            'energy':   float(facet.e_surf),
            'vertices': [[float(c) for c in v] for v in verts],
        })

    miller_energy = {
        f"({','.join(str(x) for x in hkl)})": float(e)
        for hkl, e in wulff_shape.miller_energy_dict.items()
    }

    return {
        'faces':                    faces,
        'miller_energy':            miller_energy,
        'anisotropy':               float(wulff_shape.anisotropy),
        'weighted_surface_energy':  float(wulff_shape.weighted_surface_energy),
    }


# ---------------------------------------------------------------------------
# Background job
# ---------------------------------------------------------------------------

def run_optimization_job(opt_job_id, job_id, session_id, sessions, jobs, opt_jobs):
    """
    Runs in a background thread:
      1. Optimize bulk structure
      2. Optimize each slab from jobs[job_id]['raw_slabs']
      3. Compute surface energies
      4. Build Wulff construction
    Writes progress into opt_jobs[opt_job_id].
    """
    try:
        from pymatgen.analysis.wulff import WulffShape

        opt_jobs[opt_job_id]['status'] = 'running'
        cif_path = sessions[session_id]['cif_path']

        bulk = Structure.from_file(cif_path)
        bulk = ensure_oxidation_states(bulk)
        n_bulk = len(bulk)

        raw_slabs = jobs[job_id].get('raw_slabs', [])
        n_slabs   = len(raw_slabs)
        opt_jobs[opt_job_id]['total'] = n_slabs

        # --- Bulk optimisation ---
        opt_jobs[opt_job_id]['phase'] = 'bulk'
        opt_jobs[opt_job_id]['phase_label'] = _phase_label('bulk', n_slabs)
        _, e_bulk = optimize_bulk(bulk)
        opt_jobs[opt_job_id]['e_bulk'] = e_bulk

        # --- Slab optimisations ---
        miller_best = {}   # hkl_tuple -> minimum j_m2 found so far
        slab_results = []

        for i, slab in enumerate(raw_slabs):
            phase = f'slab_{i}'
            opt_jobs[opt_job_id]['phase']       = phase
            opt_jobs[opt_job_id]['phase_label'] = _phase_label(phase, n_slabs)
            opt_jobs[opt_job_id]['current']     = i

            a_vec    = np.array(slab.lattice.matrix[0])
            b_vec    = np.array(slab.lattice.matrix[1])
            area     = float(np.linalg.norm(np.cross(a_vec, b_vec)))
            n_slab   = len(slab)

            e_slab   = optimize_slab(slab)
            ev_a2, j_m2 = surface_energy(e_slab, n_slab, e_bulk, n_bulk, area)

            hkl      = tuple(int(x) for x in slab.miller_index)
            if hkl not in miller_best or j_m2 < miller_best[hkl]:
                miller_best[hkl] = j_m2

            slab_results.append({
                'slab_index':       i,
                'hkl':              list(hkl),
                'e_surf_ev_ang2':   ev_a2,
                'e_surf_j_m2':      j_m2,
                'area':             area,
            })

            opt_jobs[opt_job_id]['completed']    = i + 1
            opt_jobs[opt_job_id]['slab_results'] = slab_results   # live update

        # --- Wulff construction ---
        wulff_data = None
        wulff_error = None
        try:
            me_dict = {hkl: max(e, 1e-6) for hkl, e in miller_best.items()}
            if me_dict:
                ws = WulffShape(bulk.lattice,
                                list(me_dict.keys()),
                                list(me_dict.values()))
                wulff_data = serialize_wulff(ws)
        except Exception as exc:
            wulff_error = str(exc)

        opt_jobs[opt_job_id].update({
            'status':       'done',
            'phase':        'done',
            'phase_label':  'Done',
            'slab_results': slab_results,
            'wulff_data':   wulff_data,
            'wulff_error':  wulff_error,
        })

    except Exception as exc:
        opt_jobs[opt_job_id]['status']    = 'error'
        opt_jobs[opt_job_id]['error']     = str(exc)
        opt_jobs[opt_job_id]['traceback'] = traceback.format_exc()
