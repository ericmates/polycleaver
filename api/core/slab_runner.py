import io
import re
import sys
import traceback

import numpy as np

import sys as _sys
import os as _os
_sys.path.insert(0, _os.path.join(_os.path.dirname(__file__), '..', '..'))

from pymatgen.core import Structure
from pymatgen.core.surface import get_symmetrically_distinct_miller_indices
from polycleaver import ChemicalAnalyzer, ensure_oxidation_states, generate_slabs_for_hkl
from api.core.structure_serializer import serialize_slab


class _ProgressCapture(io.StringIO):
    """Replaces sys.stdout during generation; parses progress print lines."""

    # Matches: "│  Unique cuts found  : 13"
    _CUTS_RE   = re.compile(r'Unique cuts found\s*:\s*(\d+)')
    # Matches: "│  ─── Cut 3 ───…"
    _CUT_RE    = re.compile(r'─── Cut (\d+) ')
    # Matches: "│     ✓ SLAB ACCEPTED.  Final dipole = …"
    _ACCEPT_RE = re.compile(r'SLAB ACCEPTED')

    def __init__(self, job_id, hkl, jobs):
        super().__init__()
        self.job_id      = job_id
        self.hkl         = list(hkl)
        self.jobs        = jobs
        self._total_cuts = None
        self._cut_done   = 0

    def _push(self, success):
        job   = self.jobs[self.job_id]
        total = self._total_cuts or 1
        pct   = round(100 * self._cut_done / total, 1)
        job['progress'].append({
            'hkl':        self.hkl,
            'percent':    pct,
            'cut':        self._cut_done,
            'total_cuts': self._total_cuts,
            'success':    success,
        })

    def write(self, s):
        m = self._CUTS_RE.search(s)
        if m:
            self._total_cuts = int(m.group(1))
            self._cut_done   = 0
            return

        m = self._CUT_RE.search(s)
        if m:
            self._cut_done = int(m.group(1)) + 1
            self._push(False)
            return

        if self._ACCEPT_RE.search(s):
            self._push(True)

    def flush(self):
        pass


def run_job(job_id, session_id, params, sessions, jobs):
    """Run slab generation in a background thread and update jobs[job_id]."""
    try:
        session  = sessions[session_id]
        cif_path = session['cif_path']
        bulk_len = session.get('bulk_len')

        # Re-load structure
        bulk = Structure.from_file(cif_path)
        bulk = ensure_oxidation_states(bulk)

        bulk_an = ChemicalAnalyzer(bulk)
        bulk_len = len(bulk)

        # Resolve miller indices
        mi    = params['miller_input']
        mode  = mi['mode']
        value = mi['value']
        if mode == 'max':
            target_indices = get_symmetrically_distinct_miller_indices(bulk, int(value))
        elif mode == 'single':
            target_indices = [tuple(int(x) for x in value)]
        else:  # 'list'
            target_indices = [tuple(int(x) for x in v) for v in value]

        jobs[job_id]['total_indices']     = len(target_indices)
        jobs[job_id]['completed_indices'] = 0
        jobs[job_id]['results']           = []
        jobs[job_id]['raw_slabs']         = []

        for idx, hkl in enumerate(target_indices):
            jobs[job_id]['current_hkl'] = list(hkl)

            old_stdout = sys.stdout
            sys.stdout = _ProgressCapture(job_id, hkl, jobs)
            try:
                slabs = generate_slabs_for_hkl(
                    bulk, bulk_an, hkl,
                    params['thickness'],
                    params['vacuum'],
                    params['dipole_tol'],
                    thickness_range=params.get('thickness_range'),
                    allow_nonstoichiometric=bool(params.get('allow_nonstoichiometric', False)),
                )
            finally:
                sys.stdout = old_stdout

            jobs[job_id]['completed_indices'] = idx + 1

            for slab in slabs:
                jobs[job_id]['raw_slabs'].append(slab)
                jobs[job_id]['results'].append(serialize_slab(slab, bulk_len=bulk_len))

        jobs[job_id]['status'] = 'done'

    except Exception as e:
        jobs[job_id]['status']    = 'error'
        jobs[job_id]['error']     = str(e)
        jobs[job_id]['traceback'] = traceback.format_exc()
