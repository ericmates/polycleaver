import json
import math
import os
import tempfile
import traceback
from concurrent.futures import ThreadPoolExecutor
from uuid import uuid4

import numpy as np
from fastapi import FastAPI, File, HTTPException, UploadFile
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, Response
from pydantic import BaseModel
from pymatgen.core import Structure
from pymatgen.core.surface import get_symmetrically_distinct_miller_indices

from api.core.structure_serializer import (
    serialize_analysis,
    serialize_structure,
)
from api.core.slab_runner import run_job
from api.core.optimizer import run_optimization_job
from polycleaver2 import ChemicalAnalyzer, ensure_oxidation_states

# ---------------------------------------------------------------------------
# App setup
# ---------------------------------------------------------------------------

app = FastAPI(title="PolyCleaver 2 — Surface Lab API")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:5173"],
    allow_methods=["*"],
    allow_headers=["*"],
)

executor     = ThreadPoolExecutor(max_workers=1)
opt_executor = ThreadPoolExecutor(max_workers=1)

# In-memory stores (single-user tool)
sessions: dict = {}
jobs: dict     = {}
opt_jobs: dict = {}


# ---------------------------------------------------------------------------
# JSON helper
# ---------------------------------------------------------------------------

class _NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, np.generic):   # catches all numpy scalars
            pv = obj.item()
            if isinstance(pv, float) and (math.isinf(pv) or math.isnan(pv)):
                return None
            return pv
        return super().default(obj)


def _json_response(data):
    return Response(
        content=json.dumps(data, cls=_NumpyEncoder),
        media_type="application/json",
    )


# ---------------------------------------------------------------------------
# Routes
# ---------------------------------------------------------------------------

@app.post("/upload")
async def upload(file: UploadFile = File(...)):
    if not file.filename.endswith(".cif"):
        raise HTTPException(400, "Only .cif files are supported.")

    session_id = str(uuid4())
    tmp_dir = tempfile.mkdtemp()
    cif_path = os.path.join(tmp_dir, f"{session_id}.cif")

    contents = await file.read()
    with open(cif_path, "wb") as f:
        f.write(contents)

    try:
        bulk = Structure.from_file(cif_path)
        bulk = ensure_oxidation_states(bulk)
    except Exception as e:
        raise HTTPException(400, f"Could not parse CIF: {e}")

    bulk_an = ChemicalAnalyzer(bulk)
    if not bulk_an.centers:
        raise HTTPException(400, "No polyanion scaffold detected in this structure.")

    bulk_json = serialize_structure(bulk, analyzer=bulk_an)
    analysis_json = serialize_analysis(bulk_an)

    sessions[session_id] = {
        "cif_path": cif_path,
        "bulk_len": len(bulk),
    }

    return _json_response({
        "session_id": session_id,
        "bulk": bulk_json,
        "analysis": analysis_json,
        "filename": file.filename,
    })


class MillerInput(BaseModel):
    mode: str   # "max" | "single" | "list"
    value: object


class GenerateParams(BaseModel):
    session_id: str
    thickness: float = 15.0
    vacuum: float = 15.0
    dipole_tol: float = 0.3
    miller_input: MillerInput
    thickness_range: list[float] | None = None      # [min, max] Å → enumerate all terminations
    allow_nonstoichiometric: bool = False


@app.post("/generate")
async def generate(params: GenerateParams):
    if params.session_id not in sessions:
        raise HTTPException(404, "Session not found. Please upload a CIF first.")

    job_id = str(uuid4())
    jobs[job_id] = {
        "status": "running",
        "progress": [],
        "results": [],
        "raw_slabs": [],
        "total_indices": 0,
        "completed_indices": 0,
        "current_hkl": None,
        "error": None,
    }

    executor.submit(
        run_job,
        job_id,
        params.session_id,
        params.model_dump(),
        sessions,
        jobs,
    )

    return {"job_id": job_id}


@app.get("/status/{job_id}")
def status(job_id: str):
    if job_id not in jobs:
        raise HTTPException(404, "Job not found.")
    job = jobs[job_id]
    # Return everything except raw_slabs (not serialisable here)
    safe = {k: v for k, v in job.items() if k != "raw_slabs"}
    return _json_response(safe)


class OptimizeParams(BaseModel):
    session_id: str
    job_id: str


@app.post("/optimize")
async def optimize(params: OptimizeParams):
    if params.session_id not in sessions:
        raise HTTPException(404, "Session not found.")
    if params.job_id not in jobs:
        raise HTTPException(404, "Generation job not found.")
    if jobs[params.job_id].get('status') != 'done':
        raise HTTPException(400, "Generation job not complete yet.")

    opt_job_id = str(uuid4())
    opt_jobs[opt_job_id] = {
        "status":       "pending",
        "phase":        "init",
        "phase_label":  "Initializing…",
        "current":      0,
        "completed":    0,
        "total":        0,
        "e_bulk":       None,
        "slab_results": [],
        "wulff_data":   None,
        "wulff_error":  None,
        "error":        None,
        "traceback":    None,
    }

    opt_executor.submit(
        run_optimization_job,
        opt_job_id,
        params.job_id,
        params.session_id,
        sessions,
        jobs,
        opt_jobs,
    )

    return {"opt_job_id": opt_job_id}


@app.get("/opt_status/{opt_job_id}")
def opt_status(opt_job_id: str):
    if opt_job_id not in opt_jobs:
        raise HTTPException(404, "Optimization job not found.")
    return _json_response(opt_jobs[opt_job_id])


@app.get("/download/{job_id}/{slab_index}")
def download(job_id: str, slab_index: int):
    if job_id not in jobs:
        raise HTTPException(404, "Job not found.")
    raw = jobs[job_id].get("raw_slabs", [])
    if slab_index >= len(raw):
        raise HTTPException(404, "Slab index out of range.")

    slab = raw[slab_index]
    hkl = "".join(str(x) for x in slab.miller_index)
    tmp = tempfile.NamedTemporaryFile(suffix=".cif", delete=False)
    tmp.close()
    slab.to(filename=tmp.name)
    return FileResponse(tmp.name, media_type="chemical/x-cif",
                        filename=f"slab_{hkl}_{slab_index}.cif")
