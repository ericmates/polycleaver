"""
Model Context Protocol (MCP) server for polycleaver.

The server exposes the bulk analysis and the slab generation as tools that an
LLM agent (any MCP client) can call on the local machine.  It communicates over
stdio, so the client starts it as a subprocess:

    polycleaver-mcp                 # installed with:  pip install "polycleaver[mcp]"
    python -m polycleaver.mcp_server

The human-readable report that the command-line interface prints is captured
and written to ``polycleaver.log`` in the output directory, because stdout is
reserved for the protocol messages.
"""

import contextlib
import io
import sys
import threading
from functools import partial
from pathlib import Path

from . import __version__

# The server class is MCPServer in the mcp SDK 2.x and FastMCP in 1.x; both
# take the same arguments for what is used here.
_IMPORT_ERROR = None
try:
    import anyio
    from mcp.server.mcpserver import MCPServer as _Server        # mcp >= 2
    from mcp.server.mcpserver.exceptions import ToolError
except ImportError:
    try:
        import anyio
        from mcp.server.fastmcp import FastMCP as _Server        # mcp 1.x
        from mcp.server.fastmcp.exceptions import ToolError
    except ImportError as exc:                 # pragma: no cover - optional extra
        _Server = None
        _IMPORT_ERROR = exc


INSTRUCTIONS = f"""\
polycleaver {__version__} builds non-polar, charge-neutral, stoichiometric surface
slabs of minerals with polyanion scaffolds (silicates, phosphates, sulfates,
hydrated and hydroxylated minerals) from a bulk structure file on this machine.

Recommended use:
1. Call analyse_structure first. Check that bulk.charge_neutral is true and that
   the scaffold centres, cations and molecules match the mineral; if they do not,
   report the problem instead of generating slabs.
2. Call generate_slabs with explicit Miller indices and a thickness (or a
   thickness_range to obtain every termination in a window). Generation may take
   from seconds to several minutes.
3. Read the returned summary: every slab record carries its CIF path, thickness,
   symmetry, residual z-dipole, net charge and the broken-bond stability metric.
   Rank 1 is the termination with the fewest broken bonds per area, which is an
   unrelaxed proxy for the surface energy, not a computed surface energy.
If no slab is found, follow the hints in the summary and state which parameter
was changed; do not raise dipole_threshold above 1 eÅ without telling the user.
The slabs are unrelaxed starting structures for DFT or machine-learned potentials.
"""


# redirect_stdout swaps sys.stdout for the whole process, so tool calls that
# run in worker threads are serialised.
_RUN_LOCK = threading.Lock()


def _capture(func, log_path=None):
    """Run *func* with stdout captured; optionally save the captured report."""
    buf = io.StringIO()
    try:
        with _RUN_LOCK, contextlib.redirect_stdout(buf):
            result = func()
    finally:
        if log_path is not None:
            Path(log_path).parent.mkdir(parents=True, exist_ok=True)
            Path(log_path).write_text(buf.getvalue())
    return result


async def _in_thread(func):
    """
    Run a blocking polycleaver call in a worker thread.

    Failures are re-raised as ToolError, so that the agent receives the reason
    (for example a missing file or a structure without a polyanion scaffold)
    rather than a generic error message.
    """
    try:
        return await anyio.to_thread.run_sync(func)
    except ToolError:
        raise
    except (FileNotFoundError, ValueError) as exc:
        raise ToolError(str(exc)) from exc
    except Exception as exc:
        raise ToolError(f"{type(exc).__name__}: {exc}") from exc


def _analyse(structure_path):
    from .cli import analyse
    path = Path(structure_path).expanduser().resolve()
    if not path.exists():
        raise FileNotFoundError(f"Structure file not found: {path}")
    return _capture(partial(analyse, path, verbose=False))


def _miller(structure_path, max_index):
    from pymatgen.core.surface import get_symmetrically_distinct_miller_indices
    from .cli import _load_bulk
    path = Path(structure_path).expanduser().resolve()
    if not path.exists():
        raise FileNotFoundError(f"Structure file not found: {path}")
    bulk = _capture(partial(_load_bulk, path))
    return [[int(x) for x in hkl]
            for hkl in get_symmetrically_distinct_miller_indices(bulk, int(max_index))]


def _generate(structure_path, output_dir, miller_indices, max_index, target_thickness,
              thickness_range, vacuum, dipole_threshold, allow_nonstoichiometric):
    from .cli import run
    path = Path(structure_path).expanduser().resolve()
    if not path.exists():
        raise FileNotFoundError(f"Structure file not found: {path}")
    outdir = Path(output_dir).expanduser().resolve()
    cfg = {
        "filename": str(path),
        "target_thickness": target_thickness,
        "thickness_range": list(thickness_range) if thickness_range else None,
        "vacuum_size": vacuum,
        "dipole_threshold": dipole_threshold,
        "allow_nonstoichiometric": bool(allow_nonstoichiometric),
    }
    if max_index is not None:
        cfg.update(mode="max", max_index=int(max_index))
    else:
        indices = miller_indices or [[0, 0, 1]]
        if len(indices) == 1:
            cfg.update(mode="single", single_hkl=list(indices[0]))
        else:
            cfg.update(mode="list", custom_list=[list(h) for h in indices])
    log_path = outdir / "polycleaver.log"
    summary = _capture(partial(run, cfg, outdir, config_label="(MCP tool call)"), log_path)
    summary["log_file"] = str(log_path)
    return summary


def build_server():
    """Create the MCP server with the polycleaver tools registered."""
    if _Server is None:
        raise ImportError(
            "The MCP server needs the 'mcp' package: pip install \"polycleaver[mcp]\""
        ) from _IMPORT_ERROR

    server = _Server("polycleaver", instructions=INSTRUCTIONS)

    @server.tool()
    async def analyse_structure(structure_path: str) -> dict:
        """Classify the atoms of a bulk structure (CIF, POSCAR, ...) without generating slabs.

        Returns the scaffold centres (Si, P, S, Al, ...), scaffold anions,
        counter-cations, molecular species (H2O, OH-, NH4+, ...) with their
        counts, the modal hydration number of each cation, the scaffold
        connectivity (isolated polyhedra or framework/sheet), and whether the
        cell is charge-neutral with the assigned oxidation states.
        """
        return await _in_thread(partial(_analyse, structure_path))

    @server.tool()
    async def list_miller_indices(structure_path: str, max_index: int = 1) -> list[list[int]]:
        """List the symmetrically distinct Miller indices of a bulk structure up to max_index."""
        return await _in_thread(partial(_miller, structure_path, max_index))

    @server.tool()
    async def generate_slabs(
        structure_path: str,
        output_dir: str,
        miller_indices: list[list[int]] | None = None,
        max_index: int | None = None,
        target_thickness: float = 15.0,
        thickness_range: list[float] | None = None,
        vacuum: float = 15.0,
        dipole_threshold: float = 0.3,
        allow_nonstoichiometric: bool = False,
    ) -> dict:
        """Generate non-polar, stoichiometric slabs and write one CIF per termination.

        Args:
            structure_path: bulk structure file on this machine.
            output_dir: directory for the CIF files, polycleaver_summary.json
                and polycleaver.log (created if needed).
            miller_indices: list of [h, k, l]; ignored when max_index is given.
                Defaults to [[0, 0, 1]].
            max_index: generate all symmetrically distinct indices up to this value.
            target_thickness: minimum slab thickness in Å (thinnest valid state is kept).
            thickness_range: [min, max] in Å; return every valid termination in
                the window (overrides target_thickness).
            vacuum: vacuum on each side of the slab in Å.
            dipole_threshold: maximum |z-dipole| in e·Å for an accepted slab.
            allow_nonstoichiometric: also keep charged slabs (files tagged _NS).

        Returns the run summary: resolved parameters, bulk analysis, and one
        record per accepted slab (rank, CIF path, formula, thickness, area,
        symmetry, z-dipole, net charge, broken bonds per Å², molecules, surface
        cation hydration), plus warnings and hints when nothing was found.
        """
        return await _in_thread(partial(
            _generate, structure_path, output_dir, miller_indices, max_index,
            target_thickness, thickness_range, vacuum, dipole_threshold,
            allow_nonstoichiometric))

    return server


def main():
    try:
        server = build_server()
    except ImportError as exc:
        print(f"polycleaver-mcp: {exc}", file=sys.stderr)
        sys.exit(1)
    server.run()                               # stdio transport


if __name__ == "__main__":
    main()
