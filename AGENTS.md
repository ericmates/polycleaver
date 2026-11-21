# AGENTS.md

Instructions for coding agents (and people) who build surface slabs with
polycleaver or work on its code. The user-facing documentation is in
[README.md](README.md); the design decisions are explained in
[docs/DESIGN_NOTES.md](docs/DESIGN_NOTES.md).

## Setting up

```bash
pip install -e ".[test,mcp]"        # or: conda env create -f environment.yml
polycleaver --version
```

## Building slabs

Prefer the command-line options and the JSON output to writing scripts:

```bash
polycleaver -s BULK.cif --analyse-only --json                    # classification only
polycleaver -s BULK.cif --hkl 1 1 0 --thickness 15 -o OUT --json  # thinnest valid slab per cut
polycleaver -s BULK.cif --hkl 1 1 0 --range 12 20 -o OUT --json   # every termination in 12-20 Å
polycleaver -s BULK.cif --max-index 1 --thickness 15 -o OUT --json
```

With `--json`, standard output holds one JSON document and the readable report
goes to standard error. The same document is written to
`OUT/polycleaver_summary.json`. Read the results from its keys; do not parse the
readable report.

Checks to make and report after every run:

1. `bulk.charge_neutral` is true and `bulk.scaffold_centres`, `bulk.cations` and
   `bulk.molecules` match the mineral. If not, stop and report the problem.
2. For the slab that is used: `stoichiometric` and `dipole_within_tolerance` are
   true; report `thickness_A`, `symmetric` and `dipole_z_eA`.
3. `rank` orders the slabs by `broken_bonds_per_A2` (lower is more stable). It is
   an estimate for unrelaxed slabs, not a surface energy.
4. If no slab is found, follow `hints` in the summary. Do not raise the dipole
   tolerance above 1 e·Å without the agreement of the user, and always state the
   value used.
5. Keep `polycleaver_summary.json` with the CIF files: it records the parameters
   and the versions of polycleaver, pymatgen, ase and numpy.

## Working on the code

- `polycleaver/slabgen.py`: chemical analysis (`ChemicalAnalyzer`), stripping and
  thinning (`SlabProcessor`), combinatorial placement and stability scoring.
- `polycleaver/cli.py`: command-line interface, `run()` and `analyse()`.
- `polycleaver/summary.py`: the JSON records. Adding keys is allowed; removing or
  renaming a key requires increasing `SUMMARY_SCHEMA_VERSION`.
- `polycleaver/mcp_server.py`: MCP server (stdio). Nothing in it may write to
  standard output, which carries the protocol messages.
- `api/` and `frontend/`: Surface Lab web interface.
- `examples/configs/`: reference runs, whose expected results are listed in the
  README and checked in `tests/test_generation.py`.

Run the tests before committing:

```bash
pytest tests/ -m "not slow"     # seconds
pytest tests/                   # full generation runs, under a minute
```

Changes to bonding cut-offs, the size of the combinatorial pool, the charge
check or the stability metric change the slabs that are produced. Such changes
need a test on the example structures and a note in `CHANGELOG.md`. Documentation
uses British spelling.
