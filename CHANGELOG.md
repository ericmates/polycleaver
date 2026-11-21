# Changelog

All notable changes to polycleaver are recorded in this file. Version numbers
follow [semantic versioning](https://semver.org).

## [1.0.0] – 2026-11-21

First stable release of the rewritten polycleaver. The code was developed during
2026 in the separate repository `ericmates/polycleaver2` under the name
PolyCleaver 2, and its history is merged into this repository. The Python
interface of 0.x is not kept, so scripts written for 0.x need to be adapted or
pinned to `polycleaver==0.0.8`.

### Added

- Chemical classification of the bulk into scaffold centres, scaffold oxygens,
  counter-cations and discrete molecules (H₂O, OH⁻, NH₃, NH₄⁺, …), with the modal
  hydration number of each cation.
- Treatment of isolated-polyhedra, chain, sheet and framework scaffolds, and of
  structural hydroxyl groups attached to their polyhedra.
- Combinatorial placement of cations and molecules that minimises the dipole
  along the surface normal, with water placed next to the cations.
- Range mode (`thickness_range`), which returns every valid termination between
  two thicknesses, and ranking of the terminations by broken bonds per Å².
- Ordering of disordered CIF files with the Ewald-ranked transformation of
  pymatgen, and a warning when the bulk cell is not charge-neutral.
- Optional non-stoichiometric slabs for reactivity studies (`allow_nonstoichiometric`).
- Command-line interface `polycleaver` driven by a JSON configuration file, with
  the options `-s/--structure`, `--hkl`, `--max-index`, `--thickness`, `--range`,
  `--vacuum`, `--dipole`, `--allow-nonstoichiometric`, `--analyse-only` and `--json`.
- Machine-readable run summary `polycleaver_summary.json`, which records the
  library versions, the parameters, the bulk classification and the properties
  of every slab.
- Python functions `polycleaver.run` and `polycleaver.analyse`, which return the
  same records.
- Model Context Protocol server `polycleaver-mcp` (extra `mcp`, compatible with
  the 1.x and 2.x series of the `mcp` SDK) and an `AGENTS.md` file for LLM agents.
- Surface Lab web interface (FastAPI backend and React/Three.js frontend), with
  optional MACE-MP relaxations and Wulff shapes.
- Conda environment file, a conda recipe for local builds, a test suite, and
  GitHub Actions workflows for the tests and for publishing to PyPI.

### Changed

- The package is configured with `pyproject.toml` and requires Python ≥ 3.10.
- The documentation is written in Markdown and built with MkDocs.

### Removed

- The 0.x Python interface (`polycleaver.core.mnx`, `mx`, `tools`, `miscellaneous`).
- Support for disulphides such as pyrite (FeS₂), whose S₂²⁻ dimers are not
  recognised as a scaffold by the new classification.

### Changes with respect to the PolyCleaver 2 development snapshot (2.0.0a1, September 2026)

- The package, the module and the command are renamed from `polycleaver2` to
  `polycleaver`.
- `run()` returns the run summary (a dict) instead of the list of CIF paths, and
  `load_config()` also accepts a dict.
- New command-line options, the JSON summary, the MCP server and `AGENTS.md`.

## [0.0.8] – 2024-03-11

Last release of the original code, described in E. Mates-Torres and A. Rimola,
*J. Appl. Cryst.* **57**, 503–508 (2024). Releases 0.0.1 to 0.0.7 were published
on PyPI between July and August 2023. The code is preserved at the tag `v0.0.8`.

[1.0.0]: https://github.com/ericmates/polycleaver/releases/tag/v1.0.0
[0.0.8]: https://github.com/ericmates/polycleaver/tree/v0.0.8
