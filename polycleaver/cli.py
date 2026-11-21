#!/usr/bin/env python3
"""
Command-line entry point for polycleaver.

    polycleaver                              # uses ./config.json
    polycleaver my_run.json                  # alternate config file
    polycleaver my_run.json -o out           # write CIFs to ./out
    polycleaver -s bulk.cif --hkl 1 1 0      # no config file needed
    polycleaver -s bulk.cif --analyse-only   # bulk chemical analysis only
    polycleaver my_run.json --json           # JSON summary on stdout, report on stderr
    python -m polycleaver ...                # equivalent

All run parameters can live in a JSON config file (see ``load_config`` in
``polycleaver.slabgen`` and the ``examples/configs`` directory); command-line
options override the values read from the file.
"""

import sys
import json
import argparse
import contextlib
import traceback
from datetime import datetime, timezone
from pathlib import Path

from pymatgen.core import Structure
from pymatgen.core.surface import get_symmetrically_distinct_miller_indices

from . import __version__
from .slabgen import (
    ChemicalAnalyzer,
    ensure_oxidation_states,
    load_config,
    generate_slabs_for_hkl,
    score_surface_stability,
    report_slab_detail,
    print_no_surfaces_diagnostic,
    _compute_bulk_cn_all,
    _deduplicate_slabs,
    _nonstoich_charges,
)
from .summary import (SUMMARY_SCHEMA_VERSION, NO_SLAB_HINTS, bulk_record, slab_record,
                      environment_versions)

SUMMARY_FILENAME = "polycleaver_summary.json"


def build_parser():
    parser = argparse.ArgumentParser(
        prog="polycleaver",
        description="polycleaver: non-polar, stoichiometric surface slabs of "
                    "complex minerals from a bulk structure.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Edit config.json to set your parameters, then run:\n"
            "  polycleaver\n\n"
            "Or point to a different config file:\n"
            "  polycleaver my_run.json\n\n"
            "Or give the parameters on the command line:\n"
            "  polycleaver -s forsterite.cif --hkl 1 1 0 --range 12 20 -o out/\n\n"
            "Example configs live in examples/configs/ of the source repository."
        ),
    )
    parser.add_argument(
        "config", nargs="?", default=None, metavar="CONFIG_FILE",
        help="JSON config file (default: ./config.json, unless --structure is given)",
    )
    parser.add_argument(
        "-o", "--outdir", default=".", metavar="DIR",
        help="directory for the CIF files and the JSON summary (default: current directory)",
    )

    ov = parser.add_argument_group("parameters (override the config file)")
    ov.add_argument("-s", "--structure", metavar="FILE",
                    help="bulk structure file (CIF, POSCAR, ... anything pymatgen reads)")
    ov.add_argument("--hkl", nargs=3, type=int, action="append", metavar=("H", "K", "L"),
                    help="Miller index; repeat the option for several indices")
    ov.add_argument("--max-index", type=int, metavar="N",
                    help="all symmetrically distinct Miller indices up to N")
    ov.add_argument("--thickness", type=float, metavar="ANGSTROM",
                    help="minimum slab thickness (target_thickness)")
    ov.add_argument("--range", nargs=2, type=float, metavar=("MIN", "MAX"), dest="thickness_range",
                    help="return every valid termination between MIN and MAX Å")
    ov.add_argument("--vacuum", type=float, metavar="ANGSTROM", help="vacuum on each side")
    ov.add_argument("--dipole", type=float, metavar="EANGSTROM",
                    help="maximum |dipole_z| in eÅ (dipole_threshold)")
    ov.add_argument("--allow-nonstoichiometric", action="store_true",
                    help="keep charged slabs (written as slab_hkl_N_NS.cif)")

    out = parser.add_argument_group("output")
    out.add_argument("--analyse-only", "--analyze-only", action="store_true",
                     help="report the bulk chemical analysis and stop")
    out.add_argument("--json", action="store_true",
                     help="print the JSON summary on stdout and send the report to stderr")
    out.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    return parser


def config_from_args(args):
    """
    Combine the config file (if any) with the command-line overrides.

    Returns ``(raw_config_dict, base_dir)`` where *base_dir* is the directory
    against which a relative structure path is resolved.
    """
    raw, base_dir = {}, Path.cwd()
    config_path = args.config
    if config_path is None and args.structure is None:
        config_path = "config.json"
    if config_path is not None:
        cfg_file = Path(config_path)
        if not cfg_file.exists():
            raise FileNotFoundError(
                f"Config file not found: {cfg_file}\n"
                f"  Copy examples/configs/forsterite_110_range.json to config.json and set "
                f"'filename', or pass the structure directly: polycleaver -s bulk.cif --hkl 1 1 0"
            )
        try:
            raw = json.loads(cfg_file.read_text())
        except json.JSONDecodeError as exc:
            raise ValueError(f"Invalid JSON in {cfg_file}: {exc}") from exc
        base_dir = cfg_file.resolve().parent

    if args.structure is not None:
        raw["filename"] = str(Path(args.structure).resolve())
    if args.hkl:
        if len(args.hkl) == 1:
            raw["mode"], raw["single_hkl"] = "single", list(args.hkl[0])
        else:
            raw["mode"], raw["custom_list"] = "list", [list(h) for h in args.hkl]
    if args.max_index is not None:
        raw["mode"], raw["max_index"] = "max", args.max_index
    if args.thickness is not None:
        # An explicit target thickness switches off a range read from the file.
        raw["target_thickness"] = args.thickness
        raw["thickness_range"] = None
    if args.thickness_range is not None:
        raw["thickness_range"] = list(args.thickness_range)
    if args.vacuum is not None:
        raw["vacuum_size"] = args.vacuum
    if args.dipole is not None:
        raw["dipole_threshold"] = args.dipole
    if args.allow_nonstoichiometric:
        raw["allow_nonstoichiometric"] = True
    return raw, base_dir


def _resolve_structure_path(filename, base_dir):
    """A relative path is tried against *base_dir* first, then the current directory."""
    cif_path = Path(filename)
    if not cif_path.is_absolute() and not cif_path.exists() and base_dir is not None:
        candidate = Path(base_dir) / cif_path
        if candidate.exists():
            cif_path = candidate
    return cif_path.resolve() if cif_path.exists() else cif_path


def _load_bulk(cif_path):
    bulk = Structure.from_file(str(cif_path))
    return ensure_oxidation_states(bulk)


def analyse(structure_path, verbose=True):
    """
    Run the bulk chemical analysis only and return its JSON record.

    The record lists the scaffold centres, scaffold anions, counter-cations,
    molecular species, hydration targets and whether the cell is charge-neutral.
    """
    cif_path = Path(structure_path)
    bulk = _load_bulk(cif_path)
    an = ChemicalAnalyzer(bulk)
    if verbose:
        an.print_detailed()
    rec = bulk_record(bulk, an, _compute_bulk_cn_all(bulk))
    rec["file"] = str(cif_path.resolve())
    return {
        "polycleaver_version": __version__,
        "schema_version": SUMMARY_SCHEMA_VERSION,
        "versions": environment_versions(),
        "bulk": rec,
    }


def run(config="config.json", outdir=".", base_dir=None, write_summary=True,
        config_label=None):
    """
    Run the full pipeline and write one CIF per accepted termination to *outdir*.

    Parameters
    ----------
    config : str, Path or dict
        JSON config file or a dict with the same keys (see ``load_config``).
    outdir : str or Path
        Output directory, created if needed.
    base_dir : str or Path, optional
        Directory against which a relative ``filename`` is resolved.  Defaults
        to the directory of the config file, or the current directory.
    write_summary : bool
        Write ``polycleaver_summary.json`` next to the CIF files.
    config_label : str, optional
        How the configuration is named in the report header.

    Returns
    -------
    dict
        The run summary (see ``polycleaver.summary``): versions, the resolved
        configuration, the bulk analysis, one record per accepted slab with
        its file path, and warnings.  Raises on configuration errors.
    """
    cfg = load_config(config)
    if base_dir is None and not isinstance(config, dict):
        base_dir = Path(config).resolve().parent
    filename                = cfg["filename"]
    target_thickness        = cfg["target_thickness"]
    vacuum_size             = cfg["vacuum_size"]
    dipole_threshold        = cfg["dipole_threshold"]
    mode                    = cfg["mode"]
    single_hkl              = cfg["single_hkl"]
    max_index               = cfg["max_index"]
    custom_list             = cfg["custom_list"]
    thickness_range         = cfg["thickness_range"]
    auto_select_stable      = cfg["auto_select_stable"]
    allow_nonstoichiometric = cfg["allow_nonstoichiometric"]

    cif_path = _resolve_structure_path(filename, base_dir)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    warnings = []

    # ------------------------------------------------------------------
    #  Load bulk structure
    # ------------------------------------------------------------------
    print("=" * 68)
    print(f"  polycleaver  —  non-polar surface slab generator  (v{__version__})")
    if config_label is None:
        config_label = "(parameters)" if isinstance(config, dict) else str(config)
    print(f"  Config : {config_label}")
    print(f"  File   : {cif_path}")
    print("=" * 68)

    bulk = _load_bulk(cif_path)

    bulk_charge = sum(site.specie.oxi_state for site in bulk)
    if abs(bulk_charge) > 0.01:
        msg = (f"The bulk cell carries a net charge of {bulk_charge:+.2f} e with the "
               f"assigned oxidation states; no stoichiometric slab can be charge-neutral. "
               f"Check the structure (missing H / OH?, partial occupancies) or add the "
               f"right entries to OXI_MAP.")
        warnings.append(msg)
        print(f"\n  WARNING: {msg}")

    # ------------------------------------------------------------------
    #  Analyse bulk chemistry
    # ------------------------------------------------------------------
    print("\n  ┌─ BULK ANALYSIS ")
    print("  │")
    bulk_an = ChemicalAnalyzer(bulk)
    bulk_an.print_detailed()
    print("  │")

    if not bulk_an.centers:
        print("  │  ERROR: No polyanion scaffold detected.")
        print("  │  The structure may be purely ionic or elemental.")
        print("  └" + "─" * 60 + "┘")
        raise ValueError("No polyanion scaffold detected in the bulk structure.")

    if bulk_an.has_shared_oxygens():
        framework_type = (
            "Framework/sheet (M–O–M bridges present) — "
            "under-coordinated surface centres will be kept intact"
        )
    else:
        framework_type = (
            "Isolated polyhedra — incomplete surface units "
            "will be removed during stripping"
        )
    print(f"  │  Framework          : {framework_type}")
    print("  └" + "─" * 60 + "┘")

    # Bulk element coordination numbers (used for stability scoring)
    bulk_cn_all = _compute_bulk_cn_all(bulk)
    if bulk_cn_all:
        cn_str = "  ".join(f"{el}: {cn:.1f}" for el, cn in sorted(bulk_cn_all.items()))
        print(f"\n  Bulk element CN      : {cn_str}")

    # ------------------------------------------------------------------
    #  Resolve target Miller indices
    # ------------------------------------------------------------------
    if mode == "single":
        target_indices = [single_hkl]
        print(f"\n  Mode  : Single index {single_hkl}")
    elif mode == "max":
        print(f"\n  Mode  : All symmetrically distinct indices up to {max_index}")
        target_indices = get_symmetrically_distinct_miller_indices(bulk, max_index)
        print(f"  Found : {len(target_indices)} distinct indices: {target_indices}")
    else:  # 'list'
        target_indices = custom_list
        print(f"\n  Mode  : Custom list: {target_indices}")

    # ------------------------------------------------------------------
    #  Generation loop
    # ------------------------------------------------------------------
    print("\n  Generation parameters:")
    if thickness_range is not None:
        print(f"    Thickness range  : {thickness_range[0]:.1f} – {thickness_range[1]:.1f} Å")
    else:
        print(f"    Target thickness : {target_thickness:.1f} Å")
    print(f"    Vacuum           : {vacuum_size:.1f} Å")
    print(f"    Dipole tolerance : {dipole_threshold}")
    print(f"    Stability sort   : {'Yes' if auto_select_stable else 'No'}")
    print(f"    Non-stoichiom.   : {'Yes — _NS files kept' if allow_nonstoichiometric else 'No — discarded'}")
    print(f"    Output directory : {outdir}")
    print()

    all_results = []
    for hkl in target_indices:
        slabs = generate_slabs_for_hkl(
            bulk, bulk_an, hkl,
            target_thickness, vacuum_size, dipole_threshold,
            thickness_range=thickness_range,
            allow_nonstoichiometric=allow_nonstoichiometric,
        )
        all_results.extend(slabs)

    # ------------------------------------------------------------------
    #  Stability sort + deduplication
    # ------------------------------------------------------------------
    if all_results:
        scores = [score_surface_stability(s, bulk_cn_all) for s in all_results]

        if auto_select_stable:
            paired = list(zip(scores, all_results))
            paired.sort(key=lambda x: (
                x[0][2],                        # broken_per_area ASC (primary)
                -int(x[1].is_symmetric()),      # symmetric first (tiebreaker)
                -x[0][0],                       # coord_saturation DESC (tiebreaker)
            ))
            scores = [sc for sc, _ in paired]
            all_results = [sl for _, sl in paired]

        all_results = _deduplicate_slabs(all_results)
        scores = [score_surface_stability(s, bulk_cn_all) for s in all_results]
    else:
        scores = []

    # ------------------------------------------------------------------
    #  Final report and CIF output
    # ------------------------------------------------------------------
    print("\n" + "=" * 68)
    print("  FINAL REPORT")
    print("=" * 68)

    if not all_results:
        print_no_surfaces_diagnostic()
        warnings.append("No slab passed the dipole and charge filters.")
    else:
        print(f"\n  Accepted {len(all_results)} unique slab(s) "
              f"across {len(target_indices)} Miller index/ices.\n")

    slab_records = []
    for i, slab in enumerate(all_results):
        stability = scores[i] if scores else None
        rank = (i + 1) if (auto_select_stable and scores) else None
        ns_charge = _nonstoich_charges.get(id(slab))
        fname = report_slab_detail(
            slab, i, len(bulk),
            target_thickness, dipole_threshold,
            stability=stability, stability_rank=rank,
            thickness_range=thickness_range,
            net_charge=ns_charge,
        )
        out_path = outdir / fname
        slab.to(filename=str(out_path))
        print(f"  ✓ Wrote  {out_path}")
        rec = slab_record(slab, fname, len(bulk), stability=stability, rank=rank,
                          net_charge=ns_charge, dipole_tol=dipole_threshold)
        rec["path"] = str(out_path.resolve())
        slab_records.append(rec)

    if slab_records:
        print(f"\n  ─── {len(slab_records)} CIF file(s) written ───")

    resolved_cfg = dict(cfg)
    resolved_cfg["filename"] = str(cif_path)
    resolved_cfg["single_hkl"] = list(single_hkl)
    resolved_cfg["custom_list"] = [list(h) for h in custom_list]
    bulk_rec = bulk_record(bulk, bulk_an, bulk_cn_all)
    bulk_rec["file"] = str(cif_path)
    summary = {
        "polycleaver_version": __version__,
        "schema_version": SUMMARY_SCHEMA_VERSION,
        "created": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "versions": environment_versions(),
        "config": resolved_cfg,
        "bulk": bulk_rec,
        "miller_indices": [[int(x) for x in hkl] for hkl in target_indices],
        "n_slabs": len(slab_records),
        "ranking": ("broken_bonds_per_A2 ascending, then symmetric first, then "
                    "coordination_saturation descending") if auto_select_stable else None,
        "slabs": slab_records,
        "warnings": warnings,
        "hints": NO_SLAB_HINTS if not slab_records else [],
        "output_dir": str(outdir.resolve()),
    }
    if write_summary:
        summary_path = outdir / SUMMARY_FILENAME
        summary["summary_file"] = str(summary_path.resolve())
        summary_path.write_text(dumps(summary))
        print(f"  ✓ Wrote  {summary_path}")
    return summary


def _json_default(obj):
    """Fallback for numpy scalars and other simple objects in the summary."""
    for attr in ("item", "tolist"):
        if hasattr(obj, attr):
            return getattr(obj, attr)()
    return str(obj)


def dumps(obj):
    """Serialise a polycleaver record to JSON text."""
    return json.dumps(obj, indent=2, ensure_ascii=False, default=_json_default)


def main(argv=None):
    args = build_parser().parse_args(argv)
    # With --json, stdout carries only the JSON document; the report goes to stderr.
    report_stream = sys.stderr if args.json else sys.stdout
    try:
        raw, base_dir = config_from_args(args)
        with contextlib.redirect_stdout(report_stream):
            if args.analyse_only:
                if not raw.get("filename"):
                    raise ValueError("--analyse-only needs a structure "
                                     "(-s FILE, or 'filename' in the config file).")
                path = _resolve_structure_path(raw["filename"], base_dir)
                result = analyse(path, verbose=True)
            else:
                label = args.config or ("config.json" if args.structure is None else None)
                if label and args.structure is not None:
                    label += " + command-line options"
                result = run(raw, args.outdir, base_dir=base_dir, config_label=label)
        if args.json:
            print(dumps(result))
    except (FileNotFoundError, ValueError) as exc:
        print(f"\n  ERROR: {exc}", file=report_stream)
        if args.json:
            print(dumps({"error": str(exc), "polycleaver_version": __version__}))
        sys.exit(1)
    except Exception as exc:
        print(f"\n  CRITICAL ERROR: {exc}", file=report_stream)
        traceback.print_exc()
        if args.json:
            print(dumps({"error": f"{type(exc).__name__}: {exc}",
                         "polycleaver_version": __version__}))
        sys.exit(1)


if __name__ == "__main__":
    main()
