#!/usr/bin/env python3
"""
Command-line entry point for PolyCleaver 2.

    polycleaver2                    # uses ./config.json
    polycleaver2 my_run.json        # alternate config file
    polycleaver2 my_run.json -o out # write CIFs to ./out
    python -m polycleaver2 …        # equivalent

All run parameters live in the JSON config file (see ``load_config`` in
``polycleaver2.slabgen`` and the ``examples/configs`` directory).
"""

import sys
import argparse
import traceback
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


def build_parser():
    parser = argparse.ArgumentParser(
        prog="polycleaver2",
        description="PolyCleaver 2 — non-polar, stoichiometric surface slab "
                    "generator for complex minerals.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Edit config.json to set your parameters, then run:\n"
            "  polycleaver2\n\n"
            "Or point to a different config file:\n"
            "  polycleaver2 my_run.json\n\n"
            "Example configs live in examples/configs/."
        ),
    )
    parser.add_argument(
        "config",
        nargs="?",
        default="config.json",
        metavar="CONFIG_FILE",
        help="path to JSON config file (default: config.json)",
    )
    parser.add_argument(
        "-o", "--outdir",
        default=".",
        metavar="DIR",
        help="directory for the output CIF files (default: current directory)",
    )
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {__version__}",
    )
    return parser


def run(config_path, outdir="."):
    """
    Run the full pipeline described by *config_path* and write CIFs to *outdir*.

    Returns the list of written CIF paths.  Raises on configuration errors.
    """
    cfg = load_config(config_path)
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

    # A relative 'filename' is resolved against the config file's directory
    # first, then against the current working directory.
    cif_path = Path(filename)
    if not cif_path.is_absolute() and not cif_path.exists():
        candidate = Path(config_path).resolve().parent / cif_path
        if candidate.exists():
            cif_path = candidate

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    #  Load bulk structure
    # ------------------------------------------------------------------
    print("=" * 68)
    print(f"  POLYCLEAVER 2  —  non-polar surface slab generator  (v{__version__})")
    print(f"  Config : {config_path}")
    print(f"  File   : {cif_path}")
    print("=" * 68)

    bulk = Structure.from_file(str(cif_path))
    bulk = ensure_oxidation_states(bulk)

    bulk_charge = sum(site.specie.oxi_state for site in bulk)
    if abs(bulk_charge) > 0.01:
        print(f"\n  WARNING: the bulk cell carries a net charge of {bulk_charge:+.2f} e with the "
              f"assigned oxidation states.\n"
              f"  No stoichiometric slab can be charge-neutral. Check the structure "
              f"(missing H / OH?, partial occupancies) or add the right entries to OXI_MAP.")

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
    else:
        print(f"\n  Accepted {len(all_results)} unique slab(s) "
              f"across {len(target_indices)} Miller index/ices.\n")

    cif_paths = []
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
        cif_paths.append(out_path)

    if cif_paths:
        print(f"\n  ─── {len(cif_paths)} CIF file(s) written ───")

    return cif_paths


def main(argv=None):
    args = build_parser().parse_args(argv)
    try:
        run(args.config, args.outdir)
    except (FileNotFoundError, ValueError) as exc:
        print(f"\n  ERROR: {exc}")
        sys.exit(1)
    except Exception as exc:
        print(f"\n  CRITICAL ERROR: {exc}")
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
