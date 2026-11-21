"""
Tests of the command-line interface, the JSON summary and the MCP server.

The analysis-only and argument tests are fast; the full runs are marked slow.
"""

import asyncio
import json
import sys
from pathlib import Path

import pytest

from polycleaver import __version__
from polycleaver.cli import build_parser, config_from_args, main

ROOT = Path(__file__).resolve().parent.parent
STRUCTURES = ROOT / "examples" / "structures"
CONFIGS = ROOT / "examples" / "configs"


def test_overrides_without_config_file(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)                    # no config.json here
    args = build_parser().parse_args(
        ["-s", str(STRUCTURES / "forsterite.cif"), "--hkl", "1", "1", "0",
         "--range", "12", "20", "--dipole", "0.5"])
    raw, _ = config_from_args(args)
    assert raw["filename"] == str(STRUCTURES / "forsterite.cif")
    assert raw["mode"] == "single" and raw["single_hkl"] == [1, 1, 0]
    assert raw["thickness_range"] == [12.0, 20.0]
    assert raw["dipole_threshold"] == 0.5


def test_overrides_on_top_of_config_file():
    args = build_parser().parse_args(
        [str(CONFIGS / "forsterite_110_range.json"), "--hkl", "0", "1", "0",
         "--hkl", "1", "0", "0", "--thickness", "10"])
    raw, base_dir = config_from_args(args)
    assert base_dir == CONFIGS
    assert raw["mode"] == "list" and raw["custom_list"] == [[0, 1, 0], [1, 0, 0]]
    assert raw["target_thickness"] == 10.0
    assert raw["thickness_range"] is None          # --thickness switches off the file's range


def test_missing_config_file_message(tmp_path, monkeypatch, capsys):
    monkeypatch.chdir(tmp_path)
    with pytest.raises(SystemExit) as exc:
        main(["--json"])
    assert exc.value.code == 1
    out = json.loads(capsys.readouterr().out)
    assert "Config file not found" in out["error"]


def test_analyse_only_json(capsys):
    main(["-s", str(STRUCTURES / "struvite.cif"), "--analyse-only", "--json"])
    captured = capsys.readouterr()
    rec = json.loads(captured.out)                 # stdout holds only the JSON document
    assert rec["polycleaver_version"] == __version__
    bulk = rec["bulk"]
    assert bulk["charge_neutral"]
    assert bulk["scaffold_centres"] == {"P": 2}
    assert bulk["cations"] == {"Mg": 2}
    assert bulk["molecules"]["H2O"] == 12
    assert bulk["cation_hydration_target"] == {"Mg": 6}
    assert bulk["scaffold_connectivity"] == "isolated_polyhedra"
    assert "Bulk composition" in captured.err      # human-readable report on stderr


@pytest.mark.slow
def test_cli_run_writes_cifs_and_summary(tmp_path, capsys):
    outdir = tmp_path / "out"
    main([str(CONFIGS / "forsterite_010.json"), "-o", str(outdir), "--json"])
    summary = json.loads(capsys.readouterr().out)
    assert summary["n_slabs"] >= 1
    assert summary["miller_indices"] == [[0, 1, 0]]
    on_disk = json.loads((outdir / "polycleaver_summary.json").read_text())
    assert on_disk["n_slabs"] == summary["n_slabs"]
    for rank, rec in enumerate(summary["slabs"], start=1):
        assert rec["rank"] == rank
        assert Path(rec["path"]).exists()
        assert abs(rec["dipole_z_eA"]) <= 0.3 and rec["dipole_within_tolerance"]
        assert rec["stoichiometric"] and rec["reduced_formula"] == "Mg2SiO4"
        assert rec["thickness_A"] >= 10.0
    bpa = [rec["broken_bonds_per_A2"] for rec in summary["slabs"]]
    assert bpa == sorted(bpa)


def test_mcp_server_tools(tmp_path):
    pytest.importorskip("mcp")
    from mcp import ClientSession, StdioServerParameters
    from mcp.client.stdio import stdio_client

    def attr(obj, *names):
        return next(getattr(obj, n) for n in names if hasattr(obj, n))

    def payload(res):
        data = attr(res, "structured_content", "structuredContent")
        if data is None:
            data = json.loads(res.content[0].text)
        return data.get("result", data) if isinstance(data, dict) else data

    async def session():
        params = StdioServerParameters(command=sys.executable,
                                       args=["-m", "polycleaver.mcp_server"])
        async with stdio_client(params) as (r, w):
            async with ClientSession(r, w) as s:
                await s.initialize()
                names = {t.name for t in (await s.list_tools()).tools}
                assert names == {"analyse_structure", "list_miller_indices", "generate_slabs"}
                res = await s.call_tool("analyse_structure",
                                        {"structure_path": str(STRUCTURES / "forsterite.cif")})
                assert not attr(res, "is_error", "isError")
                assert payload(res)["bulk"]["scaffold_centres"] == {"Si": 4}
                res = await s.call_tool("analyse_structure",
                                        {"structure_path": str(tmp_path / "missing.cif")})
                assert attr(res, "is_error", "isError")
                assert "not found" in res.content[0].text

    asyncio.run(session())
