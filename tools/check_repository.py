"""License-free integrity checks; no model is optimized or solver imported."""
import ast
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]


def main():
    checked = 0
    for path in ROOT.rglob("*.py"):
        if "results" in path.relative_to(ROOT).parts:
            continue
        ast.parse(path.read_text(encoding="utf-8-sig"), filename=str(path), feature_version=8)
        checked += 1
    spec = importlib.util.spec_from_file_location("verify_areas", ROOT / "06_validation/verify_areas.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    # Known analytic configurations exercise CSV, nested JSON, repeated n,
    # missing incumbents and out-of-square coordinates without bundled data.
    with tempfile.TemporaryDirectory() as directory:
        folder = Path(directory)
        csv_path = folder / "coordinates.csv"
        csv_path.write_text("n,i,x,y\n3,1,0,0\n3,2,1,0\n3,3,0,1\n", encoding="utf-8")
        result, triangles = module.check_file(csv_path)
        assert result[0]["minimum_area_exact_fraction"] == "1/2"
        assert result[0]["feasible_lower_bound_7dp"] == "0.5000000"
        assert len(triangles) == 1
        json_path = folder / "runs.json"
        json_path.write_text(json.dumps({"runs": [
            {"n": 3, "x": [0, 1, 0], "y": [0, 0, 1]},
            {"n": 3, "x": [0, 1, 0], "y": [0, 0, "1.000000001"]},
            {"n": 3, "x": None, "y": None},
        ]}), encoding="utf-8")
        result, triangles, skipped = module.check_inputs([json_path])
        assert len(result) == 2 and len(triangles) == 2
        assert result[0]["points_in_unit_square"]
        assert not result[1]["points_in_unit_square"]
        assert result[1]["feasible_lower_bound_7dp"] is None
        assert result[1]["minimum_area_exact_fraction"] == "1000000001/2000000000"
        # A folder prefers individual result files over aggregate copies.
        run = folder / "run"; run.mkdir()
        (run / "result.json").write_text(json.dumps(
            {"n": 3, "x": [0, 1, 0], "y": [0, 0, 1]}), encoding="utf-8")
        result, _, _ = module.check_inputs([folder])
        assert len(result) == 1
    groups = json.loads((ROOT / "03_enhancement_comparison/solver_settings.json").read_text())
    final = json.loads((ROOT / "04_final_results/solver_settings.json").read_text())
    assert groups["MIPGap"] == 1e-4 and groups["MIPGapAbs"] == 1e-10
    assert final["MIPGap"] == 1e-6 and final["MIPGapAbs"] == 0
    assert final["Threads"] == 32 and final["Symmetry"] == 2
    # --list exits before solver imports, so verify full table coverage license-free.
    def listed(relative):
        return subprocess.check_output([sys.executable, str(ROOT / relative), "--list"],
                                       text=True, encoding="utf-8")
    formulations = [json.loads(line) for line in listed("01_formulations/formulations.py").splitlines()
                    if line.startswith("{")]
    expected_formulations = {(suite, n, approach)
                             for suite in ("baseline", "shared")
                             for n in (range(3, 8) if suite == "baseline" else (7, 8))
                             for approach in ((1, 2, 3) if suite == "baseline" else (1, 3))}
    assert len(formulations) == 19
    assert {(r["suite"], r["n"], r["approach"]) for r in formulations} == expected_formulations
    solver_cases = [line.split() for line in listed("02_solver_comparison/solver_comparison.py").splitlines()
                    if line.strip()]
    assert len(solver_cases) == 20
    assert {(solver, int(n)) for solver, n in solver_cases} == {
        (solver, n) for solver in ("gurobi", "baron", "scip", "couenne") for n in range(3, 8)}
    enhancement_cases = json.loads(listed("03_enhancement_comparison/enhancement_comparison.py"))
    all_groups = {"none", "1", "2", "3", "4", "12", "13", "14", "23", "24", "34",
                  "123", "124", "134", "234", "1234"}
    assert len(enhancement_cases) == 48
    assert {(r["n"], r["groups"]) for r in enhancement_cases} == {
        (n, group) for n in (7, 8, 9) for group in all_groups}
    final_cases = json.loads(listed("04_final_results/final_results.py"))
    assert len(final_cases) == 10 and {r["n"] for r in final_cases} == set(range(3, 13))
    print("PASS: {} Python files; exact CSV/JSON configuration checks; distinct parameter profiles; "
          "19 formulation, 20 solver, 48 enhancement and 10 final cases.".format(checked))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
