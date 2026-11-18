"""License-free integrity checks; no model is optimized or solver imported."""
import ast
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main():
    checked = 0
    for path in ROOT.rglob('*.py'):
        if 'legacy' in path.parts or 'source_as_received' in path.parts:
            continue
        ast.parse(path.read_text(encoding='utf-8-sig'), filename=str(path), feature_version=8)
        checked += 1
    spec = importlib.util.spec_from_file_location('verify_areas', ROOT/'06_validation/verify_areas.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    results, triangles = module.check_file(ROOT/'data/coordinates/appendix_v15.csv')
    expected = json.loads((ROOT/'data/coordinates/appendix_expected.json').read_text())
    for result in results:
        target = expected[str(result['n'])]
        assert result['points_in_unit_square']
        assert result['minimum_area_exact_decimal'] == target['exact_minimum_area']
        assert result['one_minimizing_triple'] == target['minimizing_triple']
    assert len(results) == 7 and len(triangles) == 700
    for item in json.loads((ROOT/'data/recorded_runs/archive_integrity.json').read_text()):
        digest = hashlib.sha256((ROOT/item['path']).read_bytes()).hexdigest()
        assert digest == item['published_sha256'], item['path']
    groups = json.loads((ROOT/'03_enhancement_comparison/solver_settings.json').read_text())
    final = json.loads((ROOT/'04_final_results/solver_settings.json').read_text())
    assert groups['MIPGap'] == 1e-4 and groups['MIPGapAbs'] == 1e-10
    assert final['MIPGap'] == 1e-6 and final['MIPGapAbs'] == 0
    assert final['Threads'] == 32 and final['Symmetry'] == 2
    print('PASS: {} Python files; 700 exact triangle areas; archive hashes; distinct parameter profiles.'.format(checked))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
