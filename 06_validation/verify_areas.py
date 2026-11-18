"""Exact validation of decimal coordinates. Standard Python only; no solver.

Input CSV columns: n,i,x,y. Decimal strings are treated as exact rational
numbers. The script checks membership in the unit square and every triangle.
It validates configurations, not optimality or solver dual bounds.
"""
import argparse
import csv
from decimal import Decimal, localcontext
from fractions import Fraction
from itertools import combinations
import json
from pathlib import Path


def decimal_string(q):
    with localcontext() as ctx:
        ctx.prec = len(str(abs(q.numerator))) + len(str(q.denominator)) + 20
        return format(Decimal(q.numerator) / Decimal(q.denominator), "f")


def check_file(path):
    groups = {}
    with Path(path).open(encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f)
        if not {"n", "i", "x", "y"}.issubset(reader.fieldnames or []):
            raise ValueError("CSV must contain n,i,x,y columns")
        for row in reader:
            n, i = int(row["n"]), int(row["i"])
            if n < 3 or not 1 <= i <= n:
                raise ValueError(f"Invalid index n={n}, i={i}")
            group = groups.setdefault(n, {})
            if i in group:
                raise ValueError(f"Duplicate point index n={n}, i={i}")
            # Decimal rejects fractions such as 1/3: inputs must be finite decimals.
            xy = tuple(Decimal(row[c]) for c in ("x", "y"))
            if any(not v.is_finite() for v in xy):
                raise ValueError("Coordinates must be finite decimal numbers")
            group[i] = tuple(Fraction(v) for v in xy)
    if not groups:
        raise ValueError("Empty coordinate file")
    summaries, triangle_rows = [], []
    for n, points in sorted(groups.items()):
        if set(points) != set(range(1, n+1)):
            raise ValueError(f"Expected exactly {n} points for n={n}")
        inside = all(0 <= c <= 1 for p in points.values() for c in p)
        areas = []
        for i, j, k in combinations(range(1, n+1), 3):
            x, y = points[i]; u, v = points[j]; a, b = points[k]
            area = abs(x*(v-b) + u*(b-y) + a*(y-v))/2
            areas.append((area, (i, j, k)))
            triangle_rows.append([n, i, j, k, str(area), decimal_string(area)])
        area, triple = min(areas)
        floor_scaled = area.numerator * 10**7 // area.denominator
        lower = f"{floor_scaled // 10**7}.{floor_scaled % 10**7:07d}" if inside else None
        summaries.append({
            "n": n, "points_in_unit_square": inside,
            "minimum_area_exact_fraction": str(area),
            "minimum_area_exact_decimal": decimal_string(area),
            "one_minimizing_triple": list(triple),
            "triangle_count": len(areas),
            "feasible_lower_bound_7dp": lower,
        })
    return summaries, triangle_rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("coordinates", nargs="?", type=Path,
                        default=Path(__file__).resolve().parents[1] / "data" / "coordinates" / "appendix_v15.csv")
    parser.add_argument("--output", type=Path, default=Path("validation_results"))
    args = parser.parse_args()
    results, triangles = check_file(args.coordinates)
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output/"validation.json").write_text(json.dumps(results, indent=2), encoding="utf-8")
    with (args.output/"all_triangle_areas.csv").open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["n", "i", "j", "k", "area_exact_fraction", "area_exact_decimal"])
        w.writerows(triangles)
    for row in results:
        print(f"n={row['n']}: A(P)={row['minimum_area_exact_decimal']}; "
              f"inside_square={row['points_in_unit_square']}; "
              f"lower_bound_7dp={row['feasible_lower_bound_7dp']}")
    print(f"Results: {args.output.resolve()}")
    return 0 if all(r["points_in_unit_square"] for r in results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
