"""Check saved coordinates exactly, without a solver.

Accepts CSV (n,i,x,y), JSON containing x/y arrays, or a results directory.
Decimal coordinates are treated as exact rational numbers. This verifies
feasible configuration values, not optimality or numerical dual bounds.
"""
import argparse
import csv
from decimal import Decimal, InvalidOperation, localcontext
from fractions import Fraction
from itertools import combinations
import json
from pathlib import Path


def decimal_string(q):
    with localcontext() as ctx:
        ctx.prec = len(str(abs(q.numerator))) + len(str(q.denominator)) + 20
        return format(Decimal(q.numerator) / Decimal(q.denominator), "f")


def coordinate(value):
    if isinstance(value, bool):
        raise ValueError("Coordinates must be finite decimal numbers")
    try:
        number = Decimal(str(value))
    except InvalidOperation as error:
        raise ValueError("Coordinates must be finite decimal numbers") from error
    if not number.is_finite():
        raise ValueError("Coordinates must be finite decimal numbers")
    return Fraction(number)


def read_csv(path):
    groups = {}
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        if not {"n", "i", "x", "y"}.issubset(reader.fieldnames or []):
            raise ValueError("CSV must contain n,i,x,y columns")
        for row in reader:
            n, i = int(row["n"]), int(row["i"])
            if n < 3 or not 1 <= i <= n:
                raise ValueError("Invalid index n={}, i={}".format(n, i))
            group = groups.setdefault(n, {})
            if i in group:
                raise ValueError("Duplicate point index n={}, i={}".format(n, i))
            group[i] = tuple(coordinate(row[c]) for c in ("x", "y"))
    if not groups:
        raise ValueError("Empty coordinate file")
    return groups


def read_json(path):
    # parse_float preserves all decimal digits in the saved file.
    data = json.loads(path.read_text(encoding="utf-8-sig"), parse_float=Decimal)
    configurations = []

    def visit(value, label):
        if isinstance(value, dict):
            if "x" in value and "y" in value:
                xs, ys = value["x"], value["y"]
                if xs is None or ys is None:
                    return  # A run without an incumbent has no configuration.
                if not isinstance(xs, list) or not isinstance(ys, list):
                    raise ValueError("{}: x and y must be arrays".format(label))
                n = len(xs)
                if n < 3 or len(ys) != n:
                    raise ValueError("{}: expected matching x/y arrays of length >=3".format(label))
                if "n" in value and int(value["n"]) != n:
                    raise ValueError("{}: n does not match the coordinates".format(label))
                points = {i + 1: (coordinate(x), coordinate(y))
                          for i, (x, y) in enumerate(zip(xs, ys))}
                configurations.append((label, {n: points}))
                return
            for key, child in value.items():
                visit(child, label + "/" + str(key))
        elif isinstance(value, list):
            for i, child in enumerate(value):
                visit(child, label + "/" + str(i))

    visit(data, str(path))
    return configurations


def evaluate(groups, source=None):
    summaries, triangle_rows = [], []
    for n, points in sorted(groups.items()):
        if set(points) != set(range(1, n + 1)):
            raise ValueError("Expected exactly {} points for n={}".format(n, n))
        inside = all(0 <= c <= 1 for p in points.values() for c in p)
        areas = []
        for i, j, k in combinations(range(1, n + 1), 3):
            x, y = points[i]; u, v = points[j]; a, b = points[k]
            area = abs(x * (v - b) + u * (b - y) + a * (y - v)) / 2
            areas.append((area, (i, j, k)))
            triangle_rows.append([n, i, j, k, str(area), decimal_string(area)])
        area, triple = min(areas)
        floor_scaled = area.numerator * 10**7 // area.denominator
        lower = ("{}.{:07d}".format(floor_scaled // 10**7, floor_scaled % 10**7)
                 if inside else None)
        summaries.append({
            "source": source, "n": n, "points_in_unit_square": inside,
            "minimum_area_exact_fraction": str(area),
            "minimum_area_exact_decimal": decimal_string(area),
            "one_minimizing_triple": list(triple), "triangle_count": len(areas),
            "feasible_lower_bound_7dp": lower,
        })
    return summaries, triangle_rows


def check_file(path):
    """Validate a CSV file; retained for simple programmatic use."""
    path = Path(path)
    return evaluate(read_csv(path), str(path))


def check_inputs(paths):
    configurations = []
    skipped = []
    seen = set()
    for requested in paths:
        requested = Path(requested)
        if requested.is_dir():
            # Prefer individual run files; do not count aggregate summaries twice.
            files = sorted(requested.rglob("result.json"))
            if not files:
                files = sorted(requested.rglob("*.json")) + sorted(requested.rglob("*.csv"))
        else:
            files = [requested]
        for path in files:
            resolved = path.resolve()
            if resolved in seen:
                continue
            seen.add(resolved)
            if path.suffix.lower() == ".csv":
                configurations.append((str(path), read_csv(path)))
            elif path.suffix.lower() == ".json":
                found = read_json(path)
                configurations.extend(found)
                if not found:
                    skipped.append(str(path))
            else:
                raise ValueError("Expected a CSV/JSON file or results directory: {}".format(path))
    if not configurations:
        raise ValueError("No saved coordinate configurations found in the supplied inputs")
    summaries, triangles = [], []
    for source, groups in configurations:
        result, rows = evaluate(groups, source)
        summaries.extend(result)
        triangles.extend([[source] + row for row in rows])
    return summaries, triangles, skipped


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("coordinates", nargs="+", type=Path,
                        help="CSV/JSON coordinate files or directories of saved run results.")
    parser.add_argument("--output", type=Path, default=Path("validation_results"))
    args = parser.parse_args()
    try:
        results, triangles, skipped = check_inputs(args.coordinates)
    except (ValueError, OSError) as error:
        parser.error(str(error))
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "validation.json").write_text(json.dumps(results, indent=2), encoding="utf-8")
    with (args.output / "all_triangle_areas.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["source", "n", "i", "j", "k", "area_exact_fraction", "area_exact_decimal"])
        writer.writerows(triangles)
    for row in results:
        print("{} | n={}: A(P)={}; inside_square={}; lower_bound_7dp={}".format(
            row["source"], row["n"], row["minimum_area_exact_decimal"],
            row["points_in_unit_square"], row["feasible_lower_bound_7dp"]))
    if skipped:
        print("Skipped {} JSON files without saved coordinates.".format(len(skipped)))
    print("Results: {}".format(args.output.resolve()))
    return 0 if all(r["points_in_unit_square"] for r in results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
