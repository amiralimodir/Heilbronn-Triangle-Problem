"""Random feasible configurations for input lower bounds, n=3,...,12.

The supplied sampling order is retained: left, right, left, right, bottom,
bottom, top, top, followed by exactly n-8 interior points when n>8.
This is a lower-bound heuristic, not an optimality or upper-bound computation.
"""
import argparse
import csv
from datetime import datetime
from decimal import Decimal, localcontext
from fractions import Fraction
import json
import math
from pathlib import Path
import random
import sys
import time
import uuid


def mindist(xy, x2, y2):
    distance = math.sqrt(2)
    for x1, y1 in xy:
        distance = min(distance, ((x1 - x2)**2 + (y1 - y2)**2)**0.5)
    return distance


def sample_configuration(n, rng, min_distance, max_attempts):
    # As in the input code, the first two points are drawn without rejection.
    xy = [(0.0, rng.uniform(0, 1)), (1.0, rng.uniform(0, 1))]
    for point in range(2, n):
        accepted = False
        for attempt in range(max_attempts):
            if point == 2:
                x2, y2 = 0.0, rng.uniform(0, 1)
            elif point == 3:
                x2, y2 = 1.0, rng.uniform(0, 1)
            elif point in (4, 5):
                x2, y2 = rng.uniform(0, 1), 0.0
            elif point in (6, 7):
                x2, y2 = rng.uniform(0, 1), 1.0
            else:
                x2, y2 = rng.uniform(0, 1), rng.uniform(0, 1)
            if mindist(xy, x2, y2) >= min_distance:
                xy.append((x2, y2))
                accepted = True
                break
        if not accepted:
            return None
    return xy


def minimum_area(xy):
    minimum = 1.0
    n = len(xy)
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                area = abs(0.5 * (
                    xy[i][0] * (xy[j][1] - xy[k][1])
                    + xy[j][0] * (xy[k][1] - xy[i][1])
                    + xy[k][0] * (xy[i][1] - xy[j][1])))
                minimum = min(minimum, area)
    return minimum


def verify_configuration(xy):
    # JSON writes these same decimal representations. Treat them as exact
    # rationals when computing the reusable feasible lower bound.
    points = [(Fraction(str(x)), Fraction(str(y))) for x, y in xy]
    if not all(0 <= value <= 1 for point in points for value in point):
        raise ValueError('A generated coordinate is outside the unit square.')
    minimum, triangle = None, None
    for i in range(len(points)):
        for j in range(i + 1, len(points)):
            for k in range(j + 1, len(points)):
                area = abs((
                    points[i][0] * (points[j][1] - points[k][1])
                    + points[j][0] * (points[k][1] - points[i][1])
                    + points[k][0] * (points[i][1] - points[j][1])) / 2)
                if minimum is None or area < minimum:
                    minimum, triangle = area, (i, j, k)
    with localcontext() as context:
        context.prec = len(str(minimum.numerator)) + len(str(minimum.denominator)) + 20
        decimal_area = format(Decimal(minimum.numerator) / Decimal(minimum.denominator), 'f')
    # Integer division floors exactly, including when the area is zero.
    scaled = minimum.numerator * 10**7 // minimum.denominator
    lower_bound = '{}.{:07d}'.format(scaled // 10**7, scaled % 10**7)
    return {
        'coordinates_in_unit_square': True,
        'minimum_area_exact_fraction': str(minimum),
        'minimum_area_exact_decimal': decimal_area,
        'minimum_triangle_zero_based': list(triangle),
        'feasible_lower_bound_7dp': lower_bound,
    }


def write_json(path, values):
    temporary = path.with_name(path.name + '.tmp')
    temporary.write_text(json.dumps(values, indent=2, allow_nan=False), encoding='utf-8')
    temporary.replace(path)


def run_instance(n, args, seed, output_dir):
    rng = random.Random(seed)
    started = time.perf_counter()
    best_value, best_xy = -1.0, None
    accepted, rejected = 0, 0
    config = {
        'n': n, 'trials': args.trials, 'seed': seed,
        'minimum_distance_filter': args.min_distance,
        'max_sampling_attempts_per_point': args.max_attempts,
        'layout': 'left,right,left,right,bottom,bottom,top,top; remaining points interior',
        'python_version': sys.version,
        'purpose': 'feasible lower-bound heuristic',
    }
    write_json(output_dir / 'run_config.json', config)
    progress_interval = max(1, args.trials // 10)
    for trial in range(args.trials):
        xy = sample_configuration(n, rng, args.min_distance, args.max_attempts)
        if xy is None:
            rejected += 1
        else:
            accepted += 1
            value = minimum_area(xy)
            if value > best_value:
                best_value, best_xy = value, xy
                # Save the actual witness whenever the best sampled value improves.
                record = dict(config)
                record.update({
                    'x': [point[0] for point in best_xy],
                    'y': [point[1] for point in best_xy],
                    'sampled_minimum_area_float': best_value,
                    'best_trial_one_based': trial + 1,
                    'completed_trials': trial + 1,
                    'status': 'RUNNING',
                })
                record.update(verify_configuration(best_xy))
                write_json(output_dir / 'result.json', record)
        if not args.quiet and (trial + 1) % progress_interval == 0:
            print('n={} trials={}/{} best={}'.format(
                n, trial + 1, args.trials, best_value if best_xy is not None else None), flush=True)
    if best_xy is None:
        record = dict(config, x=None, y=None, feasible_lower_bound_7dp=None,
                      status='NO_CONFIGURATION')
    else:
        record['status'] = 'COMPLETED'
    record.update({
        'completed_trials': args.trials, 'accepted_configurations': accepted,
        'rejected_configurations': rejected, 'runtime_seconds': time.perf_counter() - started,
    })
    write_json(output_dir / 'result.json', record)
    if best_xy is not None:
        with (output_dir / 'coordinates.csv').open('w', newline='', encoding='utf-8') as stream:
            writer = csv.writer(stream)
            writer.writerow(('n', 'i', 'x', 'y'))
            for index, (x, y) in enumerate(best_xy, 1):
                writer.writerow((n, index, x, y))
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__, allow_abbrev=False)
    parser.add_argument('--n', nargs='+', type=int, choices=range(3, 13), default=list(range(3, 13)))
    parser.add_argument('--trials', type=int, default=500000, help='Attempted configurations per n; default: 500000.')
    parser.add_argument('--seed', type=int, help='Base seed; each n uses base seed + n. Saved if generated automatically.')
    parser.add_argument('--min-distance', type=float, default=0.167718)
    parser.add_argument('--max-attempts', type=int, default=10000,
                        help='Maximum rejection draws per point before discarding that trial.')
    parser.add_argument('--output', type=Path, default=Path(__file__).resolve().parent / 'results' / 'max_random')
    parser.add_argument('--quiet', action='store_true')
    parser.add_argument('--list', action='store_true', help='List selected point counts without sampling.')
    args = parser.parse_args()
    if args.trials < 1 or args.max_attempts < 1:
        parser.error('Trials and maximum attempts must be positive integers.')
    if not math.isfinite(args.min_distance) or not 0 <= args.min_distance <= math.sqrt(2):
        parser.error('Minimum distance must be finite and between 0 and sqrt(2).')
    counts = list(dict.fromkeys(args.n))
    if args.list:
        print(json.dumps([{'n': n, 'trials': args.trials} for n in counts], indent=2))
        return 0
    base_seed = args.seed if args.seed is not None else random.SystemRandom().randrange(2**63)
    batch = args.output.expanduser().resolve() / (
        datetime.now().strftime('%Y%m%d_%H%M%S_%f') + '_' + uuid.uuid4().hex[:8])
    batch.mkdir(parents=True, exist_ok=False)
    write_json(batch / 'settings.json', {
        'base_seed': base_seed, 'seed_rule': 'base_seed + n', 'n': counts,
        'trials_per_n': args.trials, 'minimum_distance_filter': args.min_distance,
        'max_sampling_attempts_per_point': args.max_attempts,
    })
    rows = []
    for n in counts:
        output_dir = batch / ('n' + str(n))
        output_dir.mkdir()
        record = run_instance(n, args, base_seed + n, output_dir)
        rows.append({key: record[key] for key in (
            'n', 'seed', 'feasible_lower_bound_7dp', 'runtime_seconds',
            'completed_trials', 'accepted_configurations', 'rejected_configurations', 'status')})
        write_json(batch / 'summary.json', rows)
        with (batch / 'summary.csv').open('w', newline='', encoding='utf-8') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        print('n={} feasible lower bound={} seed={}'.format(
            n, record['feasible_lower_bound_7dp'], base_seed + n), flush=True)
    print('Results:', batch)
    return 1 if any(row['status'] != 'COMPLETED' for row in rows) else 0


if __name__ == '__main__':
    raise SystemExit(main())
