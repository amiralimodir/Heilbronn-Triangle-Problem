"""Auxiliary Heilbronn problems in [0,1] x [0,1/M].

The model uses only y-ordering, coordinate bounds, product equalities,
and the area formulation. It does not fix points on rectangle edges.
Threshold decisions use ObjBound, never the incumbent objective.
"""
import argparse
import csv
from datetime import datetime
from decimal import Decimal, InvalidOperation
from fractions import Fraction
import json
import math
import os
from pathlib import Path
import platform
import shutil
import sys
import time
import uuid


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--n', nargs='+', type=int, default=list(range(3, 10)),
                        help='Auxiliary point counts; default: 3 through 9.')
    parser.add_argument('--m', nargs='+', type=int, default=list(range(1, 8)),
                        help='Inverse strip heights; default: 1 through 7.')
    parser.add_argument('--threshold', help='Known feasible lower bound for the original instance; decimal value.')
    parser.add_argument('--target-n', type=int,
                        help='Original point count for deriving ordered y-bounds.')
    parser.add_argument('--time-limit', type=float,
                        help='Optional per-run seconds; omitted means the solver default.')
    parser.add_argument('--output', type=Path, default=Path(__file__).resolve().parent/'results')
    parser.add_argument('--list', action='store_true', help='List planned cases without loading Gurobi.')
    args = parser.parse_args()
    if any(n < 3 for n in args.n) or any(m < 1 for m in args.m):
        parser.error('Each n must be at least 3 and each M at least 1.')
    args.n = sorted(set(args.n))
    args.m = sorted(set(args.m))
    if args.time_limit is not None and (not math.isfinite(args.time_limit) or args.time_limit <= 0):
        parser.error('The time limit must be finite and positive.')
    if args.threshold is not None:
        try:
            threshold = Decimal(args.threshold)
            if not threshold.is_finite() or not 0 < threshold <= Decimal('0.5'):
                raise ValueError
        except (InvalidOperation, ValueError):
            parser.error('The threshold must be a positive finite decimal no greater than 0.5.')
    if args.target_n is not None:
        if args.target_n < 3 or args.threshold is None:
            parser.error('--target-n requires n>=3 and --threshold.')
    return args


def heilbronn_triangle_approach1(n, M, ub, lb):
    """Build the procedural formulation; the caller records and solves it."""
    model = gp.Model('Heilbronn Triangle')
    x = model.addVars(n, vtype=GRB.CONTINUOUS, name='x', lb=0, ub=1)
    y = model.addVars(n, vtype=GRB.CONTINUOUS, name='y', lb=0, ub=1/M)
    w = model.addVars(n, n, vtype=GRB.CONTINUOUS, name='w', lb=0, ub=1)
    S = model.addVars(n, n, n, vtype=GRB.CONTINUOUS, name='S', lb=-0.5, ub=0.5)
    b = model.addVars(n, n, n, vtype=GRB.BINARY, name='b')
    z = model.addVar(vtype=GRB.CONTINUOUS, name='z', lb=lb, ub=ub)
    model.update()

    # Label ordering preserves every configuration. No boundary is fixed.
    for i in range(n-1):
        model.addConstr(y[i] <= y[i+1], name='Sort points')
    for i in range(n):
        for j in range(n):
            model.addConstr(w[i,j] == x[i]*y[j])

    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                model.addConstr(S[i,j,k] == 0.5*(w[i,j]-w[i,k]+w[j,k]-w[j,i]+w[k,i]-w[k,j]), name=f'S_constr_{i}_{j}_{k}')
                model.addConstr((1-b[i,j,k])*(ub+0.5)+S[i,j,k] >= z, name=f'linearize1_{i}_{j}_{k}')
                model.addConstr(b[i,j,k]*(ub+0.5)-S[i,j,k] >= z, name=f'linearize2_{i}_{j}_{k}')
                model.addConstr((-1)/2 <= S[i,j,k]-((1/2+lb)*b[i,j,k]))
                model.addConstr(S[i,j,k]-((1/2+lb)*b[i,j,k]) <= -lb)
                model.addConstr(S[i,j,k] <= 0.5*b[i,j,k], name='upper')
                model.addConstr(S[i,j,k] >= 0.5*(b[i,j,k]-1), name='lower')
    model.setObjective(z, GRB.MAXIMIZE)
    return model, x, y, b


def finite_attribute(model, name):
    try:
        value = float(model.getAttr(name))
        return value if math.isfinite(value) else None
    except (AttributeError, gp.GurobiError):
        return None


def write_json(path, value):
    path = Path(path)
    temporary = path.with_name(path.name+'.tmp')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False), encoding='utf-8')
    os.replace(str(temporary), str(path))


def threshold_crossed(record, threshold):
    """Strict numerical dual comparison; errors/infeasibility are not capacities."""
    if threshold is None or record.get('error'):
        return False
    usable = {'OPTIMAL', 'TIME_LIMIT', 'ITERATION_LIMIT', 'NODE_LIMIT',
              'SOLUTION_LIMIT', 'INTERRUPTED', 'USER_OBJ_LIMIT',
              'WORK_LIMIT', 'MEM_LIMIT'}
    if record.get('status') not in usable:
        return False
    bound = record.get('ObjBound')
    if (isinstance(bound, bool) or not isinstance(bound, (int, float))
            or not math.isfinite(bound) or bound < 0):
        return False
    return Decimal.from_float(float(bound)) < Decimal(str(threshold))


def derive_y_bounds(records, target_n, threshold):
    lower = [Fraction(0) for _ in range(target_n)]
    upper = [Fraction(1) for _ in range(target_n)]
    capacities = []
    for record in records:
        r, M = record['n'], record['M']
        if (record.get('model') != 'unrestricted' or M < 2 or not 3 <= r <= target_n
                or not threshold_crossed(record, threshold)):
            continue
        capacity = r-1
        capacities.append({'M': M, 'capacity': capacity, 'excluded_point_count': r,
                           'ObjBound': record['ObjBound'], 'result_file': record.get('result_file')})
        for ell in range(1, M):
            if ell*capacity >= target_n:
                break
            lower[ell*capacity] = max(lower[ell*capacity], Fraction(ell, M))
            upper[target_n-ell*capacity-1] = min(upper[target_n-ell*capacity-1], 1-Fraction(ell, M))
    # Propagate the independently proved inequalities through y-ordering.
    for i in range(1, target_n):
        lower[i] = max(lower[i], lower[i-1])
    for i in range(target_n-2, -1, -1):
        upper[i] = min(upper[i], upper[i+1])
    return {'target_n': target_n, 'threshold': threshold, 'capacities': capacities,
            'indexing': 'lists use zero-based indices; coordinate i+1 is at list index i',
            'lower_bounds': [str(value) for value in lower],
            'upper_bounds': [str(value) for value in upper],
            'consistent': all(lo <= hi for lo, hi in zip(lower, upper)),
            'qualification': 'Numerical full-rectangle dual comparisons against the supplied feasible lower bound; no independent exact certification of solver duals. No Group-2 endpoint equalities are added.'}


def solve_rectangle(n, M, args, run_dir):
    run_dir.mkdir(parents=True, exist_ok=False)
    source = Path(__file__).resolve()
    shutil.copyfile(str(source), str(run_dir/'executed_source.py'))
    model = None
    previous = Path.cwd()
    record = {'n': n, 'M': M, 'model': 'unrestricted',
              'input_lower_bound': 0, 'input_upper_bound': 1/(2*M),
              'rectangle_width': 1, 'rectangle_height': 1/M,
              'threshold': args.threshold, 'target_n': args.target_n,
              'scope': 'full_rectangle_up_to_label_order',
              'status': 'ERROR', 'ObjBound': None, 'ObjVal': None, 'x': None, 'y': None,
              'solver_runtime_seconds': None, 'elapsed_optimize_seconds': None,
              'python_version': platform.python_version(), 'gurobi_version': list(gp.gurobi.version()),
              'result_file': str(run_dir/'result.json')}
    try:
        os.chdir(str(run_dir))
        model, x, y, b = heilbronn_triangle_approach1(n, M, 1/(2*M), 0)
        model.setParam('LogFile', 'solver.log')
        if args.time_limit is not None:
            model.setParam('TimeLimit', args.time_limit)
        model.update()
        model.write('model.mps')
        model.write('parameters.prm')
        record['effective_parameters'] = {key: model.getParamInfo(key)[2] for key in (
            'TimeLimit', 'MIPGap', 'MIPGapAbs', 'Threads', 'Seed', 'NonConvex',
            'Symmetry', 'FeasibilityTol', 'IntFeasTol', 'OptimalityTol')}
        # The default time limit is infinity, which is not a JSON number.
        record['effective_parameters'] = {key: (value if not isinstance(value, float) or math.isfinite(value) else str(value))
                                          for key, value in record['effective_parameters'].items()}
        write_json(run_dir/'run_config.json', record)
        started = time.perf_counter()
        model.optimize()
        record['elapsed_optimize_seconds'] = time.perf_counter()-started
        names = ('LOADED', 'OPTIMAL', 'INFEASIBLE', 'INF_OR_UNBD', 'UNBOUNDED',
                 'CUTOFF', 'ITERATION_LIMIT', 'NODE_LIMIT', 'TIME_LIMIT',
                 'SOLUTION_LIMIT', 'INTERRUPTED', 'NUMERIC', 'SUBOPTIMAL',
                 'INPROGRESS', 'USER_OBJ_LIMIT', 'WORK_LIMIT', 'MEM_LIMIT')
        statuses = {getattr(GRB, key): key for key in names if hasattr(GRB, key)}
        record['status_code'] = int(model.Status)
        record['status'] = statuses.get(model.Status, str(model.Status))
        record['solution_count'] = int(model.SolCount)
        record['ObjBound'] = finite_attribute(model, 'ObjBound')
        record['solver_runtime_seconds'] = finite_attribute(model, 'Runtime')
        if model.SolCount:
            record['ObjVal'] = finite_attribute(model, 'ObjVal')
            record['MIPGap'] = finite_attribute(model, 'MIPGap')
            record['x'] = [x[i].X for i in range(n)]
            record['y'] = [y[i].X for i in range(n)]
            model.write('solution.sol')
    except Exception as error:
        record['error'] = repr(error)
    finally:
        try:
            if model is not None:
                model.dispose()
        except Exception as error:
            record['error'] = record.get('error', repr(error))
        finally:
            os.chdir(str(previous))
    record['threshold_crossed'] = threshold_crossed(record, args.threshold)
    record['numerical_capacity'] = (n-1 if record['threshold_crossed'] and M >= 2 else None)
    write_json(run_dir/'result.json', record)
    return record


def main(args):
    global gp, GRB
    import gurobipy as gp
    from gurobipy import GRB
    session = args.output.expanduser().resolve()/('rectangles_'+datetime.now().strftime('%Y%m%d_%H%M%S_%f')+'_'+uuid.uuid4().hex[:8])
    session.mkdir(parents=True, exist_ok=False)
    records = []
    for M in args.m:
        for n in args.n:
            run_dir = session/('M{}_n{}'.format(M, n))
            record = solve_rectangle(n, M, args, run_dir)
            records.append(record)
            write_json(session/'rectangle_results.json', records)
            fields = ('n', 'M', 'model', 'scope', 'status', 'ObjBound', 'ObjVal',
                      'solver_runtime_seconds', 'elapsed_optimize_seconds',
                      'threshold_crossed', 'numerical_capacity', 'result_file', 'error')
            with (session/'rectangle_table.csv').open('w', newline='', encoding='utf-8') as stream:
                writer = csv.DictWriter(stream, fieldnames=fields, extrasaction='ignore')
                writer.writeheader()
                writer.writerows(records)
            if args.target_n is not None:
                write_json(session/'y_bounds.json', derive_y_bounds(records, args.target_n, args.threshold))
            print('M={} n={} status={} ObjBound={} ObjVal={}'.format(M, n, record['status'], record['ObjBound'], record['ObjVal']), flush=True)
            if record['threshold_crossed']:
                print('Threshold crossed for ' + record['scope'] + '; stopping this M.', flush=True)
                break
    print('Results: ' + str(session))
    return 1 if any(record.get('error') for record in records) else 0


if __name__ == '__main__':
    arguments = parse_args()
    if arguments.list:
        for m in arguments.m:
            for n in arguments.n:
                print(json.dumps({'n': n, 'M': m, 'model': 'unrestricted',
                                  'threshold': arguments.threshold, 'time_limit': arguments.time_limit}))
        raise SystemExit(0)
    raise SystemExit(main(arguments))
