"""Complete solver-comparison table: four solvers, n=3,...,7.

LB=0, UB=0.5, requested time limit=86400; no enhancement groups.
Gurobi uses the native Approach-1 model; the other solvers retain the supplied
Pyomo model. See README.md for the formulation and timing distinctions.
"""

import argparse

SOLVERS = ("gurobi", "baron", "scip", "couenne")
PENDING_DUAL_RUNS = [
    ("couenne", 3), ("couenne", 4), ("couenne", 5),
    ("couenne", 6), ("couenne", 7), ("scip", 7),
]
ALL_RUNS = [(solver, n) for n in range(3, 8) for solver in SOLVERS]


def select_runs(args):
    selected = PENDING_DUAL_RUNS if args.missing_duals else ALL_RUNS
    return [(solver, n) for solver, n in selected
            if solver in args.solvers and n in args.n]

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run all 20 solver-table cases by default, or select a subset.",
        allow_abbrev=False)
    schedule = parser.add_mutually_exclusive_group()
    schedule.add_argument("--all", action="store_true", help="Run the full table (the default).")
    schedule.add_argument("--missing-duals", action="store_true",
                          help="Select only Couenne n=3..7 and SCIP n=7.")
    parser.add_argument("--solvers", nargs="+", choices=SOLVERS, default=list(SOLVERS))
    parser.add_argument("--n", nargs="+", type=int, choices=range(3, 8), default=list(range(3, 8)))
    parser.add_argument("--output", help="Writable ASCII-only output folder.")
    parser.add_argument("--stop-on-error", action="store_true",
                        help="Stop after a failed solve or missing dual bound; otherwise continue.")
    parser.add_argument("--list", action="store_true", help="Print selected cases without loading solvers.")
    CLI_ARGS = parser.parse_args()
    if CLI_ARGS.list:
        cases = select_runs(CLI_ARGS)
        for solver, points in cases:
            print(solver, points)
        raise SystemExit(0)

import time
import json
import csv
import math
import os
import re
import sys
import shutil
import struct
from pathlib import Path
from datetime import datetime
from decimal import Decimal
from importlib.metadata import version as package_version


TIME_LIMIT = 86400
SCIP_EXECUTABLE = r"C:\Program Files\SCIPOptSuite 10.1.0\bin\scip.exe"
COUENNE_EXECUTABLE = r"C:\Program Files\Couenne\ampl.mswin64\couenne.exe"
BARON_EXECUTABLE = r"C:\baron\baron.exe"
OUTPUT_ROOT = Path.home() / "heilbronn_solver_comparison"

# These runs deliberately retain the supplied model, including the absence
# of the extra sign-strengthening inequalities. LB=0 and UB=0.5 are unchanged.
# No feasibility, gap, or thread settings are changed.
SESSION_DIR = None
LAST_RUN_DETAILS = {}
SCRIPT_PATH = Path(__file__).resolve()
NUMBER = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eEdD][+-]?\d+)?"
BOUND_TOKEN = r"(?:" + NUMBER + r"|[+-]?(?:inf(?:inity)?|nan))"


def finite_number(value):
    try:
        number = float(str(value).replace("D", "E").replace("d", "e"))
        return number if math.isfinite(number) else None
    except (TypeError, ValueError, OverflowError):
        return None


def write_json(path, data):
    # Atomic checkpoint; an interruption cannot replace it with half a JSON file.
    path = Path(path)
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        json.dump(data, stream, indent=4, allow_nan=False)
    os.replace(str(temporary), str(path))


def printed_rounding_allowance(token, significant_digits):
    value = Decimal(token.replace("D", "E").replace("d", "e"))
    if not value.is_finite() or value.is_zero():
        return 0.0
    # Half the last printed significant digit, for the solver's %g output.
    return float(Decimal(5).scaleb(value.copy_abs().adjusted() - significant_digits))


def float_rounded_up(decimal_value):
    """Convert a finite Decimal outward, including on Python 3.8."""
    value = float(decimal_value)
    if Decimal.from_float(value) < decimal_value:
        if value == 0.0:
            return float.fromhex("0x0.0000000000001p-1022")
        bits = struct.unpack("!Q", struct.pack("!d", value))[0]
        bits += 1 if value > 0 else -1
        value = struct.unpack("!d", struct.pack("!Q", bits))[0]
    return value


def extract_couenne_bound(log_text):
    """Couenne 0.5.8: original maximize(z) is internally minimize(-z).

    The final Lower bound is the global bound for that MINIMIZATION problem.
    Hence the original maximization upper bound is its NEGATIVE, not abs()
    and not the final internal Upper bound (which is the incumbent).
    Source: Couenne release 0.5.8, BonCouenne.cpp and CouenneProblem.cpp,
    commit 7154f7a9b3cd84be378d02b483d090b76fc79ce8.
    """
    pattern = (
        r"^[ \t]*Lower bound:[ \t]*(" + BOUND_TOKEN + r")[ \t]*\r?\n"
        r"[ \t]*Upper bound:[ \t]*(" + BOUND_TOKEN + r")[^\r\n]*\r?\n"
        r"[ \t]*Branch-and-bound nodes:[^\r\n]*"
    )
    matches = list(re.finditer(pattern, log_text, re.MULTILINE | re.IGNORECASE))
    if not matches:
        return None, {"bound_error": "Couenne final bound summary was not found."}
    match = matches[-1]
    lower_token, upper_token = match.group(1), match.group(2)
    lower = finite_number(lower_token)
    upper = finite_number(upper_token)
    info = {
        "bound_source": "couenne_final_summary: minus internal Lower bound",
        "bound_raw_text": match.group(0),
        "internal_lower_token": lower_token,
        "internal_upper_token": upper_token,
        "internal_minimization_lower_bound": lower,
        "internal_minimization_incumbent": upper,
        "bound_print_significant_digits": 6,
        "bound_error": None,
    }
    if lower is None or abs(lower) >= 1e40:
        info["bound_error"] = "Couenne did not print a usable finite global lower bound."
        return None, info
    allowance = printed_rounding_allowance(lower_token, 6)
    info["bound_print_rounding_allowance"] = allowance
    if upper is not None and abs(upper) < 1e40:
        upper_allowance = printed_rounding_allowance(upper_token, 6)
        if lower > upper + allowance + upper_allowance + 1e-10:
            info["bound_error"] = "Inconsistent internal lower/upper bounds in final summary."
            return None, info
    reported_decimal = -Decimal(lower_token.replace("D", "E").replace("d", "e"))
    exact_allowance = (Decimal(0) if reported_decimal.is_zero() else
                       Decimal(5).scaleb(reported_decimal.copy_abs().adjusted() - 6))
    adjusted_decimal = reported_decimal + exact_allowance
    dual_bound = float_rounded_up(adjusted_decimal)
    # The primary field includes an allowance for the six-significant-digit
    # printed summary, which can otherwise round a numerical UB DOWNWARD.
    # This accounts ONLY for print rounding, not solver numerical errors.
    info["dual_bound_reported"] = -lower
    info["dual_bound_decimal"] = str(adjusted_decimal)
    info["dual_bound_adjustment"] = "Half the final printed significant digit, outward conversion"
    info["dual_bound_with_print_rounding_allowance"] = dual_bound
    return dual_bound, info


def extract_dual_bound(solver_name, results, log_text):
    if solver_name == "couenne":
        return extract_couenne_bound(log_text)

    if solver_name == "scip":
        matches = list(re.finditer(
            r"^[ \t]*Dual Bound[ \t]*:[ \t]*(" + BOUND_TOKEN + r")",
            log_text, re.MULTILINE | re.IGNORECASE,
        ))
        if matches:
            token = matches[-1].group(1)
            bound = finite_number(token)
            if bound is not None and abs(bound) < 1e19:
                return bound, {
                    "bound_source": "scip_final_log_dual_bound",
                    "bound_raw_text": matches[-1].group(0),
                    "bound_error": None,
                }
        if results is not None:
            bound = finite_number(getattr(results.solver, "dual_bound", None))
            if bound is not None and abs(bound) < 1e19:
                return bound, {
                    "bound_source": "results.solver.dual_bound",
                    "bound_error": None,
                }

    if results is not None and len(results.problem) > 0:
        bound = finite_number(results.problem.upper_bound)
        if bound is not None and abs(bound) < 1e19:
            return bound, {
                "bound_source": "results.problem.upper_bound",
                "bound_error": None,
            }
    return None, {"bound_error": "No finite dual bound was recovered; inspect solver.log."}


def termination_label(results, log_text):
    condition = str(results.solver.termination_condition) if results is not None else "unknown"
    # Keep the original Pyomo status separately. In ASL, a time limit can appear
    # as maxIterations, so inspect the solver's own termination message too.
    message = str(results.solver.message) if results is not None else ""
    if condition in ("optimal", "globallyOptimal"):
        return "OPTIMAL_WITHIN_SOLVER_TOLERANCES"
    message_has_limit = re.search(
        r"time limit reached|maximum (?:cpu )?time|time limit exceeded|stopped on time",
        message, re.IGNORECASE,
    )
    global_log_has_limit = re.search(
        r"Cbc0020I[^\r\n]*maximum time|SCIP Status[^\r\n]*time limit reached",
        log_text, re.IGNORECASE,
    )
    if condition == "maxTimeLimit" or message_has_limit or global_log_has_limit:
        return "TIME_LIMIT"
    return condition


def solve_heilbronn(n, ub, lb, solver_name):
    """Same supplied mathematical model; corrected output recovery."""
    import pyomo.environ as pyo
    from pyomo.opt import SolverFactory
    from pyomo.common.tempfiles import TempfileManager

    global SESSION_DIR, LAST_RUN_DETAILS
    solver_name = solver_name.lower()
    if (lb, ub) != (0, 0.5):
        raise ValueError("The solver comparison is defined for LB=0 and UB=0.5.")
    if SESSION_DIR is None:
        SESSION_DIR = OUTPUT_ROOT / datetime.now().strftime("comparison_%Y%m%d_%H%M%S_%f")
        SESSION_DIR.mkdir(parents=True, exist_ok=False)
    run_dir = SESSION_DIR / (solver_name + "_n" + str(n))
    run_dir.mkdir(parents=True, exist_ok=False)
    run_dir = run_dir.resolve()
    if not str(run_dir).isascii():
        raise RuntimeError("Set OUTPUT_ROOT to a writable path containing ASCII characters only.")
    log_path = run_dir / "solver.log"
    shutil.copyfile(str(SCRIPT_PATH), str(run_dir / "executed_source.py"))
    LAST_RUN_DETAILS = {
        "solver": solver_name, "n": n, "input_lb": lb, "input_ub": ub,
        "time_limit": TIME_LIMIT, "run_directory": str(run_dir),
        "solver_log": str(log_path), "python_version": sys.version,
        "python_executable": sys.executable, "pyomo_version": package_version("pyomo"),
        "model_variant": "Pyomo Approach 1 without sign-strengthening inequalities",
        "enhancement_groups": [],
        "time_measure": "elapsed opt.solve call, including export and solution loading",
    }
    write_json(run_dir / "run_config.json", LAST_RUN_DETAILS)

    # --------------------------------------------------------------------------
    # 1. Pyomo Model Creation
    # --------------------------------------------------------------------------
    model = pyo.ConcreteModel(name=f"Heilbronn_Triangle_{solver_name}_{n}")

    # --------------------------------------------------------------------------
    # 2. Index Sets
    # --------------------------------------------------------------------------
    model.I = pyo.RangeSet(1, n)
    model.J = pyo.RangeSet(1, n)
    model.K = pyo.RangeSet(1, n)
    model.Triangles = pyo.Set(initialize=lambda m: [(i, j, k) for i in m.I for j in m.J for k in m.K if i < j < k])

    # --------------------------------------------------------------------------
    # 3. Decision Variables
    # --------------------------------------------------------------------------
    model.x = pyo.Var(model.I, bounds=(0, 1), within=pyo.NonNegativeReals)
    model.y = pyo.Var(model.I, bounds=(0, 1), within=pyo.NonNegativeReals)
    model.w = pyo.Var(model.I, model.J, bounds=(0, 1), within=pyo.NonNegativeReals)
    model.S = pyo.Var(model.Triangles, bounds=(-0.5, 0.5), within=pyo.Reals)
    model.b = pyo.Var(model.Triangles, within=pyo.Binary)
    model.z = pyo.Var(bounds=(0, 0.5), within=pyo.NonNegativeReals)

    # --------------------------------------------------------------------------
    # 4. Objective Function
    # --------------------------------------------------------------------------
    model.objective = pyo.Objective(expr=model.z, sense=pyo.maximize)

    # --------------------------------------------------------------------------
    # 5. Constraints -- unchanged from the supplied code
    # --------------------------------------------------------------------------
    model.w_def = pyo.Constraint(model.I, model.J, rule=lambda m, i, j: m.w[i, j] == m.x[i] * m.y[j])
    model.area_calculation = pyo.ConstraintList()
    for (i, j, k) in model.Triangles:
        expr = 0.5 * (model.w[i, j] - model.w[i, k] + model.w[j, k] - model.w[j, i] + model.w[k, i] - model.w[k, j])
        model.area_calculation.add(model.S[i, j, k] == expr)
        model.area_calculation.add((1 - model.b[i, j, k]) * (ub + 0.5) + model.S[i, j, k] >= model.z)
        model.area_calculation.add(model.b[i, j, k] * (ub + 0.5) - model.S[i, j, k] >= model.z)

    # --------------------------------------------------------------------------
    # 6. Dynamic Solver Factory Setup
    # --------------------------------------------------------------------------
    print(f"\nPyomo model created. Sending to {solver_name.upper()} solver...")
    if solver_name == "scip":
        opt = SolverFactory("scip", executable=SCIP_EXECUTABLE, solver_io="nl")
        opt.options["limits/time"] = TIME_LIMIT
    elif solver_name == "couenne":
        opt = SolverFactory("couenne", executable=COUENNE_EXECUTABLE)
        with (run_dir / "couenne.opt").open("w", encoding="ascii") as stream:
            stream.write("time_limit " + str(TIME_LIMIT) + "\n")
    elif solver_name == "baron":
        opt = SolverFactory("baron", executable=BARON_EXECUTABLE)
        opt.options["MaxTime"] = TIME_LIMIT
    else:
        raise ValueError(f"Unknown solver: {solver_name}")

    LAST_RUN_DETAILS["solver_options"] = dict(opt.options)
    if solver_name == "couenne":
        LAST_RUN_DETAILS["option_file"] = {"time_limit": TIME_LIMIT}
    write_json(run_dir / "run_config.json", LAST_RUN_DETAILS)

    old_cwd = Path.cwd()
    old_tempdir = TempfileManager.tempdir
    results = None
    solve_error = None
    load_error = None
    optimal_z, optimal_x, optimal_y = None, None, None
    try:
        # Isolate Couenne options and all .nl/.sol files in the ASCII run folder.
        os.chdir(str(run_dir))
        TempfileManager.tempdir = str(run_dir)
        start_time = time.time()
        try:
            results = opt.solve(
                model, tee=True, logfile=str(log_path),
                keepfiles=True, load_solutions=False,
            )
            invalid_conditions = {
                "infeasible", "unbounded", "infeasibleOrUnbounded",
                "invalidProblem", "internalSolverError", "solverFailure", "error",
            }
            condition = str(results.solver.termination_condition)
            if len(results.solution) > 0 and condition not in invalid_conditions:
                try:
                    # Includes warning/maxIterations solutions from the ASL reader.
                    model.solutions.load_from(results)
                    values = [pyo.value(model.z, exception=False)]
                    values += [pyo.value(model.x[i], exception=False) for i in model.I]
                    values += [pyo.value(model.y[i], exception=False) for i in model.I]
                    values = [finite_number(value) for value in values]
                    if all(value is not None for value in values):
                        optimal_z = values[0]
                        optimal_x = values[1:n + 1]
                        optimal_y = values[n + 1:]
                except Exception as exc:
                    load_error = repr(exc)
        except Exception as exc:
            solve_error = repr(exc)
        finally:
            # As in the original script, includes export, execution and loading.
            optimize_time = time.time() - start_time

        # Version probing is deliberately AFTER timing, not before the solve.
        try:
            reported_version = opt.version()
            solver_version = list(reported_version) if reported_version is not None else None
        except Exception as exc:
            solver_version = None
            LAST_RUN_DETAILS["version_error"] = repr(exc)
        LAST_RUN_DETAILS["solver_version"] = solver_version
        LAST_RUN_DETAILS["solver_executable"] = str(opt.executable())
        if results is not None:
            try:
                results.write(filename=str(run_dir / "pyomo_results.yml"), format="yaml")
            except Exception as exc:
                LAST_RUN_DETAILS["results_write_error"] = repr(exc)
    finally:
        os.chdir(str(old_cwd))
        TempfileManager.tempdir = old_tempdir

    # --------------------------------------------------------------------------
    # 7. Processing Results
    # --------------------------------------------------------------------------
    log_text = log_path.read_text(encoding="utf-8", errors="replace") if log_path.exists() else ""
    dual_bound, bound_info = extract_dual_bound(solver_name, results, log_text)
    condition = str(results.solver.termination_condition) if results is not None else "unknown"
    if condition in ("infeasible", "unbounded", "infeasibleOrUnbounded"):
        dual_bound = None
        bound_info["bound_error"] = "Unexpected termination for this comparison model: " + condition

    if solver_name == "couenne" and "internal_upper_token" in bound_info:
        upper = bound_info.get("internal_minimization_incumbent")
        if upper is None or abs(upper) >= 1e40:
            # Do not mistake initial/placeholder variable values for an incumbent.
            optimal_z, optimal_x, optimal_y = None, None, None

    if solver_name == "scip":
        primal_lines = re.findall(
            r"^[ \t]*Primal Bound[ \t]*:[ \t]*(" + BOUND_TOKEN + r")([^\r\n]*)",
            log_text, re.MULTILINE | re.IGNORECASE,
        )
        if primal_lines:
            primal_value = finite_number(primal_lines[-1][0])
            no_solutions = re.search(r"\b0 solutions\b", primal_lines[-1][1])
            if primal_value is None or abs(primal_value) >= 1e19 or no_solutions:
                optimal_z, optimal_x, optimal_y = None, None, None

    warnings = []
    if dual_bound is not None and optimal_z is not None and dual_bound < optimal_z:
        difference = optimal_z - dual_bound
        allowance = bound_info.get("bound_print_rounding_allowance", 0.0)
        warnings.append("Printed dual is below incumbent by " + str(difference)
                        + "; print-rounding allowance is " + str(allowance) + ".")
        if difference > allowance + 1e-5:
            bound_info["rejected_dual_bound"] = dual_bound
            bound_info["bound_error"] = "Bound/incumbent mismatch exceeds the diagnostic tolerance."
            dual_bound = None

    LAST_RUN_DETAILS.update(bound_info)
    LAST_RUN_DETAILS.update({
        "objective": optimal_z, "dual_bound": dual_bound,
        "x": optimal_x, "y": optimal_y, "time": optimize_time,
        "solver_status": str(results.solver.status) if results is not None else "error",
        "termination_condition": condition,
        "termination_status": termination_label(results, log_text),
        "solver_message": str(results.solver.message) if results is not None else None,
        "solve_error": solve_error, "solution_load_error": load_error,
        "warnings": warnings,
    })
    write_json(run_dir / "result.json", LAST_RUN_DETAILS)
    print("\n" + "=" * 50)
    print("Solver:", solver_name.upper(), "n:", n)
    print("Termination:", LAST_RUN_DETAILS["termination_status"])
    print("Pyomo termination:", condition)
    print("Objective:", optimal_z)
    print("Dual upper bound:", dual_bound)
    print("Bound source:", bound_info.get("bound_source"))
    print("Elapsed solve-call seconds:", optimize_time)
    print("Saved run folder:", run_dir)
    if bound_info.get("bound_error"):
        print("BOUND EXTRACTION ERROR:", bound_info["bound_error"])
    if solve_error:
        print("SOLVE ERROR:", solve_error)
    if load_error:
        print("SOLUTION LOAD ERROR:", load_error)
    for warning in warnings:
        print("WARNING:", warning)
    return optimal_z, dual_bound, optimal_x, optimal_y, optimize_time


def solve_gurobi(n, ub=0.5, lb=0):
    """Native Approach 1, without enhancements; solver defaults plus TimeLimit."""
    import gurobipy as gp
    from gurobipy import GRB

    global SESSION_DIR, LAST_RUN_DETAILS
    if (lb, ub) != (0, 0.5):
        raise ValueError("The solver comparison is defined for LB=0 and UB=0.5.")
    if SESSION_DIR is None:
        SESSION_DIR = OUTPUT_ROOT / datetime.now().strftime("comparison_%Y%m%d_%H%M%S_%f")
        SESSION_DIR.mkdir(parents=True, exist_ok=False)
    run_dir = (SESSION_DIR / ("gurobi_n" + str(n))).resolve()
    run_dir.mkdir(parents=True, exist_ok=False)
    if not str(run_dir).isascii():
        raise RuntimeError("Set OUTPUT_ROOT to a writable path containing ASCII characters only.")
    shutil.copyfile(str(SCRIPT_PATH), str(run_dir / "executed_source.py"))

    model = gp.Model("Heilbronn Triangle")
    x = model.addVars(n, vtype=GRB.CONTINUOUS, name="x", lb=0, ub=1)
    y = model.addVars(n, vtype=GRB.CONTINUOUS, name="y", lb=0, ub=1)
    w = model.addVars(n, n, vtype=GRB.CONTINUOUS, name="w", lb=0, ub=1)
    S = model.addVars(n, n, n, vtype=GRB.CONTINUOUS, name="S", lb=-0.5, ub=0.5)
    b = model.addVars(n, n, n, vtype=GRB.BINARY, name="b")
    z = model.addVar(vtype=GRB.CONTINUOUS, name="z", lb=lb, ub=ub)
    model.update()
    for i in range(n):
        for j in range(n):
            model.addConstr(w[i, j] == x[i] * y[j])
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                model.addConstr(S[i, j, k] == 0.5 * (
                    w[i, j] - w[i, k] + w[j, k] - w[j, i] + w[k, i] - w[k, j]),
                    name=f"S_constr_{i}_{j}_{k}")
                model.addConstr((1 - b[i, j, k]) * (ub + 0.5) + S[i, j, k] >= z,
                                name=f"linearize1_{i}_{j}_{k}")
                model.addConstr(b[i, j, k] * (ub + 0.5) - S[i, j, k] >= z,
                                name=f"linearize2_{i}_{j}_{k}")
                model.addConstr(-0.5 <= S[i, j, k] - (0.5 + lb) * b[i, j, k])
                model.addConstr(S[i, j, k] - (0.5 + lb) * b[i, j, k] <= -lb)
    model.setObjective(z, GRB.MAXIMIZE)
    model.setParam("TimeLimit", TIME_LIMIT)
    model.setParam("LogFile", str(run_dir / "solver.log"))
    model.update()
    model.write(str(run_dir / "model.mps"))
    model.write(str(run_dir / "parameters.prm"))
    settings = {name: model.getParamInfo(name)[2] for name in (
        "TimeLimit", "MIPGap", "MIPGapAbs", "NonConvex", "Threads", "Seed",
        "FeasibilityTol", "IntFeasTol", "OptimalityTol", "Symmetry")}
    LAST_RUN_DETAILS = {
        "solver": "gurobi", "n": n, "input_lb": lb, "input_ub": ub,
        "time_limit": TIME_LIMIT, "run_directory": str(run_dir),
        "solver_log": str(run_dir / "solver.log"), "python_version": sys.version,
        "python_executable": sys.executable, "solver_version": list(gp.gurobi.version()),
        "model_variant": "native Gurobi Approach 1 with sign-strengthening inequalities",
        "enhancement_groups": [], "effective_parameters": settings,
        "time_measure": "elapsed model.optimize call",
    }
    write_json(run_dir / "run_config.json", LAST_RUN_DETAILS)
    solve_error = None
    started = time.time()
    try:
        model.optimize()
    except Exception as exc:
        solve_error = repr(exc)
    solve_time = time.time() - started
    status_names = {
        GRB.OPTIMAL: "OPTIMAL_WITHIN_SOLVER_TOLERANCES", GRB.TIME_LIMIT: "TIME_LIMIT",
        GRB.INFEASIBLE: "INFEASIBLE", GRB.INF_OR_UNBD: "INFEASIBLE_OR_UNBOUNDED",
        GRB.UNBOUNDED: "UNBOUNDED", GRB.INTERRUPTED: "INTERRUPTED",
        GRB.NUMERIC: "NUMERIC", GRB.SUBOPTIMAL: "SUBOPTIMAL",
    }
    objective, dual_bound, x_values, y_values = None, None, None, None
    try:
        if model.SolCount > 0:
            objective = finite_number(model.ObjVal)
            x_values = [x[i].X for i in range(n)]
            y_values = [y[i].X for i in range(n)]
            model.write(str(run_dir / "solution.sol"))
        if model.Status not in (GRB.INFEASIBLE, GRB.INF_OR_UNBD, GRB.UNBOUNDED):
            dual_bound = finite_number(model.ObjBound)
    except Exception as exc:
        solve_error = solve_error or repr(exc)
    LAST_RUN_DETAILS.update({
        "objective": objective, "dual_bound": dual_bound,
        "x": x_values, "y": y_values, "time": solve_time,
        "solver_runtime": finite_number(model.Runtime),
        "solver_status": model.Status,
        "termination_status": status_names.get(model.Status, str(model.Status)),
        "bound_source": "Gurobi ObjBound" if dual_bound is not None else None,
        "bound_error": None if dual_bound is not None else "No finite ObjBound was available.",
        "solve_error": solve_error,
    })
    write_json(run_dir / "result.json", LAST_RUN_DETAILS)
    model.dispose()
    print("GUROBI n=", n, "dual=", dual_bound, "time=", solve_time,
          "status=", LAST_RUN_DETAILS["termination_status"])
    return objective, dual_bound, x_values, y_values, solve_time


def solve_with_models(n_points, solver_name):
    if solver_name == "gurobi":
        return solve_gurobi(n_points, 0.5, 0)
    return solve_heilbronn(n_points, 0.5, 0, solver_name)


def write_table(path, results):
    """Table-ready dual bounds and recorded seconds; never substitute objectives."""
    rows = []
    for solver in SOLVERS:
        for entry in results[solver]:
            for n, record in entry.items():
                rows.append({"n": int(n), "solver": solver,
                             "dual_bound": record.get("dual_bound"),
                             "time_seconds": record.get("time"),
                             "termination_status": record.get("termination_status")})
    rows.sort(key=lambda row: (row["n"], SOLVERS.index(row["solver"])))
    with Path(path).open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=(
            "n", "solver", "dual_bound", "time_seconds", "termination_status"))
        writer.writeheader()
        writer.writerows(rows)


# --------------------------------------------------------------------------
# Main Loop Execution
# --------------------------------------------------------------------------
if __name__ == "__main__":
    selected_runs = select_runs(CLI_ARGS)
    if not selected_runs:
        raise SystemExit("No cases match the selected filters.")
    if CLI_ARGS.output:
        OUTPUT_ROOT = Path(CLI_ARGS.output).expanduser()
    if not str(OUTPUT_ROOT.resolve()).isascii():
        raise SystemExit("Please set OUTPUT_ROOT to an ASCII-only writable folder.")
    for solver_name, n in selected_runs:
        if solver_name == "gurobi":
            try:
                import gurobipy
            except ImportError:
                raise SystemExit("Install gurobipy before selecting Gurobi.")
            continue
        try:
            import pyomo.environ
        except ImportError:
            raise SystemExit("Install Pyomo before selecting BARON, SCIP, or Couenne.")
        executable = {"couenne": COUENNE_EXECUTABLE, "scip": SCIP_EXECUTABLE, "baron": BARON_EXECUTABLE}[solver_name]
        if not Path(executable).is_file():
            raise SystemExit("Solver executable not found: " + executable)

    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    SESSION_DIR = OUTPUT_ROOT / datetime.now().strftime("comparison_%Y%m%d_%H%M%S_%f")
    SESSION_DIR.mkdir(exist_ok=False)
    output_filename = SESSION_DIR / "solver_results.json"
    result = {solver: [] for solver in SOLVERS}
    print("Selected cases:", selected_runs)
    print("Results will be saved after EVERY run to:", output_filename)

    for solver_name, n in selected_runs:
        print(f"\n-------------------- {solver_name.upper()}, n = {n} --------------------")
        try:
            z, dual_bound, x, y, solve_time = solve_with_models(n, solver_name)
        except Exception as exc:
            # Keep the rest of a long comparison batch runnable after one failure.
            run_dir = SESSION_DIR / (solver_name + "_n" + str(n))
            run_dir.mkdir(parents=True, exist_ok=True)
            dual_bound = None
            LAST_RUN_DETAILS = {
                "solver": solver_name, "n": n, "dual_bound": None,
                "objective": None, "x": None, "y": None, "time": None,
                "termination_status": "ERROR", "solve_error": repr(exc),
                "run_directory": str(run_dir),
            }
            write_json(run_dir / "result.json", LAST_RUN_DETAILS)
            print("RUN ERROR:", repr(exc))
        # Keep the same solver -> list -> n -> record structure as the input code.
        result[solver_name].append({n: dict(LAST_RUN_DETAILS)})
        write_json(output_filename, result)
        write_table(SESSION_DIR / "solver_table.csv", result)
        if CLI_ARGS.stop_on_error and (dual_bound is None or LAST_RUN_DETAILS.get("solve_error")):
            raise SystemExit(
                "Stopped before the next run because bound recovery or solve failed. "
                "The completed output and raw solver files have been saved in: "
                + LAST_RUN_DETAILS["run_directory"]
            )

    print("\nAll selected cases attempted. Results:", output_filename)
