"""Run every case in the paper's two formulation-comparison tables.

Default: 15 unenhanced runs (n=3..7, Approaches 1/2/3 MILP P=10)
and 4 shared-enhancement runs (n=7,8, Approaches 1/3 MILP P=10).
The four procedural formulation functions are kept together in this file.
"""
import argparse
import csv
from datetime import datetime
import json
import math
import os
from pathlib import Path
import platform
import sys
import time


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--suite', choices=('all', 'baseline', 'shared'), default='all')
    parser.add_argument('--n', nargs='+', type=int, choices=range(3, 9))
    parser.add_argument('--approaches', nargs='+', type=int, choices=(1, 2, 3), default=[1, 2, 3])
    parser.add_argument('--p', type=int, default=10, help='Approach-3 binary expansion depth; paper: 10.')
    parser.add_argument('--approach3-mode', choices=('milp', 'miqcp'), default='milp',
                        help='Paper comparisons use MILP; MIQCP is an optional exact formulation.')
    parser.add_argument('--time-limit', type=float, default=86400)
    parser.add_argument('--output', type=Path, default=Path.home()/'heilbronn_formulations')
    parser.add_argument('--list', action='store_true', help='List selected cases without importing Gurobi or solving.')
    args = parser.parse_args()
    if args.p < 1 or not math.isfinite(args.time_limit) or args.time_limit <= 0:
        parser.error('P and the time limit must be positive.')
    return args


def select_experiments(args):
    cases = []
    for suite in ('baseline', 'shared'):
        if args.suite not in ('all', suite):
            continue
        for count in (range(3, 8) if suite == 'baseline' else (7, 8)):
            if args.n is not None and count not in args.n:
                continue
            for approach in ((1, 2, 3) if suite == 'baseline' else (1, 3)):
                if approach not in args.approaches:
                    continue
                settings = {'TimeLimit': args.time_limit}
                if suite == 'shared':
                    settings.update({'MIPGap': 1e-6, 'MIPGapAbs': 0.0})
                # Both tables use a genuine zero input bound; no automatic LB.
                signed_lb = 0.0
                objective_lb = (signed_lb if approach == 1 or
                                (approach == 3 and args.approach3_mode == 'miqcp') else 0.0)
                if approach == 2:
                    objective_lb = 1e-10  # Existing non-collinearity constraint.
                cases.append({
                    'suite': suite, 'n': count, 'approach': approach,
                    'approach3_mode': args.approach3_mode if approach == 3 else None,
                    'P': args.p if approach == 3 else None,
                    'input_lower_bound': 0.0, 'input_upper_bound': 0.5,
                    'effective_signed_area_lower_bound': signed_lb,
                    'effective_objective_lower_bound': objective_lb,
                    'solver_parameters': settings,
                    'table': 'tab:approach_combined' if suite == 'baseline' else 'tab:approach13_shared',
                    'paper_variant': args.time_limit == 86400 and (approach != 3 or
                                     (args.approach3_mode == 'milp' and args.p == 10)),
                })
    return cases


if __name__ == '__main__':
    ARGS = parse_args()
    CASES = select_experiments(ARGS)
    if not CASES:
        raise SystemExit('No cases match the selected suite, point counts, and approaches.')
    if ARGS.list:
        for case in CASES:
            print(json.dumps(case, sort_keys=True))
        print('Total runs:', len(CASES))
        raise SystemExit(0)

import gurobipy as gp
from gurobipy import GRB, quicksum


def y_bounds(model,y,n):
    if n == 6:

        y[1].ub = min(y[1].ub, 1/2)

        y[2].lb = max(y[2].lb, 1/5)
        y[2].ub = min(y[2].ub, 4/5)

        y[3].lb = max(y[3].lb, 1/5)
        y[3].ub = min(y[3].ub, 4/5)
    
        y[4].lb = max(y[4].lb, 1/2)

    if n == 7:

        y[1].ub = min(y[1].ub, 1/2)

        y[2].lb = max(y[2].lb, 1/6)
        y[2].ub = min(y[2].ub, 2/3)

        y[3].lb = max(y[3].lb, 1/6)
        y[3].ub = min(y[3].ub, 5/6)

        y[4].lb = max(y[4].lb, 1/3)
        y[4].ub = min(y[4].ub, 5/6)

        y[5].lb = max(y[5].lb, 1/2)
    
    if n == 8:

        y[1].ub = min(y[1].ub, 1/2)

        y[2].lb = max(y[2].lb, 1/7)
        y[2].ub = min(y[2].ub, 1/2)

        y[3].lb = max(y[3].lb, 1/7)
        y[3].ub = min(y[3].ub, 2/3)

        y[4].lb = max(y[4].lb, 1/3)
        y[4].ub = min(y[4].ub, 6/7)

        y[5].lb = max(y[5].lb, 1/2)
        y[5].ub = min(y[5].ub, 6/7)

        y[6].lb = max(y[6].lb, 1/2)
    
    if n == 9:
        y[2].lb = max(y[2].lb, 1/10)
        y[2].ub = min(y[2].ub, 1/2)

        y[3].lb = max(y[3].lb, 1/10)
        y[3].ub = min(y[3].ub, 2/3)

        y[4].lb = max(y[4].lb, 1/4)
        y[4].ub = min(y[4].ub, 3/4)

        y[5].lb = max(y[5].lb, 1/3)
        y[5].ub = min(y[5].ub, 9/10)

        y[6].lb = max(y[6].lb, 1/2)
        y[6].ub = min(y[6].ub, 9/10)
     
    bounds_data = {}

    if n == 10:
        # Based on LB approx 0.046537
        bounds_data[10] = [
            # (k points, M_fail)
            (3, 11),  # M=11 -> z=0.045 (<LB) -> Gap >= 1/11
            (4, 11),  # M=11 -> z=0.045 (<LB) -> Gap >= 1/11
            (5, 5),   # M=5  -> z=0.038 (<LB) -> Gap >= 1/5
            (6, 3),   # M=3  -> z=0.041 (<LB) -> Gap >= 1/3
            (7, 2),   # M=2  -> z=0.041 (<LB) -> Gap >= 1/2
            (8, 2)    # M=2  -> z=0.036 (<LB) -> Gap >= 1/2
        ]
    
    if n == 11:
        bounds_data[11] = [
            (3, 14),  # M=14 -> z=0.0357 (<LB)
            (4, 14),  # M=14 -> z=0.0357 (<LB)
            (5, 6),   # M=6  -> z=0.0320 (<LB)
            (6, 4),   # M=4  -> z=0.0312 (<LB)
            (7, 3),   # M=3  -> z=0.0279 (<LB)
            (8, 2)    # M=2  -> z=0.0361 (<LB)
        ]
    
    if n == 12:
        bounds_data[12] = [
            (3, 16),  # M=16 -> z=0.03125 (<LB)
            (4, 16),  # M=16 -> z=0.03125 (<LB)
            (5, 6),   # M=6  -> z=0.0321  (<LB)
            (6, 4),   # M=4  -> z=0.0312  (<LB)
            (7, 3),   # M=3  -> z=0.0279  (<LB)
            (8, 3),   # M=3  -> z=0.0241  (<LB)
            (9, 2)    # M=2  -> z~0.0274  (<LB) [Inferred]
        ]

    # --- Dynamic Logic for N=10 (and potentially others if added) ---
    if n in bounds_data:
        for k, M in bounds_data[n]:
            gap = 1.0 / M
            
            # 1. Update Absolute Bounds (y_min and y_max for specific points)
            # This logic mimics Proposition 8:
            # The k-th point from the bottom (index k-1) must be at least 'gap' from y=0
            if k - 1 < n:
                y[k - 1].lb = max(y[k - 1].lb, gap)
            
            # The k-th point from the top (index n-k) must be at least 'gap' from y=1
            if n - k >= 0:
                y[n - k].ub = min(y[n - k].ub, 1.0 - gap)

def sort_y(model,y,n):
    if n in [10,11,12]:
        # Second Group
        model.addConstr(y[0] == 0 , name = 'one point on y=0')
        model.addConstr(y[1] == 0 , name = 'one point on y=0')
        model.addConstr(y[n-2] == 1 , name = 'one point on y=1')
        model.addConstr(y[n-1] == 1 , name = 'one point on y=1')
        
        #1,n-2 - first group
        for i in range(1,n-2):
            model.addConstr(y[i] <= y[i+1] , name = 'Sort points')
    
    elif n == 9:
        # Second Group
        model.addConstr(y[0] == 0 , name = 'one point on y=0')
        model.addConstr(y[1] <= 0.01 , name = 'one point on y=0')
        # model.addConstr(y[1] >= 0.01 , name = 'one point on y=0')
        model.addConstr(y[n-2] >= 0.99 , name = 'one point on y=1')
        model.addConstr(y[n-1] == 1 , name = 'one point on y=1')
        
        #1,n-2 - first group
        for i in range(1,n-2):
            model.addConstr(y[i] <= y[i+1] , name = 'Sort points')
    
    elif n==7 or n == 8:
        # second group
        model.addConstr(y[0] == 0 , name = 'one point on y=0')
        model.addConstr(y[1] == 0 , name = 'one point on y=0')
        model.addConstr(y[n-1] == 1 , name = 'one point on y=1')
        
        # first group
        for i in range(1,n-1):
            model.addConstr(y[i] <= y[i+1] , name = 'Sort points')
    elif n>=4:
        # second group
        model.addConstr(y[0] == 0 , name = 'one point on y=0')
        model.addConstr(y[n-1] == 1 , name = 'one point on y=1')
        
        # first group
        for i in range(n-1):
            model.addConstr(y[i] <= y[i+1] , name = 'Sort points')
    else:
        # first group
        for i in range(n-1):
            model.addConstr(y[i] <= y[i+1] , name = 'Sort points')

def define_w(model,w,x,y,n):
    if n in [10,11,12]:
        for i in range(n):
            #Second Group
            model.addConstr(w[i,0] == 0)
            model.addConstr(w[i,1] == 0)
            model.addConstr(w[i,n-2] == x[i])
            model.addConstr(w[i,n-1] == x[i])

            #First Group
            #1,n-1
            for j in range(1,n-2):
                if j != 1:
                    model.addConstr(w[i,j] == x[i]*y[j])
                model.addConstr(w[i,j] <= w[i,j+1])

    elif n == 9:
        for i in range(n):
            #Second Group
            model.addConstr(w[i,0] == 0)
            # model.addConstr(w[i,1] == 0)
            # model.addConstr(w[i,n-2] == x[i])
            model.addConstr(w[i,n-1] == x[i])

            #First Group
            #1,n-1
            for j in range(1,n-1):
                model.addConstr(w[i,j] == x[i]*y[j])
                if j != n-1:
                    model.addConstr(w[i,j] <= w[i,j+1])
                
    elif n==8 or n==7:
        for i in range(n):
            # Second Group
            model.addConstr(w[i,0] == 0)
            model.addConstr(w[i,1] == 0)
            #model.addConstr(w[i,n-1] == x[i])

            # First Group
            for j in range(2,n):
                model.addConstr(w[i,j] == x[i]*y[j])
                model.addConstr(w[i,j-1] <= w[i,j])
    
    elif n>=4:
        for i in range(n):
            # Second Group
            model.addConstr(w[i,0] == 0)
            model.addConstr(w[i,n-1] == x[i])

            # First Group
            for j in range(0,n-1):
                if j != 0:
                    model.addConstr(w[i,j] == x[i]*y[j])
                model.addConstr(w[i,j] <= w[i,j+1])
    else:
        for i in range(n):
            # First Group
            for j in range(n):
                model.addConstr(w[i,j] == x[i]*y[j])
                if j != (n-1):
                    model.addConstr(w[i,j] <= w[i,j+1])

def define_phi(model,xi,phi,y,n,H):
    for i in range(n):
        for j in range(n):
            for h in range(H):
                model.addConstr( phi[i,j,h] <= xi[i,h] )
                model.addConstr( phi[i,j,h] <= y[j] )
                model.addConstr( phi[i,j,h] >= y[j] + (xi[i,h]-1) )
                model.addConstr( phi[i,j,h] >= 0 )

def define_omega(model,y,omega,ep,H,n):
    for i in range(n):
        for j in range(n):
            model.addConstr(omega[i,j] >= 0)
            model.addConstr(omega[i,j] >= (2**(-H))*y[j]+ep[i]-(2**(-H)))
            model.addConstr(omega[i,j] <= (2**(-H))*y[j])
            model.addConstr(omega[i,j] <= ep[i])

def sort_w(model,w,n):
    for i in range(n):
        model.addConstr(w[i,0] == 0)
        model.addConstr(w[i,1] == 0)
        for j in range(2,n):
            model.addConstr(w[i,j-1] <= w[i,j])

def heilbronn_triangle_approach1(n,m,m2,ub,lb,suite,run_dir,config):
    model = gp.Model("Heilbronn Triangle")
    x = model.addVars(n, vtype=GRB.CONTINUOUS, name="x", lb=0, ub=1)
    y = model.addVars(n, vtype=GRB.CONTINUOUS, name="y", lb=0, ub=1)
    w = model.addVars(n,n, vtype=GRB.CONTINUOUS, name="w", lb=0, ub=1)
    S = model.addVars(n, n, n, vtype=GRB.CONTINUOUS, name="S", lb=-0.5, ub=0.5)
    b = model.addVars(n, n, n, vtype=GRB.BINARY, name="b")
    z = model.addVar(vtype=GRB.CONTINUOUS, name="z", lb=0, ub=0.5)
    model.update()
    z.lb = lb
    z.ub = ub
    if suite == 'shared':
        sort_y(model,y,n)
        define_w(model,w,x,y,n)
        y_bounds(model,y,n)
    else:
        for i in range(n):
            for j in range(n):
                model.addConstr(w[i,j] == x[i]*y[j])

    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                model.addConstr(S[i, j, k] == 0.5 * (w[i,j] - w[i,k] + w[j,k] - w[j,i] + w[k,i] - w[k,j]), name=f"S_constr_{i}_{j}_{k}")
                model.addConstr((1 - b[i, j, k])*(ub+0.5) + S[i, j, k] >= z, name=f"linearize1_{i}_{j}_{k}")
                model.addConstr(b[i, j, k]*(ub+0.5) - S[i, j, k] >= z, name=f"linearize2_{i}_{j}_{k}")
                model.addConstr((-1)/2 <= S[i,j,k]-((1/2+lb)*b[i,j,k]))
                model.addConstr(S[i,j,k]-((1/2+lb)*b[i,j,k]) <= -lb)
    
    model.setObjective(z, GRB.MAXIMIZE)
    return solve_and_save(model,x,y,n,run_dir,config)


def heilbronn_triangle_approach2(n,m,m2,ub,lb,suite,run_dir,config):
    
    model = gp.Model("Heilbronn Triangle Quadratic")

    x = model.addVars(n, vtype=GRB.CONTINUOUS, name="x", lb=0, ub=1)
    y = model.addVars(n, vtype=GRB.CONTINUOUS, name="y", lb=0, ub=1)
    w = model.addVars(n,n, vtype=GRB.CONTINUOUS, name="w", lb=0, ub=1)
    S = model.addVars(n, n, n, vtype=GRB.CONTINUOUS, name="S", lb=-0.5, ub=0.5)
    U = model.addVars(n, n, n, vtype=GRB.CONTINUOUS, name="U", lb=0, ub=0.5)
    z = model.addVar(vtype=GRB.CONTINUOUS, name="z", lb=0 , ub=0.5)
    model.update()
    for i in range(n):
        for j in range(n):
            model.addConstr(w[i,j] == x[i]*y[j])

    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                model.addConstr(S[i,j,k] == 0.5 * (w[i,j] - w[i,k] + w[j,k] - w[j,i] + w[k,i] - w[k,j]), name=f"S_constr_{i}_{j}_{k}")
                model.addConstr(S[i,j,k]*S[i,j,k] == U[i,j,k]*U[i,j,k], name=f"quad_{i}_{j}_{k}")
                model.addConstr(U[i,j,k] >= z, name=f"U_constr_{i}_{j}_{k}")
                model.addConstr(z >= 1e-10 , name= 'Not in a line')
    
    model.setObjective(z, GRB.MAXIMIZE)
    return solve_and_save(model,x,y,n,run_dir,config)


def heilbronn_triangle_approach3_MILP(n,H,m,ub,lb,suite,run_dir,config):
    model = gp.Model("Heilbronn Triangle")

    w = model.addVars(n,n, vtype=GRB.CONTINUOUS, name="x", lb=0, ub=1)
    phi = model.addVars(n,n, H, vtype=GRB.CONTINUOUS, name="phi" , lb=0 , ub=1)
    xi = model.addVars(n, H, vtype=GRB.BINARY, name="xi")
    omega = model.addVars(n,n , vtype=GRB.CONTINUOUS, name="omega", lb=0, ub=(2**(-H)))
    ep = model.addVars(n, vtype=GRB.CONTINUOUS, name="ep", lb=0, ub=(2**(-H)))
    y = model.addVars(n, vtype=GRB.CONTINUOUS, name="y", lb=0, ub=1)
    S = model.addVars(n, n, n, vtype=GRB.CONTINUOUS, name="S", lb=-0.5, ub=0.5)
    b = model.addVars(n, n, n, vtype=GRB.BINARY, name="b")
    z = model.addVar(vtype=GRB.CONTINUOUS, name="z", lb=0, ub=0.5)
    model.update()

    if suite == 'shared':
        y_bounds(model,y,n)
        sort_y(model,y,n)
        sort_w(model,w,n)
    define_phi(model,xi,phi,y,n,H)
    define_omega(model,y,omega,ep,H,n)
    
    for i in range(n):
        for j in range(n):
            model.addConstr(w[i,j] == sum(2**(-h-1) * (phi[i,j,h]) for h in range(H))+ omega[i,j])
    
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                model.addConstr(S[i, j, k] == 0.5 * (w[i,j] - w[i,k] + w[j,k] - w[j,i] + w[k,i] - w[k,j]), name=f"S_constr_{i}_{j}_{k}")
                model.addConstr((1 - b[i, j, k])*(ub+0.5) + S[i, j, k] >= z, name=f"linearize1_{i}_{j}_{k}")
                model.addConstr(b[i, j, k]*(ub+0.5) - S[i, j, k] >= z, name=f"linearize2_{i}_{j}_{k}")
                model.addConstr((-ub) <= S[i,j,k]-((ub+lb)*b[i,j,k]))
                model.addConstr(S[i,j,k]-((ub+lb)*b[i,j,k]) <= -lb)
    
    model.setObjective(z, GRB.MAXIMIZE)
    return solve_and_save(model,None,y,n,run_dir,config,xi=xi,ep=ep,H=H)


def heilbronn_triangle_approach3_MIQCP(n,H,m,ub,lb,suite,run_dir,config):
    model = gp.Model("Heilbronn Triangle")

    w = model.addVars(n,n, vtype=GRB.CONTINUOUS, name="x", lb=0, ub=1)
    xi = model.addVars(n, H, vtype=GRB.BINARY, name="xi")
    ep = model.addVars(n, vtype=GRB.CONTINUOUS, name="ep", lb=0, ub=(2**(-H)))
    y = model.addVars(n, vtype=GRB.CONTINUOUS, name="y", lb=0, ub=1)
    S = model.addVars(n, n, n, vtype=GRB.CONTINUOUS, name="S", lb=-0.5, ub=0.5)
    b = model.addVars(n, n, n, vtype=GRB.BINARY, name="b")
    z = model.addVar(vtype=GRB.CONTINUOUS, name="z", lb=lb, ub=ub)

    model.update()
    if suite == 'shared':
        y_bounds(model,y,n)
        sort_y(model,y,n)
        sort_w(model,w,n)


    for i in range(n):
        for j in range(n):
            model.addConstr(w[i,j] == sum(2**(-h-1) * (xi[i,h]*y[j]) for h in range(H))+ (ep[i]*y[j]))
    

    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                model.addConstr(S[i, j, k] == 0.5 * (w[i,j] - w[i,k] + w[j,k] - w[j,i] + w[k,i] - w[k,j]), name=f"S_constr_{i}_{j}_{k}")
                model.addConstr((1 - b[i, j, k])*(ub+0.5) + S[i, j, k] >= z, name=f"linearize1_{i}_{j}_{k}")
                model.addConstr(b[i, j, k]*(ub+0.5) - S[i, j, k] >= z, name=f"linearize2_{i}_{j}_{k}")
                model.addConstr((-ub) <= S[i,j,k]-((ub+lb)*b[i,j,k]))
                model.addConstr(S[i,j,k]-((ub+lb)*b[i,j,k]) <= -lb)
    
    model.setObjective(z, GRB.MAXIMIZE)
    return solve_and_save(model,None,y,n,run_dir,config,xi=xi,ep=ep,H=H)


def numeric_attribute(model, name):
    try:
        value = float(model.getAttr(name))
        return value if math.isfinite(value) and abs(value) < GRB.INFINITY else None
    except (gp.GurobiError, AttributeError):
        return None


def solve_and_save(model,x,y,n,run_dir,config,xi=None,ep=None,H=None):
    # Use relative ASCII solver filenames inside each run directory. This also
    # avoids Gurobi 11 Windows log-path failures with non-ASCII parent folders.
    previous = Path.cwd()
    os.chdir(str(run_dir))
    try:
        model.setParam('LogFile', 'solver.log')
        for key,value in config['solver_parameters'].items():
            model.setParam(key,value)
        model.update()
        model.write('model.mps')
        model.write('parameters.prm')
        recorded = dict(config)
        recorded['python_version'] = platform.python_version()
        recorded['gurobi_version'] = list(gp.gurobi.version())
        recorded['effective_parameters'] = {
            key:model.getParamInfo(key)[2] for key in (
                'TimeLimit','MIPGap','MIPGapAbs','Threads','Seed','Symmetry',
                'NonConvex','FeasibilityTol','IntFeasTol','OptimalityTol')}
        Path('run_config.json').write_text(json.dumps(recorded,indent=2),encoding='utf-8')
        start=time.perf_counter()
        model.optimize()
        elapsed=time.perf_counter()-start
        status_names = {getattr(GRB,key):key for key in (
            'LOADED','OPTIMAL','INFEASIBLE','INF_OR_UNBD','UNBOUNDED','CUTOFF',
            'ITERATION_LIMIT','NODE_LIMIT','TIME_LIMIT','SOLUTION_LIMIT',
            'INTERRUPTED','NUMERIC','SUBOPTIMAL','INPROGRESS','USER_OBJ_LIMIT')
            if hasattr(GRB,key)}
        record = dict(config)
        record.update({
            'status_code':model.Status,'termination_status':status_names.get(model.Status,str(model.Status)),
            'dual_bound':numeric_attribute(model,'ObjBound') if model.Status not in
                         (GRB.INFEASIBLE,GRB.INF_OR_UNBD,GRB.UNBOUNDED) else None,
            'objective':numeric_attribute(model,'ObjVal') if model.SolCount else None,
            'runtime':numeric_attribute(model,'Runtime'),'solve_wall_time':elapsed,
            'mip_gap':numeric_attribute(model,'MIPGap') if model.SolCount else None,
            'solution_count':model.SolCount,'x':None,'y':None,
            'recomputed_minimum_area':None,'coordinates_in_unit_square':None,
        })
        if model.SolCount:
            if x is None:
                xs=[sum(2**(-h-1)*xi[i,h].X for h in range(H))+ep[i].X for i in range(n)]
            else:
                xs=[x[i].X for i in range(n)]
            ys=[y[i].X for i in range(n)]
            area=min(abs(.5*(xs[i]*(ys[j]-ys[k])+xs[j]*(ys[k]-ys[i])+xs[k]*(ys[i]-ys[j])))
                     for i in range(n) for j in range(i+1,n) for k in range(j+1,n))
            record.update({'x':xs,'y':ys,'recomputed_minimum_area':area,
                           'coordinates_in_unit_square':all(0<=v<=1 for v in xs+ys)})
            model.write('solution.sol')
        Path('result.json').write_text(json.dumps(record,indent=2,allow_nan=False),encoding='utf-8')
        print(json.dumps({key:record[key] for key in ('suite','n','approach','dual_bound','runtime','termination_status')}))
        return record
    finally:
        model.dispose()
        os.chdir(str(previous))


def main(args,cases):
    output=args.output.expanduser().resolve()/datetime.now().strftime('%Y%m%d_%H%M%S_%f')
    output.mkdir(parents=True,exist_ok=False)
    rows=[]
    for config in cases:
        n=config['n']
        approach=config['approach']
        H=config['P'] or args.p
        suite=config['suite']
        lb=config['effective_signed_area_lower_bound']
        ub=config['input_upper_bound']
        run_dir=output/('{}_n{}_approach{}{}'.format(suite,n,approach,'_'+args.approach3_mode if approach==3 else ''))
        run_dir.mkdir()
        # Neither comparison includes Group 4 packing constraints.
        # Keep the procedural function signatures; their grid arguments are unused.
        m=0
        try:
            if approach==1:
                row=heilbronn_triangle_approach1(n,m,0,ub,lb,suite,run_dir,config)
            elif approach==2:
                row=heilbronn_triangle_approach2(n,0,0,ub,lb,suite,run_dir,config)
            elif args.approach3_mode=='milp':
                row=heilbronn_triangle_approach3_MILP(n,H,m,ub,lb,suite,run_dir,config)
            else:
                row=heilbronn_triangle_approach3_MIQCP(n,H,m,ub,lb,suite,run_dir,config)
        except Exception as exc:
            row=dict(config,termination_status='ERROR',error=str(exc))
            (run_dir/'result.json').write_text(json.dumps(row,indent=2),encoding='utf-8')
            print('Run failed:',run_dir.name,str(exc),file=sys.stderr)
        rows.append(row)
        fields=('suite','n','approach','approach3_mode','P','dual_bound','runtime','solve_wall_time',
                'termination_status','objective','recomputed_minimum_area','coordinates_in_unit_square','error')
        with (output/'table_results.csv').open('w',newline='',encoding='utf-8') as stream:
            writer=csv.DictWriter(stream,fieldnames=fields,extrasaction='ignore')
            writer.writeheader()
            writer.writerows(rows)
    print('Results:',output)
    return 1 if any(row['termination_status']=='ERROR' for row in rows) else 0


if __name__ == '__main__':
    raise SystemExit(main(ARGS,CASES))
