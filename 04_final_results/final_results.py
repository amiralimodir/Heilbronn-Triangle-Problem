"""Final Results: standalone Approach 1 for n=3,...,12.

Groups 1--3 are enabled; group 4 is disabled. Additional unproved structural
restrictions are active for n=10,11,12. Edit bounds_by_n.json before running.

Run: python final_results.py --n 9
Check inputs without a solver: python final_results.py --n 9 --check
"""

import argparse
import atexit
import hashlib
import json
import math
import os
import platform
import sys
import time
from datetime import datetime
from decimal import Decimal, localcontext
from fractions import Fraction
from pathlib import Path

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

    # Fixed larger-n restrictions from the previous exploratory runs.
    # These coefficients are NOT recomputed from the input LB.
    # In particular, the n=11 and n=12 entries are not claimed to have
    # certified strip-exclusion comparisons; n=12's (9,2) is inferred.
    bounds_data = {
        10: [(3,11), (4,11), (5,5), (6,3), (7,2), (8,2)],
        11: [(3,14), (4,14), (5,6), (6,4), (7,3), (8,2)],
        12: [(3,16), (4,16), (5,6), (6,4), (7,3), (8,3), (9,2)],
    }
    if n in bounds_data:
        for k, M in bounds_data[n]:
            gap = 1.0 / M
            y[k-1].lb = max(y[k-1].lb, gap)
            y[n-k].ub = min(y[n-k].ub, 1.0-gap)
            if use_relative_y_gaps:
                for i in range(n-k+1):
                    j = i+k-1
                    model.addConstr(y[j] >= y[i]+gap,
                                    name=f"DynamicGap_N{n}_k{k}_M{M}_{i}_{j}")


def sort_y(model,y,n):
    if case_b:
        model.addConstr(y[0] == 0, name="B_bottom")
        model.addConstr(y[1] >= 0.01, name="B_deficient_bottom")
        model.addConstr(y[n-1] == 1, name="B_top")
        for i in range(n-1):
            model.addConstr(y[i] <= y[i+1], name=f"B_sort_y_{i}")
        return
    if not use_group2:
        if use_y_order:
            for i in range(n-1):
                model.addConstr(y[i] <= y[i+1], name=f"sort_y_{i}")
        return
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


def distance_points_y_0(model,x,z):
    # Here y[0]=y[1]=0 whenever bottom ordering is enabled.
    if use_bottom_order and n in [10,11,12]:
        model.addConstr(x[1] >= x[0], name="bottom_x_order")
    model.addConstr(x[0] <= 0.5)


def one_point_on_x_0_and_1(model,x,c1,c2,n):
    if case_b:
        for i in range(n):
            model.addConstr(x[i] <= 1.01-c1[i], name=f"B_left_{i}")
            model.addConstr(x[i] >= c2[i]-0.01, name=f"B_right_{i}")
        model.addConstr(c1.sum() >= 1, name="B_left_min")
        model.addConstr(c1.sum() <= 2, name="B_left_max")
        model.addConstr(c2.sum() >= 1, name="B_right_min")
        model.addConstr(c2.sum() <= 2, name="B_right_max")
        return
    if n in [10,11,12]:
        for i in range(n):
            model.addConstr(x[i] <= 1 - c1[i] , name = 'One x zero')
            model.addConstr(x[i] >= c2[i] , name = 'One x 1')
            
        model.addConstr(c1.sum() == 2, "sum_on_x0")
        model.addConstr(c2.sum() == 2, "sum_on_x1")

    elif n ==9:
        for i in range(n):
            model.addConstr(x[i] <= 1.01 - c1[i] , name = 'One x zero')
            model.addConstr(x[i] >= c2[i] - 0.01 , name = 'One x 1')
            
        model.addConstr(c1.sum() == 2, "sum_on_x0")
        # model.addConstr(c1.sum() >= 1, "sum_on_x0")
        model.addConstr(c2.sum() == 2, "sum_on_x1")
        # model.addConstr(c2.sum() >= 1, "sum_on_x1")
    
    elif n>=4:
        for i in range(n):
            model.addConstr(x[i] <= 1 - c1[i] , name = 'One x zero')
            model.addConstr(x[i] >= c2[i] , name = 'One x 1')
            
        model.addConstr(c1.sum() >= 1, "sum_on_x0")
        model.addConstr(c1.sum() <= 2, "sum_on_x0")
        model.addConstr(c2.sum() >= 1, "sum_on_x1")
        model.addConstr(c2.sum() <= 2, "sum_on_x1")


def define_w(model,w,x,y,n):
    if not use_group2:
        for i in range(n):
            for j in range(n):
                model.addConstr(w[i,j] == x[i]*y[j], name=f"product_{i}_{j}")
                if use_w_order and j < n-1:
                    model.addConstr(w[i,j] <= w[i,j+1], name=f"sort_w_{i}_{j}")
        return
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
                if use_w_order:
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
                    if use_w_order:
                        model.addConstr(w[i,j] <= w[i,j+1])
                
    elif n==8 or n==7:
        for i in range(n):
            # Second Group
            model.addConstr(w[i,0] == 0)
            model.addConstr(w[i,1] == 0)
            model.addConstr(w[i,n-1] == x[i])

            # First Group
            for j in range(1,n-1):
                if j != 1:
                    model.addConstr(w[i,j] == x[i]*y[j])
                if use_w_order:
                    model.addConstr(w[i,j] <= w[i,j+1])
    
    elif n>=4:
        for i in range(n):
            # Second Group
            model.addConstr(w[i,0] == 0)
            model.addConstr(w[i,n-1] == x[i])

            # First Group
            for j in range(0,n-1):
                if j != 0:
                    model.addConstr(w[i,j] == x[i]*y[j])
                if use_w_order:
                    model.addConstr(w[i,j] <= w[i,j+1])
    else:
        for i in range(n):
            # First Group
            for j in range(n):
                model.addConstr(w[i,j] == x[i]*y[j])
                if j != (n-1):
                    if use_w_order:
                        model.addConstr(w[i,j] <= w[i,j+1])


def one_point_each_squre(model,point_in_square,x,y,lb):
    grid_size = lb - 1e-6
    if grid_size <= 0:
        raise ValueError("Square-grid lower bound must exceed 1e-6.")
    m = math.ceil(1.0 / grid_size)
    last_grid = 1.0 - (m-1)*grid_size

    for i in range(m):
        for j in range(m):
            model.addConstr(quicksum(point_in_square[i,j,k] for k in range(n)) <= 1, name=f"Square_{i}_{j}_capacity")

    for k in range(n):
        model.addConstr(quicksum(point_in_square[i,j,k] for i in range(m) for j in range(m)) == 1, name=f"square_assignment_{k}")

    for k in range(n):
        for i in range(m):
            for j in range(m):
                x_lower = i*grid_size
                x_upper = (i+1)*grid_size if i != m-1 else i*grid_size + last_grid
                y_lower = j*grid_size
                y_upper = (j+1)*grid_size if j != m-1 else j*grid_size + last_grid
                model.addConstr(x[k] >= x_lower-(1-point_in_square[i,j,k]), name=f"square_link_x_lb_{i}_{j}_{k}")
                model.addConstr(x[k] <= x_upper+(1-point_in_square[i,j,k]), name=f"square_link_x_ub_{i}_{j}_{k}")
                model.addConstr(y[k] >= y_lower-(1-point_in_square[i,j,k]), name=f"square_link_y_lb_{i}_{j}_{k}")
                model.addConstr(y[k] <= y_upper+(1-point_in_square[i,j,k]), name=f"square_link_y_ub_{i}_{j}_{k}")


def two_point_each_rectangle(model,point_in_rectangle,y,lb):
    grid_size = 2*lb - 1e-6
    if grid_size <= 0:
        raise ValueError("Rectangle-grid lower bound must exceed 5e-7.")
    m = math.ceil(1.0 / grid_size)
    last_grid = 1.0 - (m-1)*grid_size
    orientation = "x" if point_in_rectangle[0,0].VarName.startswith("point_in_rec_x") else "y"

    for i in range(m):
        model.addConstr(quicksum(point_in_rectangle[i,k] for k in range(n)) <= 2, name=f"rectangle_{orientation}_{i}_capacity")

    for k in range(n):
        model.addConstr(quicksum(point_in_rectangle[i,k] for i in range(m)) == 1, name=f"rectangle_{orientation}_assignment_{k}")

    for k in range(n):
        for j in range(m):
            lower = j*grid_size
            upper = (j+1)*grid_size if j != m-1 else j*grid_size + last_grid
            model.addConstr(y[k] >= lower-(1-point_in_rectangle[j,k]), name=f"rectangle_{orientation}_link_lb_{j}_{k}")
            model.addConstr(y[k] <= upper+(1-point_in_rectangle[j,k]), name=f"rectangle_{orientation}_link_ub_{j}_{k}")


def split_points_x_half(model, x, n):
    is_left = model.addVars(n, vtype=GRB.BINARY, name="is_left")
    
    for i in range(n):
        # Link: is_left[i] = 1  =>  x[i] <= 0.5
        model.addConstr(x[i] <= 0.5 + (1 - is_left[i]), name=f"split_x_left_{i}")
        
        # Link: is_left[i] = 0  =>  x[i] >= 0.5
        model.addConstr(x[i] >= 0.5 - is_left[i], name=f"split_x_right_{i}")

    # Enforce the specific counts based on n
    if n == 10:
        # Exact split: 5 left, 5 right
        model.addConstr(is_left.sum() == 5, name="split_balance_n10")
        
    elif n == 11:
        # Split: Either 6 left/5 right OR 5 left/6 right
        # This means the sum of is_left must be either 5 or 6
        model.addConstr(is_left.sum() >= 5, name="split_balance_n11_lb")
        model.addConstr(is_left.sum() <= 6, name="split_balance_n11_ub")
        
    elif n == 12:
        model.addConstr(is_left.sum() == 6, name="split_balance_n12")


def heilbronn_triangle_approach1(n,m,m2,ub,lb):
    global active_model
    model = gp.Model("Heilbronn Triangle")
    active_model = model
    model.setParam("OutputFlag", 1)
    model.setParam("LogToConsole", 1)
    model.setParam("LogFile", "solver.log")
    x = model.addVars(n, vtype=GRB.CONTINUOUS, name="x", lb=0, ub=1)
    y = model.addVars(n, vtype=GRB.CONTINUOUS, name="y", lb=0, ub=1)
    w = model.addVars(n,n, vtype=GRB.CONTINUOUS, name="w", lb=0, ub=1)
    S = model.addVars(n, n, n, vtype=GRB.CONTINUOUS, name="S", lb=-0.5, ub=0.5)
    b = model.addVars(n, n, n, vtype=GRB.BINARY, name="b")
    z = model.addVar(vtype=GRB.CONTINUOUS, name="z", lb=0, ub=0.5)
    if use_squares:
        point_in_square = model.addVars(m2, m2, n, vtype=GRB.BINARY, name="point_in_square")
    if use_rectangles_y:
        point_in_rectangle_y = model.addVars(m, n, vtype=GRB.BINARY, name="point_in_rec_y")
    if use_rectangles_x:
        point_in_rectangle_x = model.addVars(m, n, vtype=GRB.BINARY, name="point_in_rec_x")
    if use_x_edges:
        c1 = model.addVars(n, vtype=GRB.BINARY, name="c1")
    if use_x_edges:
        c2 = model.addVars(n, vtype=GRB.BINARY, name="c2")

    model.update()

    # First Group
    area_lb = lb if use_objective_bounds else 0.0
    z.lb = area_lb
    z.ub = ub
    if use_reflection:
        distance_points_y_0(model,x,z)

    # First & Second Group
    sort_y(model,y,n)
    define_w(model,w,x,y,n)

    # Second Group
    if use_x_edges:
        one_point_on_x_0_and_1(model,x,c1,c2,n)

    # Fourth Group: complete occupancy formulations
    if use_rectangles_y:
        two_point_each_rectangle(model,point_in_rectangle_y,y,lb)
    if use_rectangles_x:
        two_point_each_rectangle(model,point_in_rectangle_x,x,lb)
    if use_squares:
        one_point_each_squre(model,point_in_square,x,y,lb)
    # Third Group: fixed bounds and optional larger-n relative gaps
    if use_y_bounds:
        y_bounds(model,y,n)

    # Additional restrictions requested for n=10,11,12 only.
    if n in [10,11,12]:
        if use_y_median:
            model.addConstr(y[n//2] >= 0.5, name="median_y_lower")
        if use_edge_sign_fixing:
            for k in range(2,n):
                model.addConstr(b[0,1,k] == 1, name=f"bottom_sign_{k}")
            # On the top edge, positive S requires x[n-2]>=x[n-1].
            for i in range(n-2):
                model.addConstr(b[i,n-2,n-1] == 1, name=f"top_sign_{i}")
        if use_half_balance:
            split_points_x_half(model,x,n)

    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                model.addConstr(S[i, j, k] == 0.5 * (w[i,j] - w[i,k] + w[j,k] - w[j,i] + w[k,i] - w[k,j]), name=f"S_constr_{i}_{j}_{k}")
                model.addConstr((1 - b[i, j, k])*(ub+0.5) + S[i, j, k] >= z, name=f"linearize1_{i}_{j}_{k}")
                model.addConstr(b[i, j, k]*(ub+0.5) - S[i, j, k] >= z, name=f"linearize2_{i}_{j}_{k}")
                model.addConstr((-1)/2 <= S[i,j,k]-((1/2+area_lb)*b[i,j,k]))
                model.addConstr(S[i,j,k]-((1/2+area_lb)*b[i,j,k]) <= -area_lb)
                # model.addConstr(S[i, j, k] <= 0.5*b[i,j,k] , name="upper")
                # model.addConstr(S[i, j, k] >= 0.5*(b[i,j,k]-1) , name="lower")
    
    model.setObjective(z, GRB.MAXIMIZE)

    for parameter, value in solver_settings.items():
        model.setParam(parameter, value)
    model.update()
    run_data['model_size'] = {
        'variables': int(model.NumVars), 'binary_variables': int(model.NumBinVars),
        'linear_constraints': int(model.NumConstrs),
        'quadratic_constraints': int(model.NumQConstrs)
    }
    run_data['effective_y_variable_bounds'] = [
        {'index_zero_based': i, 'lower_bound': y[i].LB, 'upper_bound': y[i].UB}
        for i in range(n)
    ]
    run_data['effective_objective_bounds'] = [z.LB, z.UB]
    try:
        fingerprint = int(model.Fingerprint)
        run_data['model_fingerprint'] = fingerprint
        run_data['model_fingerprint_hex'] = '0x%08x' % (fingerprint & 0xffffffff)
    except (gp.GurobiError, AttributeError):
        run_data['model_fingerprint'] = None
        run_data['model_fingerprint_hex'] = None
    run_data['model_build_seconds'] = time.perf_counter() - run_started_timer
    run_data['actual_parameters'] = {
        key: getattr(model.Params, key) for key in solver_settings
    }
    if WRITE_MODEL:
        model.write('model.mps')
        model.write('parameters.prm')
        run_data['model_mps_sha256'] = hashlib.sha256((output_dir / 'model.mps').read_bytes()).hexdigest()
    if BUILD_ONLY:
        run_data['status'] = 'BUILT'
        run_data['status_code'] = None
        run_data['solution_count'] = 0
        run_data['ObjBound'] = None
        run_data['ObjVal'] = None
        run_data['MIPGap'] = None
        run_data['solver_runtime_seconds'] = None
        run_data['elapsed_optimize_seconds'] = None
        result.append('Model built and exported; optimization was not started.')
        return None, None, None, None, None

    start_time = time.perf_counter()
    model.optimize()
    optimize_time = time.perf_counter() - start_time
    
    
    status_names = {
        2: "OPTIMAL", 3: "INFEASIBLE", 4: "INF_OR_UNBD", 5: "UNBOUNDED",
        6: "CUTOFF", 7: "ITERATION_LIMIT", 8: "NODE_LIMIT", 9: "TIME_LIMIT",
        10: "SOLUTION_LIMIT", 11: "INTERRUPTED", 12: "NUMERIC", 13: "SUBOPTIMAL",
        14: "INPROGRESS", 15: "USER_OBJ_LIMIT", 16: "WORK_LIMIT", 17: "MEM_LIMIT"
    }
    run_data['status_code'] = int(model.Status)
    run_data['status'] = status_names.get(int(model.Status), str(model.Status))
    run_data['solution_count'] = int(model.SolCount)
    run_data['solver_runtime_seconds'] = float(model.Runtime)
    run_data['elapsed_optimize_seconds'] = optimize_time
    run_data['model_size'] = {
        'variables': int(model.NumVars), 'binary_variables': int(model.NumBinVars),
        'linear_constraints': int(model.NumConstrs), 'quadratic_constraints': int(model.NumQConstrs)
    }
    for attribute in ['ObjBound', 'NodeCount', 'Work', 'IterCount', 'BarIterCount']:
        try:
            value = float(model.getAttr(attribute))
            run_data[attribute] = value if math.isfinite(value) else str(value)
        except (gp.GurobiError, AttributeError):
            run_data[attribute] = None
    run_data['actual_parameters'] = {
        key: getattr(model.Params, key) for key in
        ['TimeLimit', 'MIPGap', 'MIPGapAbs', 'Threads', 'Seed', 'FeasibilityTol',
         'IntFeasTol', 'OptimalityTol', 'NonConvex', 'Symmetry', 'ConcurrentMIP']
    }
    run_data['ObjVal'] = None
    run_data['MIPGap'] = None
    if model.SolCount > 0:
        run_data['ObjVal'] = float(model.ObjVal)
        try:
            value = float(model.MIPGap)
            run_data['MIPGap'] = value if math.isfinite(value) else str(value)
        except (gp.GurobiError, AttributeError):
            pass
        model.write('solution.sol')
    result.append(f"Status: {run_data['status']}")
    result.append(f"ObjBound: {run_data['ObjBound']}")
    result.append(f"ObjVal: {run_data['ObjVal']}")
    result.append(f"Runtime: {run_data['solver_runtime_seconds']}")
    result.append(f"MIPGap: {run_data['MIPGap']}")
    (output_dir / 'result.json').write_text(json.dumps(run_data, indent=2, allow_nan=False), encoding='utf-8')
    if model.SolCount > 0:
        result.append(f"Incumbent objective: {z.X}")
        optimal_z = z.X
        optimal_x = [x[i].X for i in range(n)]
        optimal_y = [y[i].X for i in range(n)]
        optimal_b = sum(b[i,j,k].X for i in range(n) for j in range(i+1,n) for k in range(j+1,n))
        return optimal_z, optimal_x, optimal_y, optimal_b, optimize_time
    else:
        result.append("No incumbent retained; see status and ObjBound.")
        return None, None, None, None, optimize_time


def plot_solution(optimal_z, optimal_x, optimal_y):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.figure(figsize=(10, 10))
    plt.scatter(optimal_x, optimal_y, c='red')

    n = len(optimal_x)
    for i in range(n):
        plt.annotate(f"{i}", (optimal_x[i], optimal_y[i]), textcoords="offset points", xytext=(5, 5), ha='center')
        
    x=optimal_x
    y=optimal_y
    z=optimal_z
    minarea=100
    best_tri = None
    
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                area =abs( 0.5 * (x[i] * (y[j] - y[k]) + x[j] * (y[k] - y[i]) + x[k] * (y[i] - y[j])) )
                if area<=minarea:
                    minarea=area
                    best_tri = (i, j, k)

    for i in range(n):
        for j in range(i + 1, n): 
            for k in range(j + 1, n):
                if (i,j,k) == best_tri:
                    trianglex=[x[i],x[j],x[k]]
                    triangley=[y[i],y[j],y[k]]
                    for t in range(3):
                        plt.plot(trianglex, triangley, 'g-')
                        
                    plt.fill(trianglex, triangley)
                    result.append(f"minimum-area triangle: {i},{j},{k}")
                    result.append(f"first point: ({x[i]},{y[i]})")
                    result.append(f"second point: ({x[j]},{y[j]})")
                    result.append(f"third point: ({x[k]},{y[k]})")
                    result.append(f"area: {minarea}")

                plt.plot([optimal_x[i], optimal_x[j]], [optimal_y[i], optimal_y[j]], 'b-')
                plt.plot([optimal_x[j], optimal_x[k]], [optimal_y[j], optimal_y[k]], 'b-')
                plt.plot([optimal_x[k], optimal_x[i]], [optimal_y[k], optimal_y[i]], 'b-')

    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Computed configuration')
    plt.grid(True)
    plt.savefig(str(output_dir / 'configuration.png'))
    plt.close()


def acquire_run_lock(package_dir):
    # A process-held lock is released even if Python is killed. It prevents
    # simultaneous model execution from this package, not unrelated server jobs.
    handle = (package_dir / '.model_run.lock').open('a+b')
    handle.seek(0, os.SEEK_END)
    if handle.tell() == 0:
        handle.write(b'0')
        handle.flush()
    handle.seek(0)
    try:
        if os.name == 'nt':
            import msvcrt
            msvcrt.locking(handle.fileno(), msvcrt.LK_NBLCK, 1)
        else:
            import fcntl
            fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
    except (OSError, BlockingIOError) as error:
        handle.close()
        raise RuntimeError('Another model from this package is running. Run final experiments sequentially.') from error
    atexit.register(handle.close)
    return handle


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Final Approach 1: groups 1--3; n=3,...,12.')
    parser.add_argument('--n', type=int, choices=range(3,13), required=True,
                        help='Number of points; one instance per invocation.')
    parser.add_argument('--build-only', action='store_true',
                        help='Build/export the model without optimizing (requires Gurobi).')
    parser.add_argument('--check', action='store_true',
                        help='Validate bounds/settings without importing Gurobi or optimizing.')
    args = parser.parse_args()

    n = args.n
    EXPERIMENT_ID = 'final_n{}_approach1'.format(n)
    experiment_kind = 'final_results'
    approach = 1
    ACTIVE_GROUPS = [1, 2, 3]
    use_group1 = True
    use_group2 = True
    use_group3 = True
    use_group4 = False
    use_objective_bounds = True
    use_y_order = True
    use_w_order = True
    use_reflection = True
    use_x_edges = n >= 4
    use_rectangles_y = False
    use_rectangles_x = False
    use_squares = False
    use_y_bounds = n >= 6
    use_half_balance = n in (10,11,12)
    case_b = False
    use_bottom_order = n in (10,11,12)
    use_relative_y_gaps = n in (10,11,12)
    use_y_median = n in (10,11,12)
    use_edge_sign_fixing = n in (10,11,12)
    WRITE_MODEL = True
    SAVE_PLOT = True
    if n in (3,4,5):
        QUALIFICATION = 'G1/G2/G3 selected; G2 has no n=3 boundary fixing and G3 has no n=3..5 table entries. G4 is disabled. Solver dual bounds retain numerical qualifications.'
        Y_BOUND_MODE = 'not_applicable'
        STRUCTURAL_SCOPE = 'boundary_and_symmetry_reductions'
        UNPROVEN_RESTRICTIONS = []
    elif n in (6,7,8):
        QUALIFICATION = 'G1/G2/G3 enabled; fixed strip-derived y bounds retain their numerical qualifications. G4 is disabled. Solver objective and independently recomputed configuration area are different quantities.'
        Y_BOUND_MODE = 'fixed_paper_table_no_dynamic_callback'
        STRUCTURAL_SCOPE = 'boundary_symmetry_and_numerical_strip_bounds'
        UNPROVEN_RESTRICTIONS = []
    elif n == 9:
        QUALIFICATION = 'Nine-point case A with epsilon=0.01, fixed y bounds and x[0]<=0.5 only. No x[1]>=x[0]. Interpretation for the unrestricted optimum depends on the complementary-class exclusion and numerical-bound validity. G4 is disabled.'
        Y_BOUND_MODE = 'fixed_paper_table_no_dynamic_callback'
        STRUCTURAL_SCOPE = 'nine_point_A_near_boundary_model'
        UNPROVEN_RESTRICTIONS = ['near_boundary_A_requires_complementary_class_exclusion']
    elif True:
        QUALIFICATION = 'Restricted exploratory run with the requested larger-n restrictions: two selected points on each exact edge, fixed and relative y-coordinate restrictions, left/right half-balance and the median-y restriction. The two bottom points are ordered left to right; bottom and top triangle signs are fixed consistently with the signed-area formula. Coverage of unrestricted optima is not claimed for these combined restrictions. ObjBound applies to this implemented model. Exact evaluation of coordinates in the square supplies a feasible lower bound for the unrestricted problem. G4 is disabled.'
        Y_BOUND_MODE = 'larger_n_exploratory_fixed_and_relative_bounds'
        STRUCTURAL_SCOPE = 'restricted_exploratory_model'
        UNPROVEN_RESTRICTIONS = ['two_points_selected_on_each_exact_edge', 'larger_n_coordinate_bounds', 'relative_y_gaps_from_table', 'left_right_half_balance', 'median_y_at_least_one_half']
    FEATURE_NAMES = ['use_group1', 'use_group2', 'use_group3', 'use_group4', 'use_objective_bounds', 'use_y_order', 'use_w_order', 'use_reflection', 'use_x_edges', 'use_rectangles_y', 'use_rectangles_x', 'use_squares', 'use_y_bounds', 'use_half_balance', 'case_b', 'use_bottom_order', 'use_relative_y_gaps', 'use_y_median', 'use_edge_sign_fixing']

    BUILD_ONLY = args.build_only
    package_dir = Path(__file__).resolve().parent
    source_path = Path(__file__).resolve()
    bounds_file = package_dir / 'bounds_by_n.json'
    settings_file = package_dir / 'solver_settings.json'
    bounds_bytes = bounds_file.read_bytes()
    settings_bytes = settings_file.read_bytes()
    all_bounds = json.loads(bounds_bytes.decode('utf-8-sig'))
    solver_settings = json.loads(settings_bytes.decode('utf-8-sig'))
    expected_parameters = {
        'TimeLimit', 'MIPGap', 'MIPGapAbs', 'Threads', 'Seed',
        'FeasibilityTol', 'IntFeasTol', 'OptimalityTol', 'NonConvex',
        'Symmetry', 'ConcurrentMIP'
    }
    if set(solver_settings) != expected_parameters:
        raise ValueError('solver_settings.json must contain exactly the documented final-run parameters.')
    if any(isinstance(value, bool) or not isinstance(value, (int, float))
           or not math.isfinite(value) for value in solver_settings.values()):
        raise ValueError('All solver settings must be finite numeric values.')

    if solver_settings['TimeLimit'] != 86400:
        raise ValueError('Final runs require TimeLimit=86400 seconds per instance.')
    if any(solver_settings[key] < 0 for key in ('MIPGap', 'MIPGapAbs')):
        raise ValueError('Gap tolerances must be nonnegative.')
    for key in ('Threads', 'Seed'):
        if not isinstance(solver_settings[key], int) or solver_settings[key] < 0:
            raise ValueError(key + ' must be a nonnegative integer.')
    if not isinstance(solver_settings['NonConvex'], int) or solver_settings['NonConvex'] not in (-1, 2):
        raise ValueError('NonConvex must be -1 or 2 for these nonconvex final models.')
    if not isinstance(solver_settings['Symmetry'], int) or solver_settings['Symmetry'] not in (-1, 0, 1, 2):
        raise ValueError('Symmetry must be one of -1, 0, 1, 2.')
    if type(solver_settings['ConcurrentMIP']) is not int or solver_settings['ConcurrentMIP'] != 1:
        raise ValueError('This final-run protocol uses ConcurrentMIP=1.')
    for key in ('FeasibilityTol', 'IntFeasTol', 'OptimalityTol'):
        maximum = 0.1 if key == 'IntFeasTol' else 0.01
        if not 1e-9 <= solver_settings[key] <= maximum:
            raise ValueError(key + ' is outside the supported tolerance range.')
    entry = all_bounds[str(n)]
    configured_lower_bound = entry.get('lower_bound')
    configured_upper_bound = entry.get('upper_bound')
    requires_bound = use_group1 or use_group3 or use_group4
    if requires_bound:
        if (isinstance(configured_lower_bound, bool)
                or not isinstance(configured_lower_bound, (int, float))
                or not math.isfinite(configured_lower_bound)
                or not 0 <= configured_lower_bound <= 0.5):
            raise ValueError(f'Fill the common feasible lower_bound for n={n} in bounds_by_n.json first.')
        if (isinstance(configured_upper_bound, bool)
                or not isinstance(configured_upper_bound, (int, float))
                or not math.isfinite(configured_upper_bound)
                or not 0 < configured_upper_bound <= 0.5
                or configured_lower_bound > configured_upper_bound):
            raise ValueError(f'Invalid common upper_bound for n={n}. Use a valid upper bound for n-1.')
    geometry_lb = float(configured_lower_bound) if requires_bound else 0.0
    lb = geometry_lb
    ub = float(configured_upper_bound) if use_objective_bounds else 0.5
    objective_lb = lb if use_objective_bounds else 0.0
    if use_group4 and lb <= 1e-6:
        raise ValueError('The geometric lower bound must exceed 1e-6 for the chosen square grid.')
    M = math.ceil(1.0/(2*lb-1e-6)) if use_group4 else 0
    M2 = math.ceil(1.0/(lb-1e-6)) if use_group4 else 0
    if args.check:
        print(json.dumps({
            'n': n, 'experiment_id': EXPERIMENT_ID,
            'groups': ACTIVE_GROUPS,
            'features': {name: globals()[name] for name in FEATURE_NAMES},
            'objective_lower_bound': objective_lb,
            'objective_upper_bound': ub,
            'requested_parameters': solver_settings,
            'qualification': QUALIFICATION,
            'validation_only': True
        }, indent=2, allow_nan=False))
        raise SystemExit(0)

    import gurobipy as gp
    from gurobipy import Model, quicksum, GRB

    run_lock = acquire_run_lock(package_dir)
    output_dir = package_dir / 'results' / (EXPERIMENT_ID + '_' + datetime.now().strftime('%Y%m%d_%H%M%S_%f'))
    output_dir.mkdir(parents=True, exist_ok=False)
    result = []
    run_started_timer = time.perf_counter()
    run_data = {
        'experiment_id': EXPERIMENT_ID, 'experiment_kind': 'final_results',
        'n': n, 'approach': 1, 'groups': ACTIVE_GROUPS,
        'processing_complete': False, 'build_only': BUILD_ONLY,
        'gurobi_version': list(gp.gurobi.version()),
        'python_version': sys.version, 'platform': platform.platform(),
        'machine': platform.node(), 'logical_cpu_count': os.cpu_count(),
        'source_sha256': hashlib.sha256(source_path.read_bytes()).hexdigest(),
        'bounds_config_sha256': hashlib.sha256(bounds_bytes).hexdigest(),
        'settings_config_sha256': hashlib.sha256(settings_bytes).hexdigest(),
        'input_lower_bound': configured_lower_bound,
        'input_upper_bound': configured_upper_bound,
        'geometry_lower_bound': geometry_lb,
        'objective_lower_bound': objective_lb, 'objective_upper_bound': ub,
        'upper_bound_reference_n': entry.get('upper_bound_reference_n'),
        'upper_bound_source': entry.get('upper_bound_source'),
        'rectangle_count_per_orientation': M, 'square_count_per_axis': M2,
        'rectangle_step': 2*lb-1e-6 if use_group4 else None,
        'rectangle_last_width': 1-(M-1)*(2*lb-1e-6) if use_group4 else None,
        'square_step': lb-1e-6 if use_group4 else None,
        'square_last_width': 1-(M2-1)*(lb-1e-6) if use_group4 else None,
        'y_bound_mode': Y_BOUND_MODE,
        'dual_bound_scope': 'implemented_model_only',
        'structural_scope': STRUCTURAL_SCOPE,
        'unproven_restrictions': UNPROVEN_RESTRICTIONS,
        'features': {name: globals()[name] for name in FEATURE_NAMES},
        'qualification': QUALIFICATION,
        'requested_parameters': solver_settings,
        'run_started_local': datetime.now().astimezone().isoformat(),
        'protocol_id': 'final_G123_all_constraints_2026-11-16',
        'reference_gurobi_version': [11, 0, 3],
        'matches_reference_gurobi_version': tuple(gp.gurobi.version()) == (11, 0, 3),
        'zero_lb_policy': 'zero_is_used_unchanged',
        'solver_file_paths': 'ASCII_relative_filenames_in_run_output_directory'
    }
    (output_dir / 'run_config.json').write_text(json.dumps(run_data, indent=2, allow_nan=False), encoding='utf-8')
    (output_dir / 'executed_source.py').write_bytes(source_path.read_bytes())
    (output_dir / 'bounds_by_n.snapshot.json').write_bytes(bounds_bytes)
    (output_dir / 'solver_settings.snapshot.json').write_bytes(settings_bytes)
    print('Experiment:', EXPERIMENT_ID, '| Groups:', ACTIVE_GROUPS or 'none', flush=True)
    print('Objective bounds:', objective_lb, ub, '| Geometric lower bound:', geometry_lb, flush=True)
    print('Output directory:', output_dir, flush=True)

    active_model = None
    original_working_directory = Path.cwd()
    try:
        os.chdir(output_dir)
        optimal_z, optimal_x, optimal_y, optimal_b, optimize_time = heilbronn_triangle_approach1(n,M,M2,ub,lb)
        result.append(f'x = {optimal_x}')
        result.append(f'y = {optimal_y}')
        result.append(f'z = {optimal_z}')
        result.append(f'b = {optimal_b}')
        result.append(f'time = {optimize_time}')
        run_data['x'] = optimal_x
        run_data['y'] = optimal_y
        run_data['z'] = optimal_z
        run_data['b_sum'] = optimal_b
        if optimal_x is not None and optimal_y is not None:
            # Recompute from the exact saved decimal coordinates, without rounding/clipping.
            xs = [Fraction(str(value)) for value in optimal_x]
            ys = [Fraction(str(value)) for value in optimal_y]
            minimum_area = None
            minimum_triangle = None
            for i in range(n):
                for j in range(i+1,n):
                    for k in range(j+1,n):
                        area = abs(xs[i]*(ys[j]-ys[k])+xs[j]*(ys[k]-ys[i])+xs[k]*(ys[i]-ys[j]))/2
                        if minimum_area is None or area < minimum_area:
                            minimum_area = area
                            minimum_triangle = (i,j,k)
            coordinates_in_square = all(0 <= value <= 1 for value in xs+ys)
            run_data['coordinates_in_unit_square'] = coordinates_in_square
            run_data['minimum_area_exact_fraction'] = str(minimum_area)
            run_data['minimum_triangle_zero_based'] = list(minimum_triangle)
            run_data['minimum_triangle_one_based'] = [i+1 for i in minimum_triangle]
            with localcontext() as decimal_context:
                decimal_context.prec = max(80, len(str(minimum_area.denominator))+len(str(abs(minimum_area.numerator)))+10)
                run_data['minimum_area_exact_decimal'] = format(Decimal(minimum_area.numerator)/Decimal(minimum_area.denominator), 'f')
            run_data['configuration_is_feasible_lower_bound'] = coordinates_in_square
            result.append('Recomputed minimum area: ' + run_data['minimum_area_exact_decimal'])
            result.append('Coordinates in unit square: ' + str(coordinates_in_square))
            # Save numeric data before plotting; plotting errors cannot discard the run.
            (output_dir / 'result.json').write_text(json.dumps(run_data, indent=2, allow_nan=False), encoding='utf-8')
            if SAVE_PLOT:
                try:
                    plot_solution(optimal_z, optimal_x, optimal_y)
                except Exception as plotting_error:
                    run_data['plotting_error'] = repr(plotting_error)
        run_data['processing_complete'] = True
    except BaseException as error:
        run_data['python_error'] = repr(error)
        raise
    finally:
        if active_model is not None:
            try:
                active_model.setParam('LogFile', '')
            except Exception as log_error:
                run_data['log_close_error'] = repr(log_error)
            try:
                active_model.dispose()
            except Exception as disposal_error:
                run_data['model_disposal_error'] = repr(disposal_error)
        os.chdir(original_working_directory)
        run_data['total_elapsed_seconds'] = time.perf_counter() - run_started_timer
        run_data['run_finished_local'] = datetime.now().astimezone().isoformat()
        (output_dir / 'result.json').write_text(json.dumps(run_data, indent=2, allow_nan=False), encoding='utf-8')
        (output_dir / 'result.txt').write_text('\n'.join(result), encoding='utf-8')
        print('\n'.join(result))
        print('Saved:', output_dir, flush=True)
        run_lock.close()
