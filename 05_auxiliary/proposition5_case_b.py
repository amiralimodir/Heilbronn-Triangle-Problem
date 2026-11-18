"""Proposition 5: nine-point complementary Case B.

n=9, LB=0, UB=0.0723767, epsilon=0.01.
MIPGap=1e-6, MIPGapAbs=0, Symmetry=2; no explicit time limit.
The x[1]>=x[0] restriction is inactive. Logs, model, parameters, status,
dual bound and coordinates are saved for each execution.
"""

from pathlib import Path
from datetime import datetime
import argparse
import hashlib
import json
import os
import platform
import shutil
import sys

SCRIPT_PATH = Path(__file__).resolve()
OUTPUT_ROOT = SCRIPT_PATH.parent / "runs"
SHOW_PLOT = False


def save_json(filename, values):
    # JSON has no Infinity/NaN literals; missing/nonfinite values stay explicit.
    def clean(value):
        if isinstance(value, float) and not math.isfinite(value):
            return None
        if isinstance(value, dict):
            return {key: clean(item) for key, item in value.items()}
        if isinstance(value, (list, tuple)):
            return [clean(item) for item in value]
        return value
    Path(filename).write_text(json.dumps(clean(values), indent=2, allow_nan=False), encoding="utf-8")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run the fixed nine-point Proposition-5 Case-B recipe; no explicit time limit.",
        allow_abbrev=False)
    parser.add_argument('--check', action='store_true',
                        help='Show the fixed Case-B inputs without importing Gurobi or optimizing.')
    args = parser.parse_args()
    if args.check:
        print(json.dumps({'n': 9, 'case': 'B', 'lower_bound': 0,
                          'upper_bound': 0.0723767, 'epsilon': 0.01,
                          'MIPGap': 1e-6, 'MIPGapAbs': 0, 'Symmetry': 2,
                          'TimeLimit': 'solver default', 'Threads': 'solver default'}, indent=2))
        raise SystemExit(0)
    import gurobipy as gp
    from gurobipy import Model,quicksum,GRB
    import matplotlib.pyplot as plt
import time
import math

# fifth group
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

            # 2. Add Relative Constraints (Gap between points)
            # "y[3] >= y[0] + 0.3" style constraints
            # Any sequence of k points (indices i to i+k-1) must span at least 'gap'
            for i in range(n - k + 1):
                j = i + k - 1
                model.addConstr(
                    y[j] >= y[i] + gap, 
                    name=f"DynamicGap_N{n}_k{k}_M{M}_{i}_{j}"
                )

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
        model.addConstr(y[1] >= 0.01 , name = 'one point on y=0')
        # model.addConstr(y[n-2] >= 0.99 , name = 'one point on y=1')
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
    # model.addConstr(x[1] >= x[0])
    model.addConstr(x[0] <= 0.5)

    # model.addConstr(x[n-2] >= x[n-1])
    # model.addConstr(x[n-1] <= 0.5)


def one_point_on_x_0_and_1(model,x,c1,c2,n):
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
            
        model.addConstr(c1.sum() <= 2, "sum_on_x0")
        model.addConstr(c2.sum() <= 2, "sum_on_x1")
    
    elif n>=4:
        for i in range(n):
            model.addConstr(x[i] <= 1 - c1[i] , name = 'One x zero')
            model.addConstr(x[i] >= c2[i] , name = 'One x 1')
            
        model.addConstr(c1.sum() >= 1, "sum_on_x0")
        model.addConstr(c1.sum() <= 2, "sum_on_x0")
        model.addConstr(c2.sum() >= 1, "sum_on_x1")
        model.addConstr(c2.sum() <= 2, "sum_on_x1")

# def x_lim(model,x,y,z,n):
#     epsilon = 1e-6
#     M = 1 / epsilon
#     for i in range(n-1):
#         model.addConstr(x[i] - x[i+1] + 2 * z <= M * (y[i+1] - y[i]))


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
            model.addConstr(w[i,n-1] == x[i])

            # First Group
            for j in range(1,n-1):
                if j != 1:
                    model.addConstr(w[i,j] == x[i]*y[j])
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
                model.addConstr(w[i,j] <= w[i,j+1])
    else:
        for i in range(n):
            # First Group
            for j in range(n):
                model.addConstr(w[i,j] == x[i]*y[j])
                if j != (n-1):
                    model.addConstr(w[i,j] <= w[i,j+1])
                

def define_phi(model,xi,phi,y):
    for i in range(n):
        for j in range(n):
            for h in range(H):
                model.addConstr( phi[i,j,h] <= xi[i,h] )
                model.addConstr( phi[i,j,h] <= y[j] )
                model.addConstr( phi[i,j,h] >= y[j] + (xi[i,h]-1) )
                model.addConstr( phi[i,j,h] >= 0 )

def define_omega(model,y,omega,ep,H):
    for i in range(n):
        for j in range(n):
            model.addConstr(omega[i,j] >= 0)
            model.addConstr(omega[i,j] >= (2**(-H))*y[j]+ep[i]-(2**(-H)))
            model.addConstr(omega[i,j] <= (2**(-H))*y[j])
            model.addConstr(omega[i,j] <= ep[i])

def sort_w(model,w):
    for i in range(n):
        model.addConstr(w[i,0] == 0)
        model.addConstr(w[i,1] == 0)
        for j in range(2,n-1):
            model.addConstr(w[i,j-1] <= w[i,j])

# fifth group
def one_point_each_squre(model,point_in_square,x,y,lb):

    grid_size = lb - 1e-6
    grid_sum = 0
    last_grid = 0
    m=0

    for i in range(50):
        if grid_sum>=1:
            grid_sum -= grid_size
            last_grid = 1 - grid_sum
            break
        grid_sum += grid_size
        m+=1


    for i in range(m):
        for j in range(m):
            model.addConstr(quicksum(point_in_square[i, j, k] for k in range(n)) <= 1, f"Square_{i}_{j}_capacity")


    for k in range(n):
        for i in range(m):
            for j in range(m):
                if j!=m-1:
                    if i!=m-1:
                        model.addConstr(point_in_square[i, j, k] * (x[k] - i * grid_size) >= 0, f"link_x_lb_{i}_{j}_{k}")
                        model.addConstr(point_in_square[i, j, k] * (x[k] - (i + 1) * grid_size) <= 0, f"link_x_ub_{i}_{j}_{k}")
                        model.addConstr(point_in_square[i, j, k] * (y[k] - j * grid_size) >= 0, f"link_y_lb_{i}_{j}_{k}")
                        model.addConstr(point_in_square[i, j, k] * (y[k] - (j + 1) * grid_size) <= 0, f"link_y_ub_{i}_{j}_{k}")
                    else:
                        model.addConstr(point_in_square[i, j, k] * (x[k] - i * grid_size) >= 0, f"link_x_lb_{i}_{j}_{k}")
                        model.addConstr(point_in_square[i, j, k] * (x[k] - (i * grid_size + last_grid)) <= 0, f"link_x_ub_{i}_{j}_{k}")
                        model.addConstr(point_in_square[i, j, k] * (y[k] - j * grid_size) >= 0, f"link_y_lb_{i}_{j}_{k}")
                        model.addConstr(point_in_square[i, j, k] * (y[k] - (j + 1) * grid_size) <= 0, f"link_y_ub_{i}_{j}_{k}")
                elif i!=m-1:
                    model.addConstr(point_in_square[i, j, k] * (x[k] - i * grid_size) >= 0, f"link_x_lb_{i}_{j}_{k}")
                    model.addConstr(point_in_square[i, j, k] * (x[k] - (i + 1) * grid_size) <= 0, f"link_x_ub_{i}_{j}_{k}")
                    model.addConstr(point_in_square[i, j, k] * (y[k] - j * grid_size) >= 0, f"link_y_lb_{i}_{j}_{k}")
                    model.addConstr(point_in_square[i, j, k] * (y[k] - (j * grid_size + last_grid)) <= 0, f"link_y_ub_{i}_{j}_{k}")
                else:
                    model.addConstr(point_in_square[i, j, k] * (y[k] - i * grid_size) >= 0, f"link_y_lb_{i}_{j}_{k}")
                    model.addConstr(point_in_square[i, j, k] * (x[k] - (i * grid_size + last_grid)) <= 0, f"link_x_ub_{i}_{j}_{k}")
                    model.addConstr(point_in_square[i, j, k] * (y[k] - j * grid_size) >= 0, f"link_y_lb_{i}_{j}_{k}")
                    model.addConstr(point_in_square[i, j, k] * (y[k] - (j * grid_size + last_grid)) <= 0, f"link_y_ub_{i}_{j}_{k}")

# fifth group
def two_point_each_rectangle(model,point_in_rectangle,y,lb):

    grid_size = (2*lb) - 1e-6
    grid_sum = 0
    last_grid = 0
    m=0

    for i in range(50):
        if grid_sum>=1:
            grid_sum -= grid_size
            last_grid = 1 - grid_sum
            break
        grid_sum += grid_size
        m+=1
    
    for i in range(m):
        model.addConstr(quicksum(point_in_rectangle[i, k] for k in range(n)) <= 2, f"rectangle_{i}_capacity")

    for k in range(n):
        for j in range(m):
            if j != m-1:
                model.addConstr(point_in_rectangle[j, k] * (y[k] - j * grid_size) >= 0, f"link_y_lb_{j}_{k}")
                model.addConstr(point_in_rectangle[j, k] * (y[k] - (j + 1) * grid_size) <= 0, f"link_y_ub_{j}_{k}")
            else:
                model.addConstr(point_in_rectangle[j, k] * (y[k] - j * grid_size) >= 0, f"link_y_lb_{j}_{k}")
                model.addConstr(point_in_rectangle[j, k] * (y[k] - (j * grid_size + last_grid)) <= 0, f"link_y_ub_{j}_{k}")                

def xi_lim(model,xi,H):
    for h in range(H):
        model.addConstr(sum(xi[i,h] for i in range(n)) <= 0.75*n)
        model.addConstr(sum(xi[i,h] for i in range(n)) >= 0.25*n)

# def enforce_convex_hull(model, x, y, n):
#     """
#     Ensures all points lie on their convex hull.
#     For any pair of points (i, j), all other points k must lie on the same side
#     of the directed line passing from i to j.
#     """
#     # Create binary variables. side[i,j,k] will capture the orientation
#     # of point k with respect to the directed line i -> j.
#     side = model.addVars(n, n, n, vtype=GRB.BINARY, name="side_var")
#     M = 2.0  # A valid 'big-M' constant for coordinates in [0, 1].

#     # Loop 1: Define the 'side' variables
#     # This loop establishes the meaning of side[i,j,k] for every possible
#     # line i->j and every other point k.
#     for i in range(n):
#         for j in range(n):
#             if i == j:
#                 continue
            
#             for k in range(n):
#                 if k == i or k == j:
#                     continue
                
#                 # This expression's sign tells us which side of line i->j point k is on.
#                 orientation_expr = (x[j] - x[i]) * (y[k] - y[i]) - (y[j] - y[i]) * (x[k] - x[i])

#                 # Link the expression's sign to the binary variable.
#                 # If orientation_expr is positive, Gurobi must set side[i,j,k] = 1.
#                 # If orientation_expr is negative, Gurobi must set side[i,j,k] = 0.
#                 model.addConstr(orientation_expr <= M * side[i, j, k], f"side_upper_{i}_{j}_{k}")
#                 model.addConstr(orientation_expr >= -M * (1 - side[i, j, k]), f"side_lower_{i}_{j}_{k}")

#     # Loop 2: Enforce the Convexity Rule
#     # This loop enforces the actual geometric property.
#     for i in range(n):
#         for j in range(n):
#             if i == j:
#                 continue

#             # For the fixed line i->j, all other points (k1, k2) must be on the same side.
#             # Therefore, their 'side' variable must be equal.
#             # We only need to loop k2 from k1+1 to avoid redundant constraints.
#             for k1 in range(n):
#                 if k1 == i or k1 == j:
#                     continue
#                 for k2 in range(k1 + 1, n):
#                     if k2 == i or k2 == j:
#                         continue
#                     model.addConstr(side[i, j, k1] == side[i, j, k2], f"consistency_{i}_{j}_{k1}_{k2}")
        

# def sum_y(model,y):
#     model.addConstr( 1 <= quicksum(y), name='lb y')
#     model.addConstr(quicksum(y) <= n-1, name='ub y')

# def sum_x(model,x):
#     model.addConstr(1 <= quicksum(x) , name = 'lb x')
#     model.addConstr(quicksum(x) <= n-1 , name= 'ub x')

# def sum_b(model,b):
#     model.addConstr( (n*(n-1)*(n-2)/(4*3*2)) <= quicksum(b), name='lb b')
#     model.addConstr( quicksum(b) <= (n*(n-1)*(n-2)/(4*2)) , name='ub b')


# def add_corner_triangle_constraints(model, n, x, y, lb):
#     lower_bound = (lb*2)**0.5

#     triangles = [
#         {"name": "bottom_left", "vertices": [(0, 0), (0, lower_bound), (lower_bound, 0)]},
#         {"name": "bottom_right", "vertices": [(1, 0), (1 - lower_bound, 0), (1, lower_bound)]},
#         {"name": "top_left", "vertices": [(0, 1), (0, 1 - lower_bound), (lower_bound, 1)]},
#         {"name": "top_right", "vertices": [(1, 1), (1 - lower_bound, 1), (1, 1 - lower_bound)]}
#     ]

#     point_in_triangle = model.addVars(len(triangles), n, vtype=GRB.BINARY, name="point_in_triangle")

#     for t_idx, triangle in enumerate(triangles):
#         vertices = triangle["vertices"]
#         x1, y1 = vertices[0]
#         x2, y2 = vertices[1]
#         x3, y3 = vertices[2]

#         model.addConstr(
#             quicksum(point_in_triangle[t_idx, k] for k in range(n)) <= 2,
#             name=f"max_points_in_{triangle['name']}_triangle"
#         )

#         for k in range(n):
#             model.addConstr(
#                 point_in_triangle[t_idx, k] * (x[k] * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2) -
#                                                y[k] * (x2 - x3) - y2 * (x3 - x1) - y3 * (x1 - x2)) >= 0,
#                 name=f"point_in_{triangle['name']}_triangle_{k}_constraint1"
#             )
#             model.addConstr(
#                 point_in_triangle[t_idx, k] * (x[k] * (y3 - y1) + x3 * (y1 - y2) + x1 * (y2 - y3) -
#                                                y[k] * (x3 - x1) - y3 * (x1 - x2) - y1 * (x2 - x3)) >= 0,
#                 name=f"point_in_{triangle['name']}_triangle_{k}_constraint2"
#             )
#             model.addConstr(
#                 point_in_triangle[t_idx, k] * (x[k] * (y1 - y2) + x1 * (y2 - y3) + x2 * (y3 - y1) -
#                                                y[k] * (x1 - x2) - y1 * (x2 - x3) - y2 * (x3 - x1)) >= 0,
#                 name=f"point_in_{triangle['name']}_triangle_{k}_constraint3"
#             )
#     return model

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


def heilbronn_triangle_approach1(n,m,m2,ub,lb,group=16):
    model = gp.Model("Heilbronn Triangle")
    x = model.addVars(n, vtype=GRB.CONTINUOUS, name="x", lb=0, ub=1)
    y = model.addVars(n, vtype=GRB.CONTINUOUS, name="y", lb=0, ub=1)
    w = model.addVars(n,n, vtype=GRB.CONTINUOUS, name="w", lb=0, ub=1)
    S = model.addVars(n, n, n, vtype=GRB.CONTINUOUS, name="S", lb=-0.5, ub=0.5)
    b = model.addVars(n, n, n, vtype=GRB.BINARY, name="b")
    z = model.addVar(vtype=GRB.CONTINUOUS, name="z", lb=0, ub=0.5)
    # point_in_square = model.addVars(m2, m2, n, vtype=GRB.BINARY, name="point_in_square")
    # point_in_rectangle_y = model.addVars(m, n, vtype=GRB.BINARY, name="point_in_rec_y")
    # point_in_rectangle_x = model.addVars(m, n, vtype=GRB.BINARY, name="point_in_rec_x")
    c1 = model.addVars(n, vtype=GRB.BINARY, name="c1")
    c2 = model.addVars(n, vtype=GRB.BINARY, name="c2")

    model.update()

    # First Group
    z.lb = lb
    z.ub = ub
    distance_points_y_0(model,x,z)

    # First & Second Group
    sort_y(model,y,n)
    define_w(model,w,x,y,n)

    # Second Group
    one_point_on_x_0_and_1(model,x,c1,c2,n)

    # Third Group
    # two_point_each_rectangle(model,point_in_rectangle_y,y,lb)
    # two_point_each_rectangle(model,point_in_rectangle_x,x,lb)
    # one_point_each_squre(model,point_in_square,x,y,lb)
    y_bounds(model,y,n)

    # x_lim(model,x,y,z,n)
    # enforce_convex_hull(model, x, y, n)

    # sum_y(model,y)
    # sum_x(model,x)
    # sum_b(model,b)

    # add_corner_triangle_constraints(model,n,x,y,lb)
    if n in [10,12]:
        model.addConstr(y[n/2]>=(1/2))
    elif n == 11:
        model.addConstr(y[(n-1)/2]>=(1/2))
    
    if n in [10,11,12]:
        for k in range(2,n):
            model.addConstr(b[0,1,k] == 1)
        for i in range(n-2):
            model.addConstr(b[i,n-2,n-1] == 1)

        split_points_x_half(model,x,n)

    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                model.addConstr(S[i, j, k] == 0.5 * (w[i,j] - w[i,k] + w[j,k] - w[j,i] + w[k,i] - w[k,j]), name=f"S_constr_{i}_{j}_{k}")
                model.addConstr((1 - b[i, j, k])*(ub+0.5) + S[i, j, k] >= z, name=f"linearize1_{i}_{j}_{k}")
                model.addConstr(b[i, j, k]*(ub+0.5) - S[i, j, k] >= z, name=f"linearize2_{i}_{j}_{k}")
                model.addConstr((-1)/2 <= S[i,j,k]-((1/2+lb)*b[i,j,k]))
                model.addConstr(S[i,j,k]-((1/2+lb)*b[i,j,k]) <= -lb)
                # model.addConstr(S[i, j, k] <= 0.5*b[i,j,k] , name="upper")
                # model.addConstr(S[i, j, k] >= 0.5*(b[i,j,k]-1) , name="lower")
    
    model.setObjective(z, GRB.MAXIMIZE)

    model.setParam('MIPGap', 1e-6)
    model.setParam('MIPGapAbs', 0.00)
    # model.setParam('FeasibilityTol', 1e-6)
    # model.setParam('IntFeasTol', 1e-6)
    model.setParam('LogFile', 'solver.log')
    # model.setParam('TimeLimit', 86000)
    model.setParam("Symmetry", 2)
    # model.setParam('ConcurrentMIP', 16)

    # x_sol = [0.18858338454351245, 1.0, 0.0, 0.6228396504398339, 1.0, 0.37716171765000717, 0.0, 0.8114187865516976]
    # y_sol = [0.0, 0.0, 0.23240576832012538, 0.23240674544937284, 0.7675914560069514, 0.7675914560069511, 1.0, 1.0]
    # for i in range(n):
    #     x[i].start = x_sol[i]
    # for i in range(n):
    #     y[i].start = y_sol[i]
    
    model.update()
    model.write("model.mps")
    model.write("parameters.prm")
    settings = {key: model.getParamInfo(key)[2] for key in (
        "TimeLimit", "MIPGap", "MIPGapAbs", "Symmetry", "NonConvex",
        "FeasibilityTol", "IntFeasTol", "OptimalityTol", "Threads", "Seed")}
    save_json("run_config.json", {
        "n": n, "input_lb": lb, "input_ub": ub, "epsilon": 0.01,
        "python_version": sys.version, "gurobi_version": list(gp.gurobi.version()),
        "platform": platform.platform(), "cpu_count": os.cpu_count(),
        "selected_effective_parameters": settings,
        "source_sha256": hashlib.sha256(SCRIPT_PATH.read_bytes()).hexdigest(),
        "recipe_basis": "supplied source plus author-confirmed LB=0 clarification",
        "time_limit_explicitly_set": False,
    })
    start_time= time.time()
    model.optimize()
    optimize_time= time.time() - start_time
    status_names = {GRB.OPTIMAL: "OPTIMAL", GRB.TIME_LIMIT: "TIME_LIMIT",
                    GRB.INFEASIBLE: "INFEASIBLE", GRB.INTERRUPTED: "INTERRUPTED"}
    save_json("result.json", {
        "status": status_names.get(model.Status, str(model.Status)),
        "status_code": model.Status, "obj_bound": model.ObjBound,
        "obj_value": model.ObjVal if model.SolCount else None,
        "mip_gap": model.MIPGap if model.SolCount else None,
        "solver_runtime": model.Runtime, "elapsed_solve_seconds": optimize_time,
        "solution_count": model.SolCount,
        "x": [x[i].X for i in range(n)] if model.SolCount else None,
        "y": [y[i].X for i in range(n)] if model.SolCount else None,
        "dual_bound_status": "floating-point bound for the implemented model",
    })
    
    
    if model.status == GRB.OPTIMAL:
        result.append(f"Optimal value: {z.X}")
        optimal_z = z.X
        optimal_x = [x[i].X for i in range(n)]
        optimal_y = [y[i].X for i in range(n)]
        optimal_b = sum(b[i,j,k].X for i in range(n) for j in range(i+1,n) for k in range(j+1,n))
        model.write("solution_optimal.sol")
        return optimal_z, optimal_x, optimal_y, optimal_b, optimize_time
    
    elif model.status == GRB.TIME_LIMIT and model.SolCount > 0:
        result.append(f"Feasible solution found under time limit (not proven optimal), value: {z.X}")
        feasible_z = z.X
        feasible_x = [x[i].X for i in range(n)]
        feasible_y = [y[i].X for i in range(n)]
        feasible_b = sum(b[i, j, k].X for i in range(n) for j in range(i + 1, n) for k in range(j + 1, n))
        model.write("solution_feasible.sol")
        return feasible_z, feasible_x, feasible_y, feasible_b, optimize_time
    
    else:
        result.append("No optimal solution found")
        return None, None, None, None, None


def heilbronn_triangle_approach2(n,m,m2,ub,lb):
    
    model = gp.Model("Heilbronn Triangle Quadratic")

    x = model.addVars(n, vtype=GRB.CONTINUOUS, name="x", lb=0, ub=1)
    y = model.addVars(n, vtype=GRB.CONTINUOUS, name="y", lb=0, ub=1)
    w = model.addVars(n,n, vtype=GRB.CONTINUOUS, name="w", lb=0, ub=1)
    S = model.addVars(n, n, n, vtype=GRB.CONTINUOUS, name="S", lb=-0.5, ub=0.5)
    U = model.addVars(n, n, n, vtype=GRB.CONTINUOUS, name="U", lb=0, ub=0.5)
    z = model.addVar(vtype=GRB.CONTINUOUS, name="z", lb=0 , ub=0.5)
    point_in_square = model.addVars(m2, m2, n, vtype=GRB.BINARY, name="point_in_square")
    point_in_rectangle_y = model.addVars(m, n, vtype=GRB.BINARY, name="point_in_rec_y")
    point_in_rectangle_x = model.addVars(m, n, vtype=GRB.BINARY, name="point_in_rec_x")
    c1 = model.addVars(n, vtype=GRB.BINARY, name="c1")
    c2 = model.addVars(n, vtype=GRB.BINARY, name="c2")


    model.update()

    #y_bounds(model,y,n)
    #sort_y(model,y)
    #distance_points_y_0(model,x,lb)
    #one_point_on_x_0_and_1(model,x,c1,c2)
    define_w(model,w,x,y,n)
    #two_point_each_rectangle(model,point_in_rectangle_y,y,lb)
    #two_point_each_rectangle(model,point_in_rectangle_x,x,lb)
    #one_point_each_squre(model,point_in_square,x,y,lb)
    # sum_y(model,y)
    # sum_x(model,x)
    # sum_b(model,b)

    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                model.addConstr(S[i,j,k] == 0.5 * (w[i,j] - w[i,k] + w[j,k] - w[j,i] + w[k,i] - w[k,j]), name=f"S_constr_{i}_{j}_{k}")
                model.addConstr(S[i,j,k]*S[i,j,k] == U[i,j,k]*U[i,j,k], name=f"quad_{i}_{j}_{k}")
                model.addConstr(U[i,j,k] >= z, name=f"U_constr_{i}_{j}_{k}")
                model.addConstr(z >= 1e-10 , name= 'Not in a line')
    
    model.setObjective(z, GRB.MAXIMIZE)
    #model.setParam('MIPGap', 1e-6)
    #model.setParam('MIPGapAbs', 0.01)
    #model.setParam('FeasibilityTol', 1e-6)
    #model.setParam('IntFeasTol', 1e-6)
    
    start_time= time.time()
    model.optimize()
    optimize_time= time.time() - start_time
    
    if model.status == GRB.OPTIMAL:
        result.append(f"Optimal value: {z.X}")
        optimal_z = z.X
        optimal_x = [x[i].X for i in range(n)]
        optimal_y = [y[i].X for i in range(n)]
        # optimal_b = sum(b[i,j,k].X for i in range(n) for j in range(i+1,n) for k in range(j+1,n))
        return optimal_z, optimal_x, optimal_y, optimize_time
    else:
        result.append("No optimal solution found")
        return None, None, None, None, None


def heilbronn_triangle_approach3_MILP(n,H,m,ub,lb):
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
    # point_in_rectangle = model.addVars(m, n, vtype=GRB.BINARY, name="point_in_rectangle")

    model.update()

    #y_bounds(model,y,n)
    #sort_y(model,y)
    #sort_w(model,w)
    define_phi(model,xi,phi,y)
    define_omega(model,y,omega,ep,H)
    
    for i in range(n):
        for j in range(n):
            model.addConstr(w[i,j] == sum(2**(-h-1) * (phi[i,j,h]) for h in range(H))+ omega[i,j])
    
    # two_point_each_rectangle(model,point_in_rectangle,y,lb)
    # sum_y(model,y)
    # sum_b(model,b)
    # xi_lim(model,xi,H)

    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                model.addConstr(S[i, j, k] == 0.5 * (w[i,j] - w[i,k] + w[j,k] - w[j,i] + w[k,i] - w[k,j]), name=f"S_constr_{i}_{j}_{k}")
                model.addConstr((1 - b[i, j, k])*(ub+0.5) + S[i, j, k] >= z, name=f"linearize1_{i}_{j}_{k}")
                model.addConstr(b[i, j, k]*(ub+0.5) - S[i, j, k] >= z, name=f"linearize2_{i}_{j}_{k}")
                model.addConstr((-ub) <= S[i,j,k]-((ub+lb)*b[i,j,k]))
                model.addConstr(S[i,j,k]-((ub+lb)*b[i,j,k]) <= -lb)
                # model.addConstr(S[i, j, k] <= 0.5*b[i,j,k] , name="upper")
                # model.addConstr(S[i, j, k] >= 0.5*(b[i,j,k]-1) , name="lower")
    
    model.setObjective(z, GRB.MAXIMIZE)
    model.setParam('TimeLimit', 86400)
    #model.setParam('MIPGap', 1e-6)
    #model.setParam('MIPGapAbs', 0.01)
    #model.setParam('FeasibilityTol', 1e-6)
    #model.setParam('IntFeasTol', 1e-6)
    
    start_time= time.time() 
    model.optimize()
    optimize_time= time.time() - start_time
    
    if model.status == GRB.OPTIMAL:
        result.append(f"Optimal value: {z.X}")
        optimal_z = z.X
        optimal_x = [ sum(2**(-h-1) * (xi[i,h].X) for h in range(H))+ (ep[i].X) for i in range(n)]
        optimal_ep = [ep[i].X for i in range(n)]
        optimal_y = [y[i].X for i in range(n)]
        optimal_b = sum(b[i,j,k].X for i in range(n) for j in range(i+1,n) for k in range(j+1,n))
        return optimal_z, optimal_x, optimal_y, optimal_b, optimal_ep, optimize_time
    else:
        result.append("No optimal solution found")
        return None, None, None, None, None, None



def heilbronn_triangle_approach3_MIQCP(n,H,m,ub,lb):
    model = gp.Model("Heilbronn Triangle")

    w = model.addVars(n,n, vtype=GRB.CONTINUOUS, name="x", lb=0, ub=1)
    # phi = model.addVars(n,n, H, vtype=GRB.CONTINUOUS, name="phi" , lb=0 , ub=1)
    xi = model.addVars(n, H, vtype=GRB.BINARY, name="xi")
    # omega = model.addVars(n,n , vtype=GRB.CONTINUOUS, name="omega", lb=0, ub=(2**(-H)))
    ep = model.addVars(n, vtype=GRB.CONTINUOUS, name="ep", lb=0, ub=(2**(-H)))
    y = model.addVars(n, vtype=GRB.CONTINUOUS, name="y", lb=0, ub=1)
    S = model.addVars(n, n, n, vtype=GRB.CONTINUOUS, name="S", lb=-0.5, ub=0.5)
    b = model.addVars(n, n, n, vtype=GRB.BINARY, name="b")
    z = model.addVar(vtype=GRB.CONTINUOUS, name="z", lb=lb, ub=ub)
    point_in_rectangle = model.addVars(m, n, vtype=GRB.BINARY, name="point_in_square")

    model.update()

    # y_bounds(model,y,n)
    # sort_y(model,y)
    # sort_w(model,w)

    for i in range(n):
        # model.addConstr(w[i,n-1] == sum(2**(-h-1) * (xi[i]) for h in range(H))+ ep[i] )
        for j in range(n):
            model.addConstr(w[i,j] == sum(2**(-h-1) * (xi[i,h]*y[j]) for h in range(H))+ (ep[i]*y[j]))
    
    # two_point_each_rectangle(model,point_in_rectangle,y,lb)
    # sum_y(model,y)
    # sum_b(model,b)
    # xi_lim(model,xi,H)

    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                model.addConstr(S[i, j, k] == 0.5 * (w[i,j] - w[i,k] + w[j,k] - w[j,i] + w[k,i] - w[k,j]), name=f"S_constr_{i}_{j}_{k}")
                model.addConstr((1 - b[i, j, k])*(ub+0.5) + S[i, j, k] >= z, name=f"linearize1_{i}_{j}_{k}")
                model.addConstr(b[i, j, k]*(ub+0.5) - S[i, j, k] >= z, name=f"linearize2_{i}_{j}_{k}")
                model.addConstr((-ub) <= S[i,j,k]-((ub+lb)*b[i,j,k]))
                model.addConstr(S[i,j,k]-((ub+lb)*b[i,j,k]) <= -lb)
                # model.addConstr(S[i, j, k] <= 0.5*b[i,j,k] , name="upper")
                # model.addConstr(S[i, j, k] >= 0.5*(b[i,j,k]-1) , name="lower")
    
    model.setObjective(z, GRB.MAXIMIZE)
    #model.setParam('MIPGap', 1e-6)
    #model.setParam('MIPGapAbs', 0.01)
    #model.setParam('FeasibilityTol', 1e-6)
    #model.setParam('IntFeasTol', 1e-6)
    
    start_time= time.time() 
    model.optimize()
    optimize_time= time.time() - start_time
    
    if model.status == GRB.OPTIMAL:
        result.append(f"Optimal value: {z.X}")
        optimal_z = z.X
        optimal_x = [ sum(2**(-h-1) * (xi[i,h].X) for h in range(H))+ (ep[i].X) for i in range(n)]
        optimal_ep = [ep[i].X for i in range(n)]
        optimal_y = [y[i].X for i in range(n)]
        optimal_b = sum(b[i,j,k].X for i in range(n) for j in range(i+1,n) for k in range(j+1,n))
        return optimal_z, optimal_x, optimal_y, optimal_b, optimal_ep, optimize_time
    else:
        result.append("No optimal solution found")
        return None, None, None, None, None, None

def plot_solution(optimal_z, optimal_x, optimal_y):
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
                    result.append(f"points optimal: {i},{j},{k}")
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
    plt.savefig('result.jpg')
    if SHOW_PLOT:
        plt.show()
    else:
        plt.close()

if __name__ == "__main__":
    result=[]
    # Author-confirmed inputs of the reported final Case-B run.
    # Zero must remain zero; the attachment's automatic log(n)/n**2 fallback
    # was not active in that run according to the author's clarification.
    n = 9
    ub = 0.0723767
    lb = 0.0
    approach = 1

    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    RUN_DIRECTORY = OUTPUT_ROOT / datetime.now().strftime("case_b_%Y%m%d_%H%M%S_%f")
    RUN_DIRECTORY.mkdir(exist_ok=False)
    if not str(RUN_DIRECTORY.resolve()).isascii():
        raise SystemExit("Set OUTPUT_ROOT to an ASCII-only writable directory.")
    shutil.copyfile(str(Path(__file__).resolve()), str(RUN_DIRECTORY / "executed_source.py"))
    os.chdir(str(RUN_DIRECTORY))

    grid_size = 2*lb - 1e-6
    M = 0
    grid_sum = 0
    for i in range(50):
        if grid_sum>=1:
            break
        grid_sum += grid_size
        M+=1

    grid_size_square = lb - 1e-6
    M2 = 0
    grid_sum2 = 0
    for i in range(50):
        if grid_sum2>=1:
            break
        grid_sum2 += grid_size_square
        M2+=1

    if approach == 1:
        optimal_z, optimal_x, optimal_y, optimal_b, optimize_time = heilbronn_triangle_approach1(n,M,M2,ub,lb)
    
        result.append(f"x = {optimal_x}")
        result.append(f"y = {optimal_y}")
        result.append(f"z = {optimal_z}")
        result.append(f"b = {optimal_b}")
        result.append(f"time = {optimize_time}")
        if optimal_x is not None and optimal_y is not None:
            plot_solution(optimal_z, optimal_x, optimal_y)

    elif approach == 2:
        optimal_z ,optimal_x, optimal_y, optimize_time = heilbronn_triangle_approach2(n,M,M2,ub,lb)
        result.append(f"x = {optimal_x}")
        result.append(f"y = {optimal_y}")
        result.append(f"z = {optimal_z}")
        result.append(f"time = {optimize_time}")
        if optimal_x is not None and optimal_y is not None:
            plot_solution(optimal_z, optimal_x, optimal_y)

    elif approach == 3:
        H = int(input('H : '))
        s = input('MIQCP or MILP ? ')
        if s == 'MILP':
            optimal_z, optimal_x, optimal_y, optimal_b, optimal_ep, optimize_time = heilbronn_triangle_approach3_MILP(n,H,M,ub,lb)
        else:
            optimal_z, optimal_x, optimal_y, optimal_b, optimal_ep, optimize_time = heilbronn_triangle_approach3_MILP(n,H,M,ub,lb)
    
        result.append(f"x = {optimal_x}")
        result.append(f"y = {optimal_y}")
        result.append(f"z = {optimal_z}")
        result.append(f"b = {optimal_b}")
        result.append(f"ep = {optimal_ep}")
        result.append(f"time = {optimize_time}")
        if optimal_x is not None and optimal_y is not None:
            plot_solution(optimal_z, optimal_x, optimal_y)

    print(result)
    Path("result.txt").write_text("\n".join(result), encoding="utf-8")
