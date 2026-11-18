import gurobipy as gp
from gurobipy import Model, quicksum, GRB
import pandas as pd
import math
from multiprocessing import Pool, cpu_count

def heilbronn_triangle_approach1(n, ub):
    model = gp.Model("Heilbronn Triangle")

    # Define variables
    x = model.addVars(n, vtype=GRB.CONTINUOUS, name="x", lb=0, ub=1)
    y = model.addVars(n, vtype=GRB.CONTINUOUS, name="y", lb=0, ub=1)
    S = model.addVars(n, n, n, vtype=GRB.CONTINUOUS, name="S", lb=-0.5, ub=0.5)
    b = model.addVars(n, n, n, vtype=GRB.BINARY, name="b")
    z = model.addVar(vtype=GRB.CONTINUOUS, name="z", lb=math.log(n)/(n**2), ub=ub)
    point_in_square = model.addVars(n, n, n, vtype=GRB.BINARY, name="point_in_square")
    
    # Constraints
    model.addConstr(y[0] == 0 , name = 'one point on y=0')
    model.addConstr(y[n-1] == 1 , name = 'one point on y=1')
    for i in range(n-1):
        model.addConstr(y[i] <= y[i+1] , name = 'Sort points')

    u = n**(-1*(8/7)-(1/2000))
    
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                model.addConstr(S[i, j, k] == 0.5 * (x[i] * (y[j] - y[k]) + x[j] * (y[k] - y[i]) + x[k] * (y[i] - y[j])), name=f"S_constr_{i}_{j}_{k}")
                model.addConstr((1 - b[i, j, k])*(u+0.5) + S[i, j, k] >= z, name=f"linearize1_{i}_{j}_{k}")
                model.addConstr(b[i, j, k]*(u+0.5) - S[i, j, k] >= z, name=f"linearize2_{i}_{j}_{k}")
                model.addConstr(S[i, j, k] <= 0.5*b[i,j,k] , name="upper")
                model.addConstr(S[i, j, k] >= 0.5*(b[i,j,k]-1) , name="lower")
    
    for i in range(n):
        for j in range(n):
            model.addConstr(quicksum(point_in_square[i, j, k] for k in range(n)) <= 1, f"Square_{i}_{j}_capacity")

    grid_size = 1.0 / n

    for k in range(n):
        for i in range(n):
            for j in range(n):
                model.addConstr(point_in_square[i, j, k] * (x[k] - i * grid_size) >= 0, f"link_x_lb_{i}_{j}_{k}")
                model.addConstr(point_in_square[i, j, k] * (x[k] - (i + 1) * grid_size) <= 0, f"link_x_ub_{i}_{j}_{k}")
                model.addConstr(point_in_square[i, j, k] * (y[k] - j * grid_size) >= 0, f"link_y_lb_{i}_{j}_{k}")
                model.addConstr(point_in_square[i, j, k] * (y[k] - (j + 1) * grid_size) <= 0, f"link_y_ub_{i}_{j}_{k}")

    model.addConstr(1 <= quicksum(x), name='lb x')
    model.addConstr(quicksum(x) <= n-1, name='ub x')
    model.addConstr(1 <= quicksum(y), name='lb y')
    model.addConstr(quicksum(y) <= n-1, name='ub y')
    
    model.addConstr( (n*(n-1)*(n-2)/(4*3*2)) <= quicksum(b), name='lb b')
    model.addConstr( quicksum(b) <= (n*(n-1)*(n-2)/(4*2)) , name='ub b')

    model.setObjective(z, GRB.MAXIMIZE)
    
    model.setParam('TimeLimit', 4000)
    model.optimize()

    if model.status == GRB.OPTIMAL:
        optimal_x = [x[i].X for i in range(n)]
        optimal_y = [y[i].X for i in range(n)]
        return ('optimal solution', (z.X, optimal_x, optimal_y))
    elif model.status == GRB.TIME_LIMIT:
        best_bound = model.ObjBound
        if model.SolCount > 0:
            incumbent_solution = model.objVal
            return (incumbent_solution, best_bound)
        else:
            return ('No incumbent', best_bound)

def process_instance(n, ub):
    return heilbronn_triangle_approach1(n, ub)

def main():
    ns = [(6,0.1924),(7,0.125),(8,0.083859),(9,0.072376),(10,0.054876),(11,0.046537),(12,0.037037),(13,0.032599),(14,0.026697),(15,0.024304),(16,0.020789),(17,0.020528),(18,0.020528),(19,0.020528),(20,0.020528)]
    N = list(range(6,21))
    data = {'N': N, 'Incumbent': [], 'BestBd': [], 'Optimal Z': [], 'Optimal X': [], 'Optimal Y': []}

    with Pool(cpu_count()) as pool:
        results = pool.starmap(process_instance, ns)

    for result in results:
        if result[0] == 'optimal solution':
            data['Incumbent'].append('-')
            data['BestBd'].append('-')
            data['Optimal Z'].append(result[1][0])
            data['Optimal X'].append(result[1][1])
            data['Optimal Y'].append(result[1][2])
        else:
            data['Incumbent'].append(result[0])
            data['BestBd'].append(result[1])
            data['Optimal Z'].append('_')
            data['Optimal X'].append('_')
            data['Optimal Y'].append('_')

    df = pd.DataFrame(data)
    df.to_excel("result_6_20_points.xlsx", index=False)

if __name__ == "__main__":
    main()
