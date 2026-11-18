import gurobipy as gp
from gurobipy import Model,quicksum,GRB
import matplotlib.pyplot as plt
import time
import math
import pandas as pd

def heilbronn_triangle_approach1(n,ub):
    model = gp.Model("Heilbronn Triangle")

    x = model.addVars(n, vtype=GRB.CONTINUOUS, name="x", lb=0, ub=1)
    y = model.addVars(n, vtype=GRB.CONTINUOUS, name="y", lb=0, ub=1)
    S = model.addVars(n, n, n, vtype=GRB.CONTINUOUS, name="S", lb=-0.5, ub=0.5)
    b = model.addVars(n, n, n, vtype=GRB.BINARY, name="b")
    z = model.addVar(vtype=GRB.CONTINUOUS, name="z", lb=math.log(n)/(n**2), ub=ub)
    point_in_square = model.addVars(n, n, n, vtype=GRB.BINARY, name="point_in_square")
    
    model.update()
    #model.addConstr(x[2] == 0 , name = 'one point on x=0')
    #model.addConstr(x[1] == 1 , name = 'one point on y=0')
    model.addConstr(y[0] == 0 , name = 'one point on y=0')
    for i in range(n-1):
        model.addConstr(y[i] <= y[i+1] , name = 'Sort points')
    model.addConstr(y[n-1] == 1 , name = 'one point on y=1')
    
    u=n**(-1*(8/7)-(1/2000))

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


            

    model.addConstr(1 <=quicksum(x) , name = 'lb x')
    model.addConstr(quicksum(x) <= n-1 , name= 'ub x')
    model.addConstr( 1 <= quicksum(y), name='lb y')
    model.addConstr(quicksum(y) <= n-1, name='ub y')
    
    model.addConstr( n*(n-1)*(n-2)/(4*3*2) <= quicksum(b), name='lb b')
    model.addConstr(quicksum(b) <= n*(n-1)*(n-2)/(4*2) , name='ub b')



    # if n%2 == 0:
    #     model.addConstr(n/4 <=quicksum(x) , name = 'lb x')
    #     model.addConstr(quicksum(x) <= 3*n/4 , name= 'ub x')
    #     model.addConstr( n/4 <= quicksum(y), name='lb y')
    #     model.addConstr(quicksum(y) <= 3*n/4, name='ub y')
    # else:
    #     model.addConstr((n-1)/4 <=quicksum(x) , name = 'lb x')
    #     model.addConstr(quicksum(x) <= (3*n+1)/4 , name= 'ub x')
    #     model.addConstr( (n-1)/4 <= quicksum(y), name='lb y')
    #     model.addConstr(quicksum(y) <= (3*n+1)/4, name='ub y')

    model.setObjective(z, GRB.MAXIMIZE)
    
    #model.setParam('MIPGap', 1e-6)
    #model.setParam('MIPGapAbs', 0.01)
    #model.setParam('FeasibilityTol', 1e-6)
    #model.setParam('IntFeasTol', 1e-6)
    #model.setParam(GRB.Param.Threads,96)
    model.setParam('TimeLimit', 4000)

    model.optimize()
    
    if model.status == GRB.OPTIMAL:
        optimal_x = [x[i].X for i in range(n)]
        optimal_y = [y[i].X for i in range(n)]
        return ('optimal solution',(z.X,optimal_x,optimal_y))
    elif model.status == GRB.TIME_LIMIT:
        best_bound = model.ObjBound
        if model.SolCount > 0:
            incumbent_solution = model.objVal
            return (incumbent_solution,best_bound)
        else:
            return('No incumbent',best_bound)


data={'N':[6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]}
ns=[6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
ans=heilbronn_triangle_approach1(6,0.1924)
data['Incumbent']=['-']
data['BestBd']=['-']
data['Optimal Z']=[ans[1][0]]
data['Optimal X']=[ans[1][1]]
data['Optimal Y']=[ans[1][2]]

for i in range(1,len(ns)):
    if ans[0] == 'optimal solution':
        ans=heilbronn_triangle_approach1(ns[i],ans[1][0])
    elif ans[0] == 'No incumbent':
        ans=heilbronn_triangle_approach1(ns[i],ans[1])
    else:
        ans=heilbronn_triangle_approach1(ns[i],ans[0])
    
    if ans[0] == 'optimal solution':
        data['Incumbent'].append('-')
        data['BestBd'].append('-')
        data['Optimal Z'].append(ans[1][0])
        data['Optimal X'].append(ans[1][1])
        data['Optimal Y'].append(ans[1][2])
    else:
        data['Incumbent'].append(ans[0])
        data['BestBd'].append(ans[1])
        data['Optimal Z'].append('_')
        data['Optimal X'].append('_')
        data['Optimal Y'].append('_')

print(data)
df = pd.DataFrame(data)
df.to_excel("result 6 - 20 points.xlsx", index=False)