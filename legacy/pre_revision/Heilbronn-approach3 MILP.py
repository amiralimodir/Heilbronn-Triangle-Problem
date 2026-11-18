import gurobipy as gp
from gurobipy import Model,quicksum,GRB
import matplotlib.pyplot as plt
import time
import math

result=[]
    
def heilbronn_triangle(n,H):
    model = gp.Model("Heilbronn Triangle")

    w = model.addVars(n,n, vtype=GRB.CONTINUOUS, name="x", lb=0, ub=1)
    phi = model.addVars(n,n, H, vtype=GRB.CONTINUOUS, name="phi" , lb=0 , ub=1)
    xi = model.addVars(n, H, vtype=GRB.BINARY, name="xi")
    omega = model.addVars(n,n , vtype=GRB.CONTINUOUS, name="omega", lb=0, ub=(2**(-H)))
    ep = model.addVars(n, vtype=GRB.CONTINUOUS, name="ep", lb=0, ub=(2**(-H)))
    y = model.addVars(n, vtype=GRB.CONTINUOUS, name="y", lb=0, ub=1)
    S = model.addVars(n, n, n, vtype=GRB.CONTINUOUS, name="S", lb=-0.5, ub=0.5)
    b = model.addVars(n, n, n, vtype=GRB.BINARY, name="b")
    z = model.addVar(vtype=GRB.CONTINUOUS, name="z", lb=0, ub=1/2)

    model.update()

    for i in range(n):
        for j in range(n):
            for h in range(H):
                model.addConstr( phi[i,j,h] <= xi[i,h] )
                model.addConstr( phi[i,j,h] <= y[j] )
                model.addConstr( phi[i,j,h] >= y[j] + (xi[i,h]-1) )
                model.addConstr( phi[i,j,h] >= 0 )
            
            model.addConstr(omega[i,j] >= 0)
            model.addConstr(omega[i,j] >= (2**(-H))*y[j]+ep[i]-(2**(-H)))
            model.addConstr(omega[i,j] <= (2**(-H))*y[j])
            model.addConstr(omega[i,j] <= ep[i])

            model.addConstr(w[i,j] == sum(2**(-h) * (phi[i,j,h]) for h in range(H))+ omega[i,j])
    
    #u=n**(-1*(8/7)-(1/2000))
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                model.addConstr(S[i, j, k] == 0.5 * ((w[i,j] - w[i,k]) - (w[j,i] - w[j,k]) + (w[k,i] - w[k,j])), name=f"S_constr_{i}_{j}_{k}")
                model.addConstr((1 - b[i, j, k]) + S[i, j, k] >= z, name=f"linearize1_{i}_{j}_{k}")
                model.addConstr(b[i, j, k] - S[i, j, k] >= z, name=f"linearize2_{i}_{j}_{k}")
                model.addConstr(S[i, j, k] <= 0.5*b[i,j,k] , name="upper")
                model.addConstr(S[i, j, k] >= 0.5*(b[i,j,k]-1) , name="lower")
    
    model.setObjective(z, GRB.MAXIMIZE)
    
    start_time= time.time()
    
    model.optimize()

    optimize_time= time.time() - start_time
    
    if model.status == GRB.OPTIMAL:
        result.append(f"Optimal value: {z.X}")
        optimal_z = z.X
        return optimal_z, optimize_time
    else:
        result.append("No optimal solution found")
        return None, None

n = int(input())
optimal_z , optimize_time = heilbronn_triangle(n,5)
result.append(f"time = {optimize_time}")

with open('result.text' , 'w') as file:
    for item in result:
        file.write(f"{item}\n")

for item in result:
    print(item)