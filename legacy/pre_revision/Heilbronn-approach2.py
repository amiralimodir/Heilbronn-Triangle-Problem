import gurobipy as gp
from gurobipy import GRB,quicksum
import matplotlib.pyplot as plt
import time

def heilbronn_triangle(n):
    model = gp.Model("Heilbronn Triangle Quadratic")

    x = model.addVars(n, vtype=GRB.CONTINUOUS, name="x", lb=0, ub=1)
    y = model.addVars(n, vtype=GRB.CONTINUOUS, name="y", lb=0, ub=1)
    S = model.addVars(n, n, n, vtype=GRB.CONTINUOUS, name="S", lb=-0.5, ub=0.5)
    U = model.addVars(n, n, n, vtype=GRB.CONTINUOUS, name="U", lb=0, ub=0.5)
    z = model.addVar(vtype=GRB.CONTINUOUS, name="z", lb=0, ub=0.5)

    model.update()
    model.addConstr(x[0] == 0 , name = 'one point on x=0')
    model.addConstr(x[1] == 1 , name = 'one point on y=0')
    model.addConstr(y[2] == 0 , name = 'one point on x=1')
    model.addConstr(y[3] == 1 , name = 'one point on y=1')

    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                model.addConstr(S[i,j,k] == 0.5 * (x[i]*(y[j] - y[k]) + x[j]*(y[k] - y[i]) + x[k]*(y[i] - y[j])), name=f"S_constr_{i}_{j}_{k}")
                model.addConstr(S[i,j,k]*S[i,j,k] == U[i,j,k]*U[i,j,k], name=f"quad_{i}_{j}_{k}")
                model.addConstr(U[i,j,k] >= z, name=f"U_constr_{i}_{j}_{k}")
                model.addConstr(z >= 1e-10 , name= 'Not in a line')

    model.addConstr(1 <=quicksum(x) , name = 'lb x')
    model.addConstr(quicksum(x) <= n-1 , name= 'ub x')
    model.addConstr( 1 <= quicksum(y), name='lb y')
    model.addConstr(quicksum(y) <= n-1, name='ub y')

    model.setObjective(z, GRB.MAXIMIZE)

    start_time= time.time()
    model.optimize()
    optimize_time= time.time() - start_time

    if model.status == GRB.OPTIMAL:
        print(f"Optimal value: {z.X}")
        optimal_x = [x[i].X for i in range(n)]
        optimal_y = [y[i].X for i in range(n)]
        return optimal_x, optimal_y, optimize_time
    else:
        print("No optimal solution found")
        return None, None

def plot_solution(optimal_x, optimal_y):
    plt.figure(figsize=(10, 10))
    plt.scatter(optimal_x, optimal_y, c='red')

    n = len(optimal_x)
    for i in range(n):
        plt.annotate(f"{i}", (optimal_x[i], optimal_y[i]), textcoords="offset points", xytext=(5, 5), ha='center')
    x=optimal_x
    y=optimal_y
    minarea=100
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                area =abs( 0.5 * (x[i] * (y[j] - y[k]) + x[j] * (y[k] - y[i]) + x[k] * (y[i] - y[j])) )
                if area<=minarea:
                    minarea=area
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                area =abs( 0.5 * (x[i] * (y[j] - y[k]) + x[j] * (y[k] - y[i]) + x[k] * (y[i] - y[j])) )
                if(area == minarea):
                    trianglex=[x[i],x[j],x[k]]
                    triangley=[y[i],y[j],y[k]]
                    for t in range(3):
                        plt.plot(trianglex, triangley, 'g-')
                    plt.fill(trianglex, triangley)
                    print(i,j,k)
                    print((x[i],y[i]))
                    print((x[j],y[j]))
                    print((x[k],y[k]))
                    print(area)

                plt.plot([optimal_x[i], optimal_x[j]], [optimal_y[i], optimal_y[j]], 'b-')
                plt.plot([optimal_x[j], optimal_x[k]], [optimal_y[j], optimal_y[k]], 'b-')
                plt.plot([optimal_x[k], optimal_x[i]], [optimal_y[k], optimal_y[i]], 'b-')


    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Optimal Points and Triangles')
    plt.grid(True)
    plt.show()

n = int(input())
optimal_x, optimal_y, optimize_time = heilbronn_triangle(n)
print('x = ',optimal_x)
print('y = ',optimal_y)
print('time = ',optimize_time)
if optimal_x is not None and optimal_y is not None:
    plot_solution(optimal_x, optimal_y)