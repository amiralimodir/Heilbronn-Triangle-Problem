import gurobipy as gp
from gurobipy import Model,quicksum,GRB
import matplotlib.pyplot as plt
import time
import math

result=[]

    

def heilbronn_triangle(n):
    model = gp.Model("Heilbronn Triangle")

    x = model.addVars(n, vtype=GRB.CONTINUOUS, name="x", lb=0, ub=1)
    y = model.addVars(n, vtype=GRB.CONTINUOUS, name="y", lb=0, ub=1)
    S = model.addVars(n, n, n, vtype=GRB.CONTINUOUS, name="S", lb=-0.5, ub=0.5)
    b = model.addVars(n, n, n, vtype=GRB.BINARY, name="b")
    z = model.addVar(vtype=GRB.CONTINUOUS, name="z", lb=0.07184691488468573, ub=n**(-1*(8/7)-(1/2000)))
    point_in_square = model.addVars(n, n, n, vtype=GRB.BINARY, name="point_in_square")
    
    model.update()
    # model.addConstr(x[2] == 0 , name = 'one point on x=0')
    # model.addConstr(x[1] == 1 , name = 'one point on y=0')
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
    
    model.setParam('MIPGap', 1e-6)
    model.setParam('MIPGapAbs', 0.01)
    model.setParam('FeasibilityTol', 1e-6)
    model.setParam('IntFeasTol', 1e-6)
    
    start_time= time.time()
    
    model.optimize()

    optimize_time= time.time() - start_time
    
    if model.status == GRB.OPTIMAL:
        result.append(f"Optimal value: {z.X}")
        optimal_z = z.X
        optimal_x = [x[i].X for i in range(n)]
        optimal_y = [y[i].X for i in range(n)]
        return optimal_z, optimal_x, optimal_y, optimize_time
    else:
        result.append("No optimal solution found")
        return None, None

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
                    result.append(f"{i},{j},{k}")
                    result.append(f"({x[i]},{y[i]})")
                    result.append(f"({x[j]},{y[j]})")
                    result.append(f"({x[k]},{y[k]})")
                    result.append(f"{area}")

                plt.plot([optimal_x[i], optimal_x[j]], [optimal_y[i], optimal_y[j]], 'b-')
                plt.plot([optimal_x[j], optimal_x[k]], [optimal_y[j], optimal_y[k]], 'b-')
                plt.plot([optimal_x[k], optimal_x[i]], [optimal_y[k], optimal_y[i]], 'b-')


    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Optimal Points and Triangles')
    plt.grid(True)
    plt.savefig('result.jpg')
    plt.show()

n = int(input())
optimal_z ,optimal_x, optimal_y, optimize_time = heilbronn_triangle(n)
result.append(f"x = {optimal_x}")
result.append(f"y = {optimal_y}")
result.append(f"time = {optimize_time}")
if optimal_x is not None and optimal_y is not None:
    plot_solution(optimal_z, optimal_x, optimal_y)



with open('result.text' , 'w') as file:
    for item in result:
        file.write(f"{item}\n")

for item in result:
    print(item)
