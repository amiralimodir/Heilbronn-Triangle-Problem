import gurobipy as gp
from gurobipy import Model,quicksum,GRB
import matplotlib.pyplot as plt
import time
import math



def sort_y(model,y,M):
    model.addConstr(y[0] == 0 , name = 'one point on y=0')
    model.addConstr(y[1] == 0 , name = 'one point on y=0')
    model.addConstr(y[n-1] == 1/M , name = 'one point on y=1')
    
    for i in range(1,n-1):
        model.addConstr(y[i] <= y[i+1] , name = 'Sort points')

def one_point_on_x_0_and_1(model,M,x,c1,c2):
    for i in range(n):
         model.addConstr(x[i] <= 1- c1[i] , name = 'One x zero')
         model.addConstr(x[i] >= (c2[i])/M , name = 'One x 1')
         
    model.addConstr(quicksum(c1) == 1)
    model.addConstr(quicksum(c2) == 1)

def define_w(model,M,w,x,y):
    for i in range(n):
        model.addConstr(w[i,0] == 0)
        model.addConstr(w[i,n-1] == x[i]/M)
        for j in range(1,n-1):
            model.addConstr(w[i,j] == x[i]*y[j])
            model.addConstr(w[i,j-1] <= w[i,j])

def heilbronn_triangle_approach1(n,M,ub,lb):
    model = gp.Model("Heilbronn Triangle")

    x = model.addVars(n, vtype=GRB.CONTINUOUS, name="x", lb=0, ub=1)
    y = model.addVars(n, vtype=GRB.CONTINUOUS, name="y", lb=0, ub=1/M)
    w = model.addVars(n,n, vtype=GRB.CONTINUOUS, name="w", lb=0, ub=1)
    S = model.addVars(n, n, n, vtype=GRB.CONTINUOUS, name="S", lb=-0.5, ub=0.5)
    b = model.addVars(n, n, n, vtype=GRB.BINARY, name="b")
    z = model.addVar(vtype=GRB.CONTINUOUS, name="z", lb=lb, ub=ub)
    c1 = model.addVars(n, vtype=GRB.BINARY, name="c1")
    c2 = model.addVars(n, vtype=GRB.BINARY, name="c2")

    model.update()

    #y_bounds(model,y,yb)
    sort_y(model,y,M)
    one_point_on_x_0_and_1(model,M,x,c1,c2)
    define_w(model,M,w,x,y)

    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                model.addConstr(S[i, j, k] == 0.5 * (w[i,j] - w[i,k] + w[j,k] - w[j,i] + w[k,i] - w[k,j]), name=f"S_constr_{i}_{j}_{k}")
                model.addConstr((1 - b[i, j, k])*(ub+0.5) + S[i, j, k] >= z, name=f"linearize1_{i}_{j}_{k}")
                model.addConstr(b[i, j, k]*(ub+0.5) - S[i, j, k] >= z, name=f"linearize2_{i}_{j}_{k}")
                model.addConstr((-1)/2 <= S[i,j,k]-((1/2+lb)*b[i,j,k]))
                model.addConstr(S[i,j,k]-((1/2+lb)*b[i,j,k]) <= -lb)
                model.addConstr(S[i, j, k] <= 0.5*b[i,j,k] , name="upper")
                model.addConstr(S[i, j, k] >= 0.5*(b[i,j,k]-1) , name="lower")
    
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
        optimal_b = sum(b[i,j,k].X for i in range(n) for j in range(i+1,n) for k in range(j+1,n))
        return optimal_z, optimal_x, optimal_y, optimal_b, optimize_time
    else:
        result.append("No optimal solution found")
        return None, None, None, None, None

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
                    result.append(f"points optimal: {i},{j},{k}")
                    result.append(f"first point: ({x[i]},{y[i]})")
                    result.append(f"second point: ({x[j]},{y[j]})")
                    result.append(f"third point: ({x[k]},{y[k]})")
                    result.append(f"area: {area}")

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

result=[]
b=int(input())

if b==0:
    n = int(input('n: '))
    M = int(input('M: '))
    ub = 1/(2*M)
    lb = 0


    optimal_z, optimal_x, optimal_y, optimal_b, optimize_time = heilbronn_triangle_approach1(n,M,ub,lb)

    result.append(f"x = {optimal_x}")
    result.append(f"y = {optimal_y}")
    result.append(f"z = {optimal_z}")
    result.append(f"b = {optimal_b}")
    result.append(f"time = {optimize_time}")
    if optimal_x is not None and optimal_y is not None:
        plot_solution(optimal_z, optimal_x, optimal_y)

else:
    ans=[]
    
    optimal_z=1
    for i in range(2,6):
        n=3
        while optimal_z >= 0.083859:
            optimal_z, optimal_x, optimal_y, optimal_b, optimize_time = heilbronn_triangle_approach1(n,i,1/(2*i),0)
            ans.append((i,n,optimal_z))
            n+=1
        print('M:',i,', n:',n)
    
    print(ans)
