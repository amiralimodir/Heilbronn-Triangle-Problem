import gurobipy as gp
from gurobipy import Model,quicksum,GRB
import matplotlib.pyplot as plt
import time
import math

result=[]

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
    model.addConstr(y[n-1] == 1 , name = 'one point on y=1')
    for i in range(n-1):
        model.addConstr(y[i] <= y[i+1] , name = 'Sort points')
    
    # equal_ij = model.addVars(n, n, vtype=GRB.BINARY, name="equal_ij")
    # equal_ik = model.addVars(n, n, vtype=GRB.BINARY, name="equal_ik")
    # equal_jk = model.addVars(n, n, vtype=GRB.BINARY, name="equal_jk")
    # d_ki = model.addVars(n, n, vtype=GRB.BINARY, name="d_ki")

    # M = 1e6  # A large constant for the big-M method
    # epsilon = 1e-6  # A small tolerance value

    # for i in range(n):
    #     for j in range(i + 1, n):
    #         model.addConstr(x[i] - x[j] <= epsilon + M * (1 - equal_ij[i, j]), name=f"x_eq_{i}_{j}_ub")
    #         model.addConstr(x[j] - x[i] <= epsilon + M * (1 - equal_ij[i, j]), name=f"x_eq_{i}_{j}_lb")

    #     for k in range(i + 1, n):
    #         model.addConstr(x[i] - x[k] <= epsilon + M * (1 - equal_ik[i, k]), name=f"x_eq_{i}_{k}_ub")
    #         model.addConstr(x[k] - x[i] <= epsilon + M * (1 - equal_ik[i, k]), name=f"x_eq_{i}_{k}_lb")

    #         model.addConstr(x[j] - x[k] <= epsilon + M * (1 - equal_jk[j, k]), name=f"x_eq_{j}_{k}_ub")
    #         model.addConstr(x[k] - x[j] <= epsilon + M * (1 - equal_jk[j, k]), name=f"x_eq_{j}_{k}_lb")
    
    # for k in range(n):
    #     for i in range(k + 1, n):
    #         model.addConstr(x[k] - x[i] >= epsilon - M * (1 - d_ki[k, i]), name=f"d_ki_lb_{k}_{i}")
    #         model.addConstr(x[i] - x[k] >= epsilon - M * d_ki[k, i], name=f"d_ki_ub_{k}_{i}")


    # for i in range(n):
    #     for j in range(i + 1, n):
    #         for k in range(j + 1, n):
    #             # x[i] = x[j] and x[k] > x[i] -> b[i,j,k] = 0
    #             model.addConstr(b[i, j, k] <= 1 - equal_ij[i, j] + d_ki[k, i], name=f"b_0_{i}_{j}_{k}")

    #             # x[i] = x[j] and x[k] < x[i] -> b[i,j,k] = 1
    #             model.addConstr(b[i, j, k] >= equal_ij[i, j] - d_ki[k, i], name=f"b_1_{i}_{j}_{k}")

    #             # x[i] = x[k] and x[j] > x[i] -> b[i,j,k] = 0
    #             model.addConstr(b[i, j, k] <= 1 - equal_ik[i, k] + d_ki[j, i], name=f"b_2_{i}_{j}_{k}")

    #             # x[i] = x[k] and x[j] < x[i] -> b[i,j,k] = 1
    #             model.addConstr(b[i, j, k] >= equal_ik[i, k] - d_ki[j, i], name=f"b_3_{i}_{j}_{k}")

    #             # x[j] = x[k] and x[i] > x[j] -> b[i,j,k] = 1
    #             model.addConstr(b[i, j, k] >= equal_jk[j, k] - d_ki[i, j], name=f"b_4_{i}_{j}_{k}")

    #             # x[j] = x[k] and x[i] < x[j] -> b[i,j,k] = 0
    #             model.addConstr(b[i, j, k] <= 1 - equal_jk[j, k] + d_ki[i, j], name=f"b_5_{i}_{j}_{k}")

    u=n**(-1*(8/7)-(1/2000))
    
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                model.addConstr(S[i, j, k] == 0.5 * (x[i] * (y[j] - y[k]) + x[j] * (y[k] - y[i]) + x[k] * (y[i] - y[j])), name=f"S_constr_{i}_{j}_{k}")
                model.addConstr((1 - b[i, j, k])*(u+0.5) + S[i, j, k] >= z, name=f"linearize1_{i}_{j}_{k}")
                model.addConstr(b[i, j, k]*(u+0.5) - S[i, j, k] >= z, name=f"linearize2_{i}_{j}_{k}")
                model.addConstr(S[i, 
                
                j, k] <= 0.5*b[i,j,k] , name="upper")
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
    
    model.addConstr( (n*(n-1)*(n-2)/(4*3))*0.9 <= sum(b[i,j,k] for i in range(n) for j in range(i+1,n) for k in range(j+1,n)), name='lb b')
    model.addConstr(sum(b[i,j,k] for i in range(n) for j in range(i+1,n) for k in range(j+1,n)) <= (n*(n-1)*(n-2)/(4*3))*1.1 , name='ub b')

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
    
    start_time= time.time()
    
    model.optimize()

    optimize_time= time.time() - start_time
    
    if model.status == GRB.OPTIMAL:
        result.append(f"Optimal value: {z.X}")
        optimal_z = z.X
        optimal_x = [x[i].X for i in range(n)]
        optimal_y = [y[i].X for i in range(n)]
        optimal_b = sum(b[i,j,k].X for i in range(n) for j in range(i+1,n) for k in range(j+1,n))
        return optimal_z, optimal_b, optimal_x, optimal_y, optimize_time
    else:
        result.append("No optimal solution found")
        return None, None, None, None
    
def heilbronn_triangle_approach2(n):
    
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
        result.append(f"Optimal value: {z.X}")
        optimal_z = z.X
        optimal_x = [x[i].X for i in range(n)]
        optimal_y = [y[i].X for i in range(n)]
        return optimal_z,optimal_x, optimal_y, optimize_time
    else:
        result.append("No optimal solution found")
        return None, None, None, None

def heilbronn_triangle_approach3_MILP(n,H,ub):
    model = gp.Model("Heilbronn Triangle")

    w = model.addVars(n,n, vtype=GRB.CONTINUOUS, name="x", lb=0, ub=1)
    phi = model.addVars(n,n, H, vtype=GRB.CONTINUOUS, name="phi" , lb=0 , ub=1)
    xi = model.addVars(n, H, vtype=GRB.BINARY, name="xi")
    omega = model.addVars(n,n , vtype=GRB.CONTINUOUS, name="omega", lb=0, ub=(2**(-H)))
    ep = model.addVars(n, vtype=GRB.CONTINUOUS, name="ep", lb=0, ub=(2**(-H)))
    y = model.addVars(n, vtype=GRB.CONTINUOUS, name="y", lb=0, ub=1)
    S = model.addVars(n, n, n, vtype=GRB.CONTINUOUS, name="S", lb=-0.5, ub=0.5)
    b = model.addVars(n, n, n, vtype=GRB.BINARY, name="b")
    z = model.addVar(vtype=GRB.CONTINUOUS, name="z", lb=math.log(n)/(n**2), ub=ub)

    model.update()
    
    model.addConstr(y[0] == 0 , name = 'one point on y=0')
    for i in range(n-1):
        model.addConstr(y[i] <= y[i+1] , name = 'Sort points')
        
    model.addConstr(y[n-1] == 1 , name = 'one point on y=1')
    
    for i in range(n):
        model.addConstr(w[i,0] == 0 , name = 'one point on w=0')
        for j in range(n-1):
            model.addConstr(w[i,j] <= w[i,j+1] , name = 'Sort points')

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

            model.addConstr(w[i,j] == sum(2**(-h-1) * (phi[i,j,h]) for h in range(H))+ omega[i,j])
    
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                model.addConstr(S[i, j, k] == 0.5 * ((w[i,j] - w[i,k]) - (w[j,i] - w[j,k]) + (w[k,i] - w[k,j])), name=f"S_constr_{i}_{j}_{k}")
                model.addConstr((1 - b[i, j, k]) + S[i, j, k] >= z, name=f"linearize1_{i}_{j}_{k}")
                model.addConstr(b[i, j, k] - S[i, j, k] >= z, name=f"linearize2_{i}_{j}_{k}")
                model.addConstr(S[i, j, k] <= 0.5*b[i,j,k] , name="upper")
                model.addConstr(S[i, j, k] >= 0.5*(b[i,j,k]-1) , name="lower")
    
    model.addConstr( 1 <= quicksum(y), name='lb y')
    model.addConstr(quicksum(y) <= n-1, name='ub y')
    
    model.addConstr( (n*(n-1)*(n-2)/(4*3))-1 <= quicksum(b), name='lb b')
    model.addConstr(quicksum(b) <= (n*(n-1)*(n-2)/(4*3))+1 , name='ub b')
    
    for h in range(H):
        model.addConstr(sum(xi[i,h] for i in range(n)) <= 0.75*n)
        model.addConstr(sum(xi[i,h] for i in range(n)) >= 0.25*n)
    
    model.setObjective(z, GRB.MAXIMIZE)
    
    start_time= time.time()
    
    model.optimize()

    optimize_time= time.time() - start_time
    
    if model.status == GRB.OPTIMAL:
        result.append(f"Optimal value: {z.X}")
        optimal_z = z.X
        optimal_x = [ w[i,n-1].X/y[n-1].X for i in range(n)]
        optimal_ep = [ep[i].X for i in range(n)]
        optimal_y = [y[i].X for i in range(n)]
        return optimal_z, optimal_x,optimal_ep,optimal_y,optimize_time
    else:
        result.append("No optimal solution found")
        return None, None


def heilbronn_triangle_approach3_MIQCP(n,H,ub):
    model = gp.Model("Heilbronn Triangle")

    w = model.addVars(n,n, vtype=GRB.CONTINUOUS, name="x", lb=0, ub=1)
    phi = model.addVars(n,n, H, vtype=GRB.CONTINUOUS, name="phi" , lb=0 , ub=1)
    xi = model.addVars(n, H, vtype=GRB.BINARY, name="xi")
    omega = model.addVars(n,n , vtype=GRB.CONTINUOUS, name="omega", lb=0, ub=(2**(-H)))
    ep = model.addVars(n, vtype=GRB.CONTINUOUS, name="ep", lb=0, ub=(2**(-H)))
    y = model.addVars(n, vtype=GRB.CONTINUOUS, name="y", lb=0, ub=1)
    S = model.addVars(n, n, n, vtype=GRB.CONTINUOUS, name="S", lb=-0.5, ub=0.5)
    b = model.addVars(n, n, n, vtype=GRB.BINARY, name="b")
    z = model.addVar(vtype=GRB.CONTINUOUS, name="z", lb=math.log(n)/(n**2), ub=ub)

    model.update()
    
    model.addConstr(y[0] == 0 , name = 'one point on y=0')
    for i in range(n-1):
        model.addConstr(y[i] <= y[i+1] , name = 'Sort points')
        
    model.addConstr(y[n-1] == 1 , name = 'one point on y=1')
    
    for i in range(n):
        model.addConstr(w[i,0] == 0 , name = 'one point on w=0')
        for j in range(n-1):
            model.addConstr(w[i,j] <= w[i,j+1] , name = 'Sort points')

    for i in range(n):
        for j in range(n):
            for h in range(H):
                model.addConstr( phi[i,j,h] == xi[i,h] * y[j]) 

            model.addConstr(omega[i,j] == ep[i] * y[j])

            model.addConstr(w[i,j] == sum(2**(-h-1) * (phi[i,j,h]) for h in range(H))+ omega[i,j])
    
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                model.addConstr(S[i, j, k] == 0.5 * ((w[i,j] - w[i,k]) - (w[j,i] - w[j,k]) + (w[k,i] - w[k,j])), name=f"S_constr_{i}_{j}_{k}")
                model.addConstr((1 - b[i, j, k]) + S[i, j, k] >= z, name=f"linearize1_{i}_{j}_{k}")
                model.addConstr(b[i, j, k] - S[i, j, k] >= z, name=f"linearize2_{i}_{j}_{k}")
                model.addConstr(S[i, j, k] <= 0.5*b[i,j,k] , name="upper")
                model.addConstr(S[i, j, k] >= 0.5*(b[i,j,k]-1) , name="lower")
    
    model.addConstr( 1 <= quicksum(y), name='lb y')
    model.addConstr(quicksum(y) <= n-1, name='ub y')
    
    model.addConstr( (n*(n-1)*(n-2)/(4*3))-1 <= quicksum(b), name='lb b')
    model.addConstr(quicksum(b) <= (n*(n-1)*(n-2)/(4*3))+1 , name='ub b')
    
    for h in range(H):
        model.addConstr(sum(xi[i,h] for i in range(n)) <= 0.75*n)
        model.addConstr(sum(xi[i,h] for i in range(n)) >= 0.25*n)
    
    model.setObjective(z, GRB.MAXIMIZE)
    
    start_time= time.time()
    
    model.optimize()

    optimize_time= time.time() - start_time
    
    if model.status == GRB.OPTIMAL:
        result.append(f"Optimal value: {z.X}")
        optimal_z = z.X
        optimal_x = [ w[i,n-1].X/y[n-1].X for i in range(n)]
        optimal_ep = [ep[i].X for i in range(n)]
        optimal_y = [y[i].X for i in range(n)]
        return optimal_z, optimal_x,optimal_ep,optimal_y,optimize_time
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

n = int(input('n: '))
ub = float(input('Upper bound ? '))
m = int(input('Approch? '))

if m == 1:
    optimal_z, optimal_b,optimal_x, optimal_y, optimize_time = heilbronn_triangle_approach1(n,ub)
    print(optimal_b)
    result.append(f"x = {optimal_x}")
    result.append(f"y = {optimal_y}")
    result.append(f"time = {optimize_time}")
    if optimal_x is not None and optimal_y is not None:
        plot_solution(optimal_z, optimal_x, optimal_y)

if m == 2:
    optimal_z ,optimal_x, optimal_y, optimize_time = heilbronn_triangle_approach2(n)
    result.append(f"x = {optimal_x}")
    result.append(f"y = {optimal_y}")
    result.append(f"time = {optimize_time}")
    if optimal_x is not None and optimal_y is not None:
        plot_solution(optimal_z, optimal_x, optimal_y)

elif m == 3:
    H = int(input('H: '))
    s = input('MIQCP or MILP ? ')
    if s == 'MILP':
        optimal_z, optimal_x,optimal_ep,optimal_y,optimize_time = heilbronn_triangle_approach3_MILP(n,H,ub)
    else:
        optimal_z, optimal_x,optimal_ep,optimal_y,optimize_time = heilbronn_triangle_approach3_MILP(n,H,ub)
    #print(optimal_ep, 'epsilon')
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
