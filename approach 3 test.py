import gurobipy as gp
from gurobipy import Model,quicksum,GRB
import matplotlib.pyplot as plt
import time
import math


def my_callback(model, where):
    if where == GRB.Callback.MIPNODE:
        rel_gap = model.cbGet(GRB.Callback.MIPNODE_OBJBST)
        bound = model.cbGet(GRB.Callback.MIPNODE_OBJBND)
        if bound != 0:
            rel_gap = abs((bound - rel_gap) / bound)
        else:
            rel_gap = float("inf")

        if rel_gap == 0.99:
            x_vals = [model.cbGetNodeRel(model.getVarByName(f"x[{i},{n-1}]")) for i in range(n)]
            y_vals = [model.cbGetNodeRel(model.getVarByName(f"y[{i}]")) for i in range(n)]
            z_val = model.cbGetNodeRel(model.getVarByName("z"))
            
            print(f"Solution with gap {rel_gap*100:.2f}%: z = {z_val}")
            print(f"x: {x_vals}")
            print(f"y: {y_vals}")
            
            result.append(f"Gap: {rel_gap*100:.2f}%")
            result.append(f"Intermediate z: {z_val}")
            result.append(f"Intermediate x: {x_vals}")
            result.append(f"Intermediate y: {y_vals}")
    elif where == GRB.Callback.MIPSOL:
        pass

def y_bounds(model,y,yb):
    for i in range(2,n-2):
        y[i].UB= yb[i-2][1]
        y[i].LB= yb[i-2][0]

def sort_y(model,y):
    model.addConstr(y[0] == 0 , name = 'one point on y=0')
    model.addConstr(y[1] == 0 , name = 'one point on y=0')
    model.addConstr(y[n-1] == 1 , name = 'one point on y=1')
    
    for i in range(1,n-1):
        model.addConstr(y[i] <= y[i+1] , name = 'Sort points')
    
def distance_points_y_0(model,x,m):
    model.addConstr(x[1]-x[0] >= 1/(m))

def one_point_on_x_0_and_1(model,x,c1,c2):
    for i in range(n):
         model.addConstr(x[i] <= 1- c1[i] , name = 'One x zero')
         model.addConstr(x[i] >= c2[i] , name = 'One x 1')
         
    model.addConstr(quicksum(c1) == 1)
    model.addConstr(quicksum(c2) == 1)

def define_w(model,w,x,y):
    for i in range(n):
        model.addConstr(w[i,0] == 0)
        model.addConstr(w[i,1] == 0)
        model.addConstr(w[i,n-1] == x[i])
        for j in range(2,n-1):
            model.addConstr(w[i,j] == x[i]*y[j])
            model.addConstr(w[i,j-1] <= w[i,j])

def define_phi(model,xi,phi,y):
    for i in range(n):
        for j in range(2,n):
            for h in range(H):
                model.addConstr( phi[i,j,h] <= xi[i,h] )
                model.addConstr( phi[i,j,h] <= y[j] )
                model.addConstr( phi[i,j,h] >= y[j] + (xi[i,h]-1) )
                model.addConstr( phi[i,j,h] >= 0 )

def define_omega(model,y,omega,ep,H):
    for i in range(n):
        for j in range(2,n):
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
    


def one_point_each_squre(model,point_in_square,x,y,m):
    for i in range(m):
        for j in range(m):
            model.addConstr(quicksum(point_in_square[i, j, k] for k in range(n)) <= 1, f"Square_{i}_{j}_capacity")

    grid_size = 1.0 / m

    for k in range(n):
        for i in range(m):
            for j in range(m):
                model.addConstr(point_in_square[i, j, k] * (x[k] - i * grid_size) >= 0, f"link_x_lb_{i}_{j}_{k}")
                model.addConstr(point_in_square[i, j, k] * (x[k] - (i + 1) * grid_size) <= 0, f"link_x_ub_{i}_{j}_{k}")
                model.addConstr(point_in_square[i, j, k] * (y[k] - j * grid_size) >= 0, f"link_y_lb_{i}_{j}_{k}")
                model.addConstr(point_in_square[i, j, k] * (y[k] - (j + 1) * grid_size) <= 0, f"link_y_ub_{i}_{j}_{k}")
    
    for i in range(m):
        model.addConstr(quicksum(point_in_square[i,j, k] for k in range(n) for j in range(m)) <= 2, f"recrangle_{i}_capacity")
    
    # for j in range(m):
    #     model.addConstr(quicksum(point_in_square[i,j, k] for k in range(n) for i in range(m)) <= 2, f"rectangle_{i}_capacity")

def one_point_each_rectangle(model,point_in_rectangle,y,m):
    for i in range(m):
        model.addConstr(quicksum(point_in_rectangle[i, k] for k in range(n)) <= 2, f"rectangle_{i}_capacity")

    grid_size = 1.0 / m

    for k in range(n):
        for j in range(m):
            model.addConstr(point_in_rectangle[j, k] * (y[k] - j * grid_size) >= 0, f"link_y_lb_{j}_{k}")
            model.addConstr(point_in_rectangle[j, k] * (y[k] - (j + 1) * grid_size) <= 0, f"link_y_ub_{j}_{k}")

def xi_lim(model,xi,H):
    for h in range(H):
        model.addConstr(sum(xi[i,h] for i in range(n)) <= 0.75*n)
        model.addConstr(sum(xi[i,h] for i in range(n)) >= 0.25*n)

# def sum_y(model,y):
#     model.addConstr( 1 <= quicksum(y), name='lb y')
#     model.addConstr(quicksum(y) <= n-1, name='ub y')

# def sum_x(model,x):
#     model.addConstr(1 <= quicksum(x) , name = 'lb x')
#     model.addConstr(quicksum(x) <= n-1 , name= 'ub x')

def sum_b(model,b):
    model.addConstr( (n*(n-1)*(n-2)/(4*3*2)) <= quicksum(b), name='lb b')
    model.addConstr( quicksum(b) <= (n*(n-1)*(n-2)/(4*2)) , name='ub b')



def heilbronn_triangle_approach3_MILP(n,H,m,ub,lb,yb):
    model = gp.Model("Heilbronn Triangle")

    w = model.addVars(n,n, vtype=GRB.CONTINUOUS, name="x", lb=0, ub=1)
    phi = model.addVars(n,n, H, vtype=GRB.CONTINUOUS, name="phi" , lb=0 , ub=1)
    xi = model.addVars(n, H, vtype=GRB.BINARY, name="xi")
    omega = model.addVars(n,n , vtype=GRB.CONTINUOUS, name="omega", lb=0, ub=(2**(-H)))
    ep = model.addVars(n, vtype=GRB.CONTINUOUS, name="ep", lb=0, ub=(2**(-H)))
    y = model.addVars(n, vtype=GRB.CONTINUOUS, name="y", lb=0, ub=1)
    S = model.addVars(n, n, n, vtype=GRB.CONTINUOUS, name="S", lb=-0.5, ub=0.5)
    b = model.addVars(n, n, n, vtype=GRB.BINARY, name="b")
    z = model.addVar(vtype=GRB.CONTINUOUS, name="z", lb=lb, ub=ub)
    point_in_rectangle = model.addVars(m, n, vtype=GRB.BINARY, name="point_in_square")

    model.update()

    y_bounds(model,y,yb)
    sort_y(model,y)
    sort_w(model,w)
    define_phi(model,xi,phi,y)
    define_omega(model,y,omega,ep,H)
    
    for i in range(n):
        for j in range(2,n):
            model.addConstr(w[i,j] == sum(2**(-h-1) * (phi[i,j,h]) for h in range(H))+ omega[i,j])
    
    #one_point_each_rectangle(model,point_in_rectangle,y,m)
    # sum_y(model,y)
    # sum_b(model,b)
    xi_lim(model,xi,H)

    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                model.addConstr(S[i, j, k] == 0.5 * (w[i,j] - w[i,k] + w[j,k] - w[j,i] + w[k,i] - w[k,j]), name=f"S_constr_{i}_{j}_{k}")
                model.addConstr((1 - b[i, j, k])*(ub+0.5) + S[i, j, k] >= z, name=f"linearize1_{i}_{j}_{k}")
                model.addConstr(b[i, j, k]*(ub+0.5) - S[i, j, k] >= z, name=f"linearize2_{i}_{j}_{k}")
                model.addConstr(S[i, j, k] <= 0.5*b[i,j,k] , name="upper")
                model.addConstr(S[i, j, k] >= 0.5*(b[i,j,k]-1) , name="lower")
    
    model.setObjective(z, GRB.MAXIMIZE)
    #model.setParam('MIPGap', 1e-6)
    #model.setParam('MIPGapAbs', 0.01)
    #model.setParam('FeasibilityTol', 1e-6)
    #model.setParam('IntFeasTol', 1e-6)

    model.setParam('Method', 1)
    model.setParam('MIPFocus', 1)
    model.setParam('Cuts', 2)
    model.setParam('Heuristics', 0.5)
    
    start_time= time.time() 
    model.update()
    print(f"Quadratic constraints: {model.NumQConstrs}")
    model.printStats()
    model.optimize(my_callback)
    optimize_time= time.time() - start_time
    
    if model.status == GRB.OPTIMAL:
        result.append(f"Optimal value: {z.X}")
        optimal_z = z.X
        optimal_x = [ w[i,n-1].X/y[n-1].X for i in range(n)]
        optimal_ep = [ep[i].X for i in range(n)]
        optimal_y = [y[i].X for i in range(n)]
        optimal_b = sum(b[i,j,k].X for i in range(n) for j in range(i+1,n) for k in range(j+1,n))
        return optimal_z, optimal_x, optimal_y, optimal_b, optimal_ep, optimize_time
    else:
        result.append("No optimal solution found")
        return None, None, None, None, None, None



def heilbronn_triangle_approach3_MIQCP(n,H,m,ub,lb,yb):
    model = gp.Model("Heilbronn Triangle")

    w = model.addVars(n,n, vtype=GRB.CONTINUOUS, name="x", lb=0, ub=1)
    #phi = model.addVars(n,n, H, vtype=GRB.CONTINUOUS, name="phi" , lb=0 , ub=1)
    xi = model.addVars(n, H, vtype=GRB.BINARY, name="xi")
    #omega = model.addVars(n,n , vtype=GRB.CONTINUOUS, name="omega", lb=0, ub=(2**(-H)))
    ep = model.addVars(n, vtype=GRB.CONTINUOUS, name="ep", lb=0, ub=(2**(-H)))
    y = model.addVars(n, vtype=GRB.CONTINUOUS, name="y", lb=0, ub=1)
    S = model.addVars(n, n, n, vtype=GRB.CONTINUOUS, name="S", lb=-0.5, ub=0.5)
    b = model.addVars(n, n, n, vtype=GRB.BINARY, name="b")
    z = model.addVar(vtype=GRB.CONTINUOUS, name="z", lb=lb, ub=ub)
    point_in_rectangle = model.addVars(m, n, vtype=GRB.BINARY, name="point_in_square")

    model.update()

    y_bounds(model,y,yb)
    sort_y(model,y)
    sort_w(model,w)

    for i in range(n):
        model.addConstr(w[i,n-1] == sum(2**(-h-1) * (xi[i]) for h in range(H))+ ep[i] )
        for j in range(2,n-1):
            model.addConstr(w[i,j] == sum(2**(-h-1) * (xi[i,h]*y[j]) for h in range(H))+ (ep[i]*y[j]))
    
    one_point_each_rectangle(model,point_in_rectangle,y,m)
    # sum_y(model,y)
    # sum_b(model,b)
    xi_lim(model,xi,H)

    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                model.addConstr(S[i, j, k] == 0.5 * (w[i,j] - w[i,k] + w[j,k] - w[j,i] + w[k,i] - w[k,j]), name=f"S_constr_{i}_{j}_{k}")
                model.addConstr((1 - b[i, j, k])*(ub+0.5) + S[i, j, k] >= z, name=f"linearize1_{i}_{j}_{k}")
                model.addConstr(b[i, j, k]*(ub+0.5) - S[i, j, k] >= z, name=f"linearize2_{i}_{j}_{k}")
                model.addConstr(S[i, j, k] <= 0.5*b[i,j,k] , name="upper")
                model.addConstr(S[i, j, k] >= 0.5*(b[i,j,k]-1) , name="lower")
    
    model.setObjective(z, GRB.MAXIMIZE)
    #model.setParam('MIPGap', 1e-6)
    #model.setParam('MIPGapAbs', 0.01)
    #model.setParam('FeasibilityTol', 1e-6)
    #model.setParam('IntFeasTol', 1e-6)

    model.setParam('NonConvex', 2)
    model.setParam('Method', 2)
    model.setParam('MIQCPMethod', 0)
    
    start_time= time.time() 
    model.update()
    print(f"Quadratic constraints: {model.NumQConstrs}")
    model.printStats()
    model.optimize()
    optimize_time= time.time() - start_time
    
    if model.status == GRB.OPTIMAL:
        result.append(f"Optimal value: {z.X}")
        optimal_z = z.X
        optimal_x = [ w[i,n-1].X/y[n-1].X for i in range(n)]
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
n = int(input('n: '))
ub = float(input('Upper bound ? '))
lb = float(input('Lower bound ? '))
if lb == 0:
    lb=math.log(n)/(n**2)

ML=[4,6,7,10,11]
if n>=6 and n<=10:
    M = ML[n-6]
else:
    M=int(input('m: '))
    
YB= [[(0.1666,0.6667),(0.1666,0.8334),(0.3333,0.8334),(0.3333,1)],[(0.142857,0.714286),(0.142857,0.714286),(0.285714,0.857143),(0.285714,0.857143),(0.428571,1)]]
if n>=7 and n<=8:
    yb=YB[n-7]
else:
    yb=[(0, 1) for _ in range(n-3)]

H = int(input('H : '))
s = input('MIQCP or MILP ? ')
if s == 'MILP':
    optimal_z, optimal_x, optimal_y, optimal_b, optimal_ep, optimize_time = heilbronn_triangle_approach3_MILP(n,H,M,ub,lb,yb)
else:
    optimal_z, optimal_x, optimal_y, optimal_b, optimal_ep, optimize_time = heilbronn_triangle_approach3_MILP(n,H,M,ub,lb,yb)

result.append(f"x = {optimal_x}")
result.append(f"y = {optimal_y}")
result.append(f"z = {optimal_z}")
result.append(f"b = {optimal_b}")
result.append(f"ep = {optimal_ep}")
result.append(f"time = {optimize_time}")
if optimal_x is not None and optimal_y is not None:
    plot_solution(optimal_z, optimal_x, optimal_y)


