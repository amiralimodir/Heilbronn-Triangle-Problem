import gurobipy as gp
from gurobipy import Model,quicksum,GRB
import matplotlib.pyplot as plt
import time
import math

def heilbronn_triangle_approach1(n,upperbound):
    model = gp.Model("Heilbronn Triangle")

    x = model.addVars(n, vtype=GRB.CONTINUOUS, name="x", lb=0, ub=1)
    y = model.addVars(n, vtype=GRB.CONTINUOUS, name="y", lb=0, ub=1)
    S = model.addVars(n, n, n, vtype=GRB.CONTINUOUS, name="S", lb=-0.5, ub=0.5)
    b = model.addVars(n, n, n, vtype=GRB.BINARY, name="b")
    z = model.addVar(vtype=GRB.CONTINUOUS, name="z", lb=math.log(n)/(n**2), ub=upperbound)
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
    model.setParam('TimeLimit', 1000)

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
                    print(i,',',j,',',k)
                    print('(',x[i],',',y[i],')')
                    print('(',x[j],',',y[j],')')
                    print('(',x[k],',',y[k],')')
                    print('Area:',area)

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
    

upper_bounds=[0.1924,0.125,0.083859,0.072376,0.054876,0.046537,0.037037,0.032599,0.026697,0.024304,0.020789,0.020528,0.020528,0.020528,0.020528]
ns=[6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
ans=[]
#for i in range(len(ns)):
#    ans.append(heilbronn_triangle_approach1(ns[i],upper_bounds[i]))

#print(ans)


ans=heilbronn_triangle_approach1(8, 0.083859)
print(ans)
if ans[0]=='optimal solution':
    plot_solution(ans[1][0], ans[1][1], ans[1][2])
