from pyomo.environ import *
import matplotlib.pyplot as plt
import math
import time

def heilbronn_triangle_approach3_MILP(n, H, m, ub, lb, yb):
    model = ConcreteModel()

    # Sets
    model.i = RangeSet(1, n)  # Points
    model.j = RangeSet(1, n)  # Secondary points
    model.k = RangeSet(1, n)  # Tertiary points
    model.h = RangeSet(1, H)  # Binary precision
    model.m = RangeSet(1, m)  # Grid squares

    # Variables
    model.w = Var(model.i, model.j, bounds=(0, 1), domain=Reals)
    model.phi = Var(model.i, model.j, model.h, bounds=(0, 1), domain=Reals)
    model.xi = Var(model.i, model.h, domain=Binary)
    model.omega = Var(model.i, model.j, bounds=(0, 2**(-H)), domain=Reals)
    model.ep = Var(model.i, bounds=(0, 2**(-H)), domain=Reals)
    model.y = Var(model.j, bounds=(0, 1), domain=Reals)
    model.S = Var(model.i, model.j, model.k, bounds=(-0.5, 0.5), domain=Reals)
    model.b = Var(model.i, model.j, model.k, domain=Binary)
    model.z = Var(bounds=(lb, ub), domain=Reals)

    # Constraints
    # Bounds on y
    def y_bounds_lower_rule(model, j):
        if j > 2 and j < n - 2:
            return model.y[j] >= yb[j - 2][0]
        return Constraint.Skip

    def y_bounds_upper_rule(model, j):
        if j >= 2 and j < n - 2:
            return model.y[j] <= yb[j - 2][1]
        return Constraint.Skip

    model.y_bounds_lower = Constraint(model.j, rule=y_bounds_lower_rule)
    model.y_bounds_upper = Constraint(model.j, rule=y_bounds_upper_rule)

    # Sorting constraints for y
    def sort_y_rule(model, j):
        if j < n:
            return model.y[j] <= model.y[j + 1]
        return Constraint.Skip
    
    def y_0_rule(model, j):
        if j==0:
            return model.y[j] == 0
        return Constraint.Skip
    
    def y_1_rule(model, j):
        if j==1:
            return model.y[j] == 0
        return Constraint.Skip

    def y_n_rule(model, j):
        if j==n-1:
            return model.y[j] == 1
        return Constraint.Skip
    
    model.sort_y = Constraint(model.j, rule=sort_y_rule)
    model.y_0 = Constraint(model.j, rule=y_0_rule)
    model.y_1 = Constraint(model.j, rule=y_1_rule)
    model.y_n = Constraint(model.j, rule=y_n_rule)

    # Sorting constraints for w
    def sort_w_rule(model, i, j):
        if j > 2 and j < n - 1:
            return model.w[i, j - 1] <= model.w[i, j]
        return Constraint.Skip

    def w_i_0_rule(model, i, j):
        if j==0:
            return model.w[i, j] == 0
        return Constraint.Skip
    
    def w_i_1_rule(model, i, j):
        if j==1:
            return model.w[i, j] == 0
        return Constraint.Skip

    model.sort_w = Constraint(model.i, model.j, rule=sort_w_rule)
    model.w_i_0 = Constraint(model.i, model.j, rule=w_i_0_rule)
    model.w_i_1 = Constraint(model.i, model.j, rule=w_i_1_rule)


    # Define phi constraints
    def define_phi_rule_upper_xi(model, i, j, h):
        return model.phi[i, j, h] <= model.xi[i, h]

    def define_phi_rule_upper_y(model, i, j, h):
        return model.phi[i, j, h] <= model.y[j]

    def define_phi_rule_lower_y_xi(model, i, j, h):
        return model.phi[i, j, h] >= model.y[j] + (model.xi[i, h] - 1)

    def define_phi_rule_nonnegativity(model, i, j, h):
        return model.phi[i, j, h] >= 0

    model.define_phi_upper_xi = Constraint(model.i, model.j, model.h, rule=define_phi_rule_upper_xi)
    model.define_phi_upper_y = Constraint(model.i, model.j, model.h, rule=define_phi_rule_upper_y)
    model.define_phi_lower_y_xi = Constraint(model.i, model.j, model.h, rule=define_phi_rule_lower_y_xi)
    model.define_phi_nonnegativity = Constraint(model.i, model.j, model.h, rule=define_phi_rule_nonnegativity)

    # Define omega constraints
    def define_omega_rule_lower_bound(model, i, j):
        return model.omega[i, j] >= 0

    def define_omega_rule_first_lower_bound(model, i, j):
        return model.omega[i, j] >= (2**(-H)) * model.y[j] + model.ep[i] - (2**(-H))

    def define_omega_rule_upper_bound_y(model, i, j):
        return model.omega[i, j] <= (2**(-H)) * model.y[j]

    def define_omega_rule_upper_bound_ep(model, i, j):
        return model.omega[i, j] <= model.ep[i]

    model.define_omega_lower_bound = Constraint(model.i, model.j, rule=define_omega_rule_lower_bound)
    model.define_omega_first_lower_bound = Constraint(model.i, model.j, rule=define_omega_rule_first_lower_bound)
    model.define_omega_upper_bound_y = Constraint(model.i, model.j, rule=define_omega_rule_upper_bound_y)
    model.define_omega_upper_bound_ep = Constraint(model.i, model.j, rule=define_omega_rule_upper_bound_ep)

    # Limits on xi
    def xi_lim_rule_upper(model, h):
        return sum(model.xi[i, h] for i in model.i) <= 0.75 * n

    def xi_lim_rule_lower(model, h):
        return sum(model.xi[i, h] for i in model.i) >= 0.25 * n

    model.xi_lim_upper = Constraint(model.h, rule=xi_lim_rule_upper)
    model.xi_lim_lower = Constraint(model.h, rule=xi_lim_rule_lower)

    # Define w
    def w_definition_rule(model, i, j):
        return model.w[i, j] == sum(2**(-h-1) * model.phi[i, j, h] for h in model.h) + model.omega[i, j]

    model.w_definition = Constraint(model.i, model.j, rule=w_definition_rule)

    # Triangle area constraints
    def triangle_area_rule(model, i, j, k):
        if i < j and j < k:
            return model.S[i, j, k] == 0.5 * (
                model.w[i, j] - model.w[i, k] + model.w[j, k] - model.w[j, i] + model.w[k, i] - model.w[k, j])
        return Constraint.Skip

    def triangle_area_linearize_1_rule(model, i, j, k):
        if i < j and j < k:
            return (1 - model.b[i, j, k]) * (ub + 0.5) + model.S[i, j, k] >= model.z
        return Constraint.Skip

    def triangle_area_linearize_2_rule(model, i, j, k):
        if i < j and j < k:
            return model.b[i, j, k] * (ub + 0.5) - model.S[i, j, k] >= model.z
        return Constraint.Skip

    def triangle_area_upper_bound_rule(model, i, j, k):
        if i < j and j < k:
            return model.S[i, j, k] <= 0.5 * model.b[i, j, k]
        return Constraint.Skip

    def triangle_area_lower_bound_rule(model, i, j, k):
        if i < j and j < k:
            return model.S[i, j, k] >= 0.5 * (model.b[i, j, k] - 1)
        return Constraint.Skip

    model.triangle_area = Constraint(model.i, model.j, model.k, rule=triangle_area_rule)
    model.triangle_area_linearize_1 = Constraint(model.i, model.j, model.k, rule=triangle_area_linearize_1_rule)
    model.triangle_area_linearize_2 = Constraint(model.i, model.j, model.k, rule=triangle_area_linearize_2_rule)
    model.triangle_area_upper_bound = Constraint(model.i, model.j, model.k, rule=triangle_area_upper_bound_rule)
    model.triangle_area_lower_bound = Constraint(model.i, model.j, model.k, rule=triangle_area_lower_bound_rule)

    # Objective function
    model.objective = Objective(expr=model.z, sense=maximize)

    # Solve
    solver = SolverFactory('baron')
    results = solver.solve(model, tee=True)

    # Extract results
    optimal_z = model.z.value
    optimal_x = [model.w[i, n - 1].value / model.y[n - 1].value for i in range(1, n + 1)]
    optimal_y = [model.y[i].value for i in range(1, n + 1)]
    optimal_b = sum(model.b[i, j, k].value for i in range(1, n + 1) for j in range(i + 1, n + 1) for k in range(j + 1, n + 1))
    optimal_ep = [model.ep[i].value for i in range(1, n + 1)]

    return optimal_z, optimal_x, optimal_y, optimal_b, optimal_ep


def plot_solution(optimal_z, optimal_x, optimal_y):
    plt.figure(figsize=(10, 10))
    plt.scatter(optimal_x, optimal_y, c='red')

    n = len(optimal_x)
    for i in range(n):
        plt.annotate(f"{i}", (optimal_x[i], optimal_y[i]), textcoords="offset points", xytext=(5, 5), ha='center')

    x = optimal_x
    y = optimal_y
    z = optimal_z
    minarea = 100

    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                area = abs(0.5 * (x[i] * (y[j] - y[k]) + x[j] * (y[k] - y[i]) + x[k] * (y[i] - y[j])))
                if area <= minarea:
                    minarea = area

    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                area = abs(0.5 * (x[i] * (y[j] - y[k]) + x[j] * (y[k] - y[i]) + x[k] * (y[i] - y[j])))
                if area == minarea:
                    trianglex = [x[i], x[j], x[k]]
                    triangley = [y[i], y[j], y[k]]
                    plt.plot(trianglex, triangley, 'g-')
                    plt.fill(trianglex, triangley)

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


# Inputs
n = int(input('n: '))
ub = float(input('Upper bound ? '))
lb = float(input('Lower bound ? '))
if lb == 0:
    lb = math.log(n) / (n ** 2)

ML = [4, 6, 7, 10, 11]
if n >= 6 and n <= 10:
    M = ML[n - 6]
else:
    M = int(input('m: '))

YB = [[(0.1666, 0.6667), (0.1666, 0.8334), (0.3333, 0.8334), (0.3333, 1)],
      [(0.142857, 0.714286), (0.142857, 0.714286), (0.285714, 0.857143), (0.285714, 0.857143), (0.428571, 1)]]
if n >= 7 and n <= 8:
    yb = YB[n - 7]
else:
    yb = [(0, 1) for _ in range(n - 3)]

H = int(input('H : '))
optimal_z, optimal_x, optimal_y, optimal_b, optimal_ep = heilbronn_triangle_approach3_MILP(n, H, M, ub, lb, yb)

if optimal_x is not None and optimal_y is not None:
    plot_solution(optimal_z, optimal_x, optimal_y)