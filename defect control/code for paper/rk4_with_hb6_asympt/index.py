# %%
from math import sin, sqrt, exp, cos
import matplotlib.pyplot as plt
from rk4 import (
    # rk_defect_control_static_alpha_func_call, 
    # rk_defect_control_perfect_first_step_smooth, 
    rk_defect_control_perfect_first_step
)
# %%
def create_t_eval(start, end, num_points = 100):
    res = [start]
    h = (end - start)/(num_points - 1)

    for _ in range(num_points - 1):
        res.append(
            res[-1] + h
        )
    return res

# %%
def experiment(model, y0, t_span, solution):
    t_eval = create_t_eval(t_span[0], t_span[1])
    tol = 1e-6
    (res, sol, first_deriv, derivs) = rk_defect_control_perfect_first_step(model, t_span, y0[0], tol, solution)

    # print("derivs", derivs)
    print("integration complete\n======================================================\n")

    # ====================================== figure of rk6 vs rk6_interps vs rk45
    # plt.figure()
    xs = [x[0] for x in res]
    # ys = [x[1] for x in res]
    # plt.plot(xs, ys, label="rk6")

    # plots of where the end of the steps occured to look at the interp
    # for this_x in xs:
    #     plt.axvline(x=this_x) 

    computed_solutions = [sol(x) for x in t_eval]
    # plt.plot(t_eval, computed_solutions, label="rk6_interpolated")

    actual_solutions = solution(t_eval)
    # plt.plot(t_eval, actual_solutions, label="solution")

    # # removed rk45 plt.title("solution vs rk45 vs rk6 vs rk6_interpolated")
    # plt.title("solution vs rk6_interpolated")
    # plt.xlabel("t")
    # plt.ylabel('y')
    # plt.legend(loc="upper right")
    # plt.show()
    # ====================================== end figure of rk6 vs rk6_interps vs rk45

    # ====================================== global error
    plt.figure()
    error = [abs(computed_solution - actual_solution) for (computed_solution, actual_solution) in zip(computed_solutions, actual_solutions)]
    # if tol > 1e-10:
    #     for this_x in xs:
    #         plt.axvline(x=this_x) 
    plt.plot(t_eval, error, label="global error")
    # plt.title("global error")
    # plt.xticks(xs)
    plt.xlabel("t")
    plt.ylabel("error")
    # plt.legend(loc="upper right")
    plt.show()
    # ====================================== end of global error

    # ====================================== figure of satisfying global defect
    actual_f_evals = [model(x, solution([x])[0])[0] for x in t_eval]
    hb_prime_evals = [first_deriv(x) for x in t_eval]
    # plt.figure()
    # plt.plot(t_eval, actual_f_evals)
    # plt.plot(t_eval, hb_prime_evals)
    # plt.title(f"first derivative for tol={tol}")
    # plt.show()

    defects = [abs(actual_f_eval - hb_prime_eval) for (actual_f_eval, hb_prime_eval) in zip(actual_f_evals, hb_prime_evals)]
    plt.figure()
    plt.plot(t_eval, defects, label="global defect")
    # plt.xticks(xs)    
    plt.ylabel("defect")
    plt.xlabel("t")
    # plt.title(f"global defect for tol={tol}")
    # plt.legend()    
    plt.show()
    
    # ====================================== end figure of satisfying global defect

    # ====================================== figure of defect shape

    minimum_step_size = 0.08
    defects_small_steps = []
    plt.figure()
    for (x_i_minus_1, x_i, x_i_plus_1, hb) in derivs:
        num_points = 100
        pts_to_sample = create_t_eval(x_i, x_i_plus_1, num_points)
        defects = []
        for pt in pts_to_sample:
            y = solution([pt])[0]
            f_eval  = model(pt, y)[0]
            hb_prime_eval = hb.prime(pt)
            defects.append( abs(hb_prime_eval - f_eval) )
        maximum_defect = max(defects)
        # minimum_defect = min(defects)
        # plot_vals = [(defect - minimum_defect) / (maximum_defect - minimum_defect) for defect in defects]
        plot_vals = [ defect/ maximum_defect for defect in defects]
        #plt.plot(xs, defects, label=f"x_{str(x_i_minus_1)}_{str(x_i_plus_1)}")
        x_axis = [i/(num_points - 1) for i in range(num_points)]

        if (x_i_plus_1 - x_i < minimum_step_size):
            defects_small_steps.append([x_axis, plot_vals])
            continue
        plt.plot(x_axis, plot_vals, label=f"x_{str(x_i_minus_1)}_{str(x_i_plus_1)}")
    # plt.title("plot of defects")
    plt.xlabel(r"$x_i$ to $x_{i+1}$")
    # plt.ylabel('defect/(max_defect on x_i to x_i_plus_1)')
    plt.ylabel('scaled defects')
    # plt.legend()
    plt.show()
    # ====================================== end figure of defect

    # ====================================== defects on small steps
    plt.figure()
    for [x_axis, plot_vals] in defects_small_steps:
        plt.plot(x_axis, plot_vals)
    # plt.title("plot of defects on small step sizes")
    plt.xlabel(r"$x_i$ to $x_{i+1}$")
    # plt.ylabel('defect/(max_defect on x_i to x_i_plus_1)')
    plt.ylabel('scaled defects')
    plt.show()


# %%
t_span_1 = [0, 10]
y0_1 = [1]

def model1(t, y):
    return [(-1/2) * y**3]

def solution1(t):
    return [1/sqrt(1+x) for x in t]

experiment(model1, y0_1, t_span_1, solution1)

# # %%
# t_span_2 = [0, 10]
# y0_2 = [1]

# def model2(t, y):
#     return [-2*t*y**2]

# def solution2(t):
#     return [1/(1+x**2) for x in t]

# experiment(model2, y0_2, t_span_2, solution2)

# %%
t_span_3 = [0, 10]
y0_3 = [1]

def model3(t, y):
    return [(1/4)*y*(1 - y/20)]

def solution3(t):
    return [20 / ( 1 + 19 * exp(-x/4) ) for x in t]

experiment(model3, y0_3, t_span_3, solution3)

# # %%
# t_span_4 = [0, 10]
# y0_4 = [0]

# def model4(t, y):
#     return [100 * (sin(t) - y)]
#     # return [10 * (sin(t) - y)]

# def solution4(t):
#     return [( 100 * ( exp(-100 * x) - cos(x) ) +  10000 * sin(x) ) / 10001 for x in t]
#     # return [( 10 * ( exp(-10 * x) - cos(x) ) +  100 * sin(x) ) / 101 for x in t]

# experiment(model4, y0_4, t_span_4, solution4)

# # %%
# t_span_5 = [0, 10]
# y0_5 = [2]

# def model5(t, y):
#     return [(15 * cos(10 * t))/y]

# def solution5(t):
#     return [sqrt(3*sin(10*x) + 4) for x in t]

# experiment(model5, y0_5, t_span_5, solution5)

# # %%
# t_span_6 = [0, 10]
# y0_6 = [1]

# def model6(t, y):
#     return [-y/(t+1)]

# def solution6(t):
#     return [1/(x+1) for x in t]

# experiment(model6, y0_6, t_span_6, solution6)

# %%
t_span_7 = [0, 10]
y0_7 = [1]

def model7(t, y):
    alpha = 0.1
    return [ -alpha*y - exp(-alpha*t)*sin(t)]

def solution7(t):
    alpha = 0.1
    return [exp(-alpha*x)*cos(x) for x in t]

experiment(model7, y0_7, t_span_7, solution7)

# # %%
# t_span_11 = [0, 10]
# y0_11 = [1]

# def model11(t, y):
#     return [-2*y + t]

# def solution11(t):
#     return [1/4 * (-1 + 5 * exp(-2 * x) + 2 * x) for x in t]

# experiment(model11, y0_11, t_span_11, solution11)

# # %%
# # THE PROBLEMS BELOW CANNOT BE DONE YET
# # CANNOT BE DONE AS MY CURRENT rk6 does not handle a vector for the ys
# ### ======================================================================


# # Jeff cash test set first one

# t_span_8 = [0, 10]
# eps = 0.1
# a = exp(-1/eps)
# y0_8 = [1, a/(eps*(-1+a))]

# def model8(t, y):
#     return [y[1], y[0]/eps]

# def solution8(t):
#     # THE experiment method calculates error on "computed[0]"
#     # so we can only verify the error of y[0] there
#     return [(1-exp(x/eps)*a)/(1-a) for x in t]

# experiment(model8, y0_8, t_span_8, solution8)

# ## the results were extremely bad. So i wanted to see the solution
# plt.figure()
# plt.plot(t_span_8, solution8(t_span_8))



# # %%
# # Jeff cash test set second one
# t_span_9 = [0, 10]
# eps = 0.1
# y0_9 = [1, -1/sqrt(eps)]

# def model9(t, y):
#     return [y[1], (y[0] + y[0]**2 - exp(-2*t/sqrt(eps)))/eps]

# def solution9(t):
#     # THE experiment method calculates error on "computed[0]"
#     # so we can only verify the error of y[0] there
#     return [exp(-x/sqrt(eps)) for x in t]

# experiment(model9, y0_9, t_span_9, solution9)

# # %%
# # Wolfram Alpha first problem

# t_span_10 = [0, 10]
# y0_10 = [1, 2]

# def model10(t, y):
#     return [y[1], -3*y[0] + 2*cos(4*t)]

# def solution10(t):
#     # THE experiment method calculates error on "computed[0]"
#     # so we can only verify the error of y[0] there
#     s = sqrt(3)
#     return [(26*s*sin(s*x) - 6*cos(4*x) + 45*cos(s*x))/39 for x in t]

# experiment(model10, y0_10, t_span_10, solution10)


