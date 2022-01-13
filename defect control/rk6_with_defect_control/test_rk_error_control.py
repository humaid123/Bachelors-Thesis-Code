
from rk6 import rk_error_control
from math import exp, floor
import matplotlib.pyplot as plt

def create_t_eval(start, end, num_points = 100):
    res = [start]
    h = (end - start)/ (num_points - 1)

    for _ in range(num_points - 1):
        res.append(
            res[-1] + h
        )
    return res

t_span_3 = [0, 10]
y0_3 = [1]

def model3(t, y):
    return [(1/4)*y*(1 - y/20)]

def solution3(t):
    return [20 / ( 1 + 19 * exp(-x/4) ) for x in t]

t_eval = create_t_eval(t_span_3[0], t_span_3[1])
actual_solutions = solution3(t_eval)


for tol in [1e-3, 1e-6, 1e-9, 1e-12]:
    res, sol, first_deriv, derivs = rk_error_control(model3, t_span_3, y0_3[0], tol)
    computed_solutions = [sol(x) for x in t_eval]
    
    plt.figure()
    plt.plot(t_eval, actual_solutions, label="actual solution")
    plt.plot(t_eval, computed_solutions, label="computed solution")
    #for (x, y) in res:
    #    plt.axvline(x=x)
    plt.title(f"solutions for tol={tol}")
    plt.legend()
    plt.show()


    error = [(actual_solution - computed_solution) for actual_solution, computed_solution in zip(actual_solutions, computed_solutions)]
    plt.figure()
    plt.plot(t_eval, error)
    plt.title(f"error for tol={tol}")
    plt.show()


    actual_f_evals = [model3(x, solution3([x])[0] )[0] for x in t_eval]
    hb_prime_evals = [first_deriv(x) for x in t_eval]
    plt.figure()
    plt.plot(t_eval, actual_f_evals)
    plt.plot(t_eval, hb_prime_evals)
    plt.title(f"first derivative for tol={tol}")
    plt.show()

    defects = [abs(actual_f_eval - hb_prime_eval) for (actual_f_eval, hb_prime_eval) in zip(actual_f_evals, hb_prime_evals)]
    plt.figure()
    plt.plot(t_eval, defects)
    plt.title(f"global defect for tol={tol}")
    plt.show()
    
    # defect graphs
    # we only pick 10 of the derivs
    pick_every = (len(derivs) // 10) + 1
    # print("pick every", pick_every)
    plotted_derivs = []
    for i in range(0, len(derivs), pick_every):
        plotted_derivs.append(
            derivs[i]
        )



    plt.figure()
    for (x_i_minus_1, x_i_plus_1, hb) in plotted_derivs:
        num_points = 100
        pts_to_sample = create_t_eval(x_i_minus_1, x_i_plus_1, num_points)
        defects = []
        for pt in pts_to_sample:
            y = solution3([pt])[0]
            f_eval  = model3(pt, y)[0]
            hb_prime_eval = hb.prime(pt)
            defects.append(hb_prime_eval - f_eval)
        maximum_defect = max(defects)
        minimum_defect = min(defects)
        plot_vals = [(defect - minimum_defect) / (maximum_defect - minimum_defect) for defect in defects]
        #plt.plot(xs, defects, label=f"x_{str(x_i_minus_1)}_{str(x_i_plus_1)}")

        x_axis = [i/(num_points - 1) for i in range(num_points)]
        plt.plot(x_axis, plot_vals, label=f"x_{str(x_i_minus_1)}_{str(x_i_plus_1)}")
    plt.title("plot of defects")
    plt.xlabel("step scaled between 0 and 2")
    plt.ylabel('defect')
    # plt.legend()
    plt.show()

print("H CANNOT GET TOO SMALL OR WE LOSE ALL THE ADVANTAGES we wanted....")
print("it seems we stop getting any more accurated for 1e-3 or smaller step sizes....")

