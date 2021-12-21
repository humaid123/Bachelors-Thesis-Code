
from rk6 import rk_error_control
from math import exp
import matplotlib.pyplot as plt


def create_t_eval(start, end, num_points = 100):
    res = [start]
    h = (end - start)/num_points

    for _ in range(num_points):
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
    res, sol = rk_error_control(model3, t_span_3, y0_3[0], tol)
    computed_solutions = [sol(x) for x in t_eval]

    plt.figure()
    plt.plot(t_eval, actual_solutions)
    plt.plot(t_eval, computed_solutions)
    for (x, y) in res:
        plt.axvline(x=x)
    plt.show()


    error = [(actual_solution - computed_solution) for actual_solution, computed_solution in zip(actual_solutions, computed_solutions)]
    plt.figure()
    plt.plot(t_eval, error)
    plt.show()
