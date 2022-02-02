
from math import exp, log2
import matplotlib.pyplot as plt
from HB4 import HB

def create_t_eval(start, end, num_points = 100):
    res = [start]
    h = (end - start)/(num_points - 1)

    for _ in range(num_points - 1):
        res.append(
            res[-1] + h
        )
    return res

class Monitor:
    def __init__(self) -> None:
        self.different_values_alpha = set()
        self.n_steps=0
        self.n_successful_steps=0
    def print(self):
        print("alpha values", list(self.different_values_alpha))
        print("n_steps", self.n_steps)
        print("n_successful_steps", self.n_successful_steps)

def experiment(model, solution, t_span, y0):
    monitor = Monitor()
    the_hs = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8]
    for h in the_hs:
        # the_xs = [1, 2, 5, 4.633, 9]
        the_xs = [1]
        for x_i in the_xs:
            x_i_plus_1 = x_i + h

            y_i         = solution([ x_i ])[0]
            y_i_plus_1  = solution([ x_i_plus_1 ])[0]

            f_i         = model(x_i        , y_i)[0]
            f_i_plus_1  = model(x_i_plus_1 , y_i_plus_1)[0]


            this_hb = HB( 
                x_i, x_i_plus_1, 
                y_i, f_i, 
                y_i_plus_1, f_i_plus_1
            )

            t_eval = create_t_eval(x_i, x_i_plus_1)
            defects = []
            errors = []
            errors_horner = []
            defects_horner = []
            for x_eval in t_eval:
                actual_solution = solution([x_eval])[0]
                calculated_solution = this_hb.eval(x_eval)
                calculated_solution_horner = this_hb.eval_horner(x_eval)
                actual_prime = model(x_eval, actual_solution)[0]
                calculated_prime = this_hb.prime(x_eval)
                calculated_prime_horner = this_hb.prime_horner(x_eval)

                error = abs(actual_solution - calculated_solution)
                error_horner = abs(actual_solution - calculated_solution_horner)

                defect = abs(actual_prime - calculated_prime) 
                defect_horner = abs(actual_prime - calculated_prime_horner) 

                defects.append(defect)
                defects_horner.append(defect_horner)

                errors.append(error)  
                errors_horner.append(error_horner)  

               
            plt.figure()
            plt.plot(t_eval, defects)
            plt.title(f"defects at h={h}")
            plt.xlabel("from x_i to x_i_plus_1")
            plt.ylabel("defect")
            plt.show()

            plt.figure()
            plt.plot(t_eval, defects_horner)
            plt.title(f"defects of horner form at h={h}")
            plt.xlabel("from x_i to x_i_plus_1")
            plt.ylabel("defect")
            plt.show()


            plt.figure()
            plt.plot(t_eval, errors)
            plt.title(f"errors at h={h}")
            plt.xlabel("from x_i to x_i_plus_1")
            plt.ylabel("error")
            plt.show()

            plt.figure()
            plt.plot(t_eval, errors_horner)
            plt.title(f"errors of horner form at h={h}")
            plt.xlabel("from x_i to x_i_plus_1")
            plt.ylabel("error")
            plt.show()

            


def model1(t, y):
    return [(1/4)*y*(1-y/20)]

def solution1(t):
    return [20/(1 + 19*exp(-x/4)) for x in t]
t_span1 = [0, 10]
y0_1 = [1]
experiment(model1, solution1, t_span1, y0_1)