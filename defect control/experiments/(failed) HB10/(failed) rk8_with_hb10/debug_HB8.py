

from math import exp, log2
import matplotlib.pyplot as plt
from HB8_second_scheme import HB
from rk8 import one_step

class Monitor:
    def __init__(self) -> None:
        self.different_values_alpha = set()
        self.different_values_beta = set()
        self.n_steps=0
        self.n_successful_steps=0
    def print(self):
        print("alpha values", list(self.different_values_alpha))
        print("beta values", list(self.different_values_beta))
        print("n_steps", self.n_steps)
        print("n_successful_steps", self.n_successful_steps)

def take_step(func, xn, yn, f_start, h):
    x_i_plus_1 = xn + h
    (k, yn_plus_1, yn_plus_1_higher_order) = one_step(func, xn, yn, f_start, h)
    return (x_i_plus_1, yn_plus_1)

def experiment(model, solution, t_span, y0):
    monitor = Monitor()

    alpha = 1
    beta = 1

    for disturbance in [1e-6, 1e-7, 1e-8, 1e-9]:
        the_hs = [4, 2, 1, 1/2, 1/(2**2), 1/(2**3), 1/(2**4), 1/(2**5)]
        for h in the_hs:
            for x0 in [0, 0.1, 1, 2, 5, 4.633, 10]:
                (x_i_minus_2, y_i_minus_2)       =    x0, solution([x0])[0]
                f_i_minus_2 = model(x_i_minus_2, y_i_minus_2)[0]

                (x_i_minus_1, y_i_minus_1)       =    x_i_minus_2 + alpha*h , solution([ x_i_minus_2 + alpha*h  ])[0]
                x_i_minus_1_rk, y_i_minus_1_rk   = take_step(model, x_i_minus_2, y_i_minus_2, f_i_minus_2, alpha*h)
                f_i_minus_1_rk = model(x_i_minus_1_rk, y_i_minus_1_rk)[0] 

                (x_i_minus_0_5, y_i_minus_0_5)   =    x_i_minus_1 + h, solution([ x_i_minus_1 + h ])[0]
                x_i_minus_0_5_rk, y_i_minus_0_5_rk   = take_step(model, x_i_minus_1_rk, y_i_minus_1_rk, f_i_minus_1_rk, h)
                f_i_minus_0_5_rk = model(x_i_minus_0_5_rk, y_i_minus_0_5_rk)[0] 

                (x_i, y_i)                       =    x_i_minus_0_5 + h, solution([ x_i_minus_0_5 + h ])[0]
                x_i_rk, y_i_rk   = take_step(model, x_i_minus_0_5_rk, y_i_minus_0_5_rk, f_i_minus_0_5_rk, h)
                f_i_rk = model(x_i_rk, y_i_rk)[0] 

                (x_i_plus_1, y_i_plus_1)         =    x_i + beta*h, solution([x_i + beta*h])[0]
                x_i_plus_1_rk, y_i_plus_1_rk   = take_step(model, x_i_rk, y_i_rk, f_i_rk, beta*h)
                f_i_plus_1_rk = model(x_i_plus_1_rk, y_i_plus_1_rk)[0] 

                f_i_minus_2   =   model(x_i_minus_2, y_i_minus_2)[0]
                f_i_minus_1   =   model(x_i_minus_1, y_i_minus_1)[0]
                f_i_minus_0_5 =   model(x_i_minus_0_5, y_i_minus_0_5)[0]
                f_i           =   model(x_i, y_i)[0]   
                f_i_plus_1    =   model(x_i_plus_1, y_i_plus_1)[0]

                this_hb = HB(
                    x_i_minus_2, x_i_minus_1, x_i, x_i_plus_1,
                    y_i_minus_2, f_i_minus_2,
                    y_i_minus_1, f_i_minus_1,
                    y_i, f_i,
                    y_i_plus_1, f_i_plus_1,
                    monitor
                )

                this_hb_rk = HB(
                    x_i_minus_2, x_i_minus_1_rk, x_i_rk, x_i_plus_1_rk,
                    y_i_minus_2, f_i_minus_2,
                    y_i_minus_1_rk, f_i_minus_1_rk,
                    y_i_rk, f_i_rk,
                    y_i_plus_1_rk, f_i_plus_1_rk,
                    monitor
                )

                x_i_minus_0_5_next_step = x_i + (x_i_plus_1 - x_i) / 2
                eval_next_step = this_hb.eval(x_i_minus_0_5_next_step)
                actual_next_step = solution([ x_i_minus_0_5_next_step ])[0]
                error = abs(eval_next_step - actual_next_step)


                x_i_minus_0_5_next_step_rk = x_i_rk + (x_i_plus_1_rk - x_i_rk) / 2
                eval_next_step_rk = this_hb_rk.eval(x_i_minus_0_5_next_step_rk)
                actual_next_step_rk = solution([ x_i_minus_0_5_next_step_rk ])[0]
                error_rk = abs(eval_next_step_rk - actual_next_step_rk)

                print(x0, h, error, error_rk, abs(y_i_minus_1_rk - y_i_minus_1),sep="\t")


def model1(t, y):
    return [(1/4)*y*(1-y/20)]

def solution1(t):
    return [20/(1 + 19*exp(-x/4)) for x in t]
t_span1 = [0, 10]
y0_1 = [1]

experiment(model1, solution1, t_span1, y0_1)