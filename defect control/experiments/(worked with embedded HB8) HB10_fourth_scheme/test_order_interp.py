from math import exp, log2
import matplotlib.pyplot as plt
from HB10_fourth_scheme import HB
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

def experiment(model, solution, t_span, y0):
    monitor = Monitor()
    log_max_error_at_each_h = []
    log_max_defect_at_each_h = []
    log_max_error_horner_at_each_h = []
    log_max_defect_horner_at_each_h = []
    log_hs = []

    alpha = 4
    beta = 1/4

    the_hs = [2, 1, 1/2, 1/(2**2), 1/(2**3), 1/(2**4), 1/(2**5), 1/(2**6), 1/(2**7), 1/(2**8), 1/(2**9), 1/(2**10), 1/(2**11)]
    for h in the_hs:
        max_all_errors_at_this_h = float("-inf")
        max_all_defects_at_this_h = float("-inf")
        max_all_errors_horner_at_this_h = float("-inf")
        max_all_defects_horner_at_this_h = float("-inf")
        for x0 in [0, 0.1, 1, 2, 5, 4.633, 10]:
            (x_i_minus_2, y_i_minus_2)       =    x0, solution([x0])[0]
            (x_i_minus_1, y_i_minus_1)       =    x_i_minus_2 + alpha*h , solution([ x_i_minus_2 + alpha*h  ])[0]
            (x_i_minus_0_5, y_i_minus_0_5)   =    x_i_minus_1 + h, solution([ x_i_minus_1 + h ])[0]
            (x_i, y_i)                       =    x_i_minus_0_5 + h, solution([ x_i_minus_0_5 + h ])[0]
            (x_i_plus_1, y_i_plus_1)         =    x_i + beta*h, solution([x_i + beta*h])[0]

            f_i_minus_2   =   model(x_i_minus_2, y_i_minus_2)[0]
            f_i_minus_1   =   model(x_i_minus_1, y_i_minus_1)[0]
            f_i_minus_0_5 =   model(x_i_minus_0_5, y_i_minus_0_5)[0]
            f_i           =   model(x_i, y_i)[0]   
            f_i_plus_1    =   model(x_i_plus_1, y_i_plus_1)[0]

            this_hb = HB(
        x_i_minus_2, x_i_minus_1, x_i_minus_0_5, x_i, x_i_plus_1,
        y_i_minus_2, f_i_minus_2,
        y_i_minus_1, f_i_minus_1,
        y_i_minus_0_5, f_i_minus_0_5,
        y_i, f_i,
        y_i_plus_1, f_i_plus_1,
        monitor
            )


            max_error = float("-inf")
            max_defect = float("-inf")
            max_error_horner = float("-inf")
            max_defect_horner = float("-inf")
            # beta can go up to 3 => for trial in [0.1, 0.3, 0.7, 1.2, 1.5, 1.7, 2, 2.3, 2.5, 2.7, 3]:
            # beta can go up to 1 => 
            for trial in [0.1, 0.2, 0.3, 0.4, 0.6, 0.7, 0.8, 1]:
                x_eval = x_i + trial * h
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

                max_error = max(max_error, error)
                max_defect = max(max_defect, defect)

                max_error_horner = max(max_error_horner, error_horner)
                max_defect_horner = max(max_defect_horner, defect_horner)
            max_all_errors_at_this_h = max(max_all_errors_at_this_h, max_error)
            max_all_defects_at_this_h = max(max_all_defects_at_this_h, max_defect)

            max_all_errors_horner_at_this_h = max(max_all_errors_horner_at_this_h, max_error_horner)
            max_all_defects_horner_at_this_h = max(max_all_defects_horner_at_this_h, max_defect_horner)


        log_max_error_at_each_h.append(  log2( max_all_errors_at_this_h ))
        log_max_defect_at_each_h.append( log2(max_all_defects_at_this_h) )

        log_max_error_horner_at_each_h.append(  log2( max_all_errors_horner_at_this_h ))
        log_max_defect_horner_at_each_h.append( log2(max_all_defects_horner_at_this_h) )

        log_hs.append( -log2(h) )


    plt.figure()
    plt.plot(log_hs, log_max_error_at_each_h)
    plt.show()

    plt.figure()
    plt.plot(log_hs, log_max_defect_at_each_h)
    plt.show()

    plt.figure()
    plt.plot(log_hs, log_max_error_horner_at_each_h)
    plt.show()

    plt.figure()
    plt.plot(log_hs, log_max_defect_horner_at_each_h)
    plt.show()

    maximum = 0
    for i in range(len(log_max_error_at_each_h) - 1):
        value = log_max_error_at_each_h[i] - log_max_error_at_each_h[i + 1]
        # print("value", value)
        maximum = max(maximum, value)
    print("order of method", maximum)
    maximum = 0
    for i in range(len(log_max_defect_at_each_h) - 1):
        value = log_max_defect_at_each_h[i] - log_max_defect_at_each_h[i + 1]
        # print("value", value)
        maximum = max(maximum, value)
    print("order of prime", maximum)

    maximum = 0
    for i in range(len(log_max_error_horner_at_each_h) - 1):
        value = log_max_error_horner_at_each_h[i] - log_max_error_horner_at_each_h[i + 1]
        # print("value", value)
        maximum = max(maximum, value)
    print("order of method in horner form", maximum)
    maximum = 0
    for i in range(len(log_max_defect_horner_at_each_h) - 1):
        value = log_max_defect_horner_at_each_h[i] - log_max_defect_horner_at_each_h[i + 1]
        # print("value", value)
        maximum = max(maximum, value)
    print("order of prime in horner form", maximum)


def model1(t, y):
    return [(1/4)*y*(1-y/20)]

def solution1(t):
    return [20/(1 + 19*exp(-x/4)) for x in t]
t_span1 = [0, 10]
y0_1 = [1]

experiment(model1, solution1, t_span1, y0_1)