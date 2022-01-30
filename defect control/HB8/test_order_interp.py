from math import exp, log2
import matplotlib.pyplot as plt
from HB8 import HB

def experiment(model, solution, t_span, y0):
    log_max_error_at_each_h = []
    log_hs = []

    alpha = 1
    beta = 1

    the_hs = [2, 1, 1/2, 1/(2**2), 1/(2**3), 1/(2**4), 1/(2**5), 1/(2**6), 1/(2**7), 1/(2**8), 1/(2**9), 1/(2**10), 1/(2**11)]
    for h in the_hs:
        max_all_errors_at_this_h = float("-inf")
        for x_i in [0, 0.1, 1, 2, 5, 4.633, 10]:
            x_i_plus_1 = x_i + beta * h
            x_i_minus_1 = x_i - h
            x_i_minus_2 = x_i - (1 + alpha) * h

            y_i         = solution([ x_i ])[0]
            y_i_plus_1  = solution([ x_i_plus_1 ])[0]
            y_i_minus_1 = solution([ x_i_minus_1 ])[0]
            y_i_minus_2 = solution([ x_i_minus_2 ])[0]

            f_i         = model(x_i        , y_i)[0]
            f_i_plus_1  = model(x_i_plus_1 , y_i_plus_1)[0]
            f_i_minus_1 = model(x_i_minus_1, y_i_minus_1)[0]
            f_i_minus_2 = model(x_i_minus_2, y_i_minus_2)[0]

            this_hb = HB(
                x_i_minus_2, x_i_minus_1, x_i, x_i_plus_1,
                y_i_minus_2, f_i_minus_2,
                y_i_minus_1, f_i_minus_1,
                y_i, f_i,
                y_i_plus_1, f_i_plus_1
            )

            max_error = float("-inf")
            for trial in [0.1, 0.3, 0.7, 1.2, 1.5, 1.7]:
                x_eval = x_i_minus_1 + trial * h
                error = this_hb.eval(x_eval) - solution([x_eval])[0]
                max_error = max(max_error, error)

            max_all_errors_at_this_h = max(max_all_errors_at_this_h, max_error)

        print("max all errors at this h", max_all_errors_at_this_h)
        log_max_all_errors_at_this_h = log2( max_all_errors_at_this_h )
        print("log10 max all errors at this h", log_max_all_errors_at_this_h)

        log_max_error_at_each_h.append( log_max_all_errors_at_this_h )
        log_hs.append( -log2(h) )


    plt.figure()
    plt.plot(log_hs, log_max_error_at_each_h)
    plt.show()

    maximum = 0
    for i in range(len(log_max_error_at_each_h) - 1):
        value = log_max_error_at_each_h[i] - log_max_error_at_each_h[i + 1]
        print("value", value)
        maximum = max(maximum, value)
    print("maximum", maximum)

def model1(t, y):
    return [(1/4)*y*(1-y/20)]

def solution1(t):
    return [20/(1 + 19*exp(-x/4)) for x in t]
t_span1 = [0, 10]
y0_1 = [1]

experiment(model1, solution1, t_span1, y0_1)