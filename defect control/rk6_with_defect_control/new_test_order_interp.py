from math import exp, log10
import matplotlib.pyplot as plt

from HB import HB

def model(t, y):
    return [(1/4)*y*(1-y/20)]

def solution(t):
    return [20/(1 + 19*exp(-x/4)) for x in t]
t_span = [0, 10]
y0 = [1]

log_max_error_at_each_h = []
log_hs = []

the_hs = [1e1, 1e0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10]
for h in the_hs:
    max_all_errors_at_this_h = float("-inf")
    for x_i in [0, 0.1, 1, 2, 5, 4.633, 10]:
        x_i_plus_1 = x_i + h
        x_i_minus_1 = x_i - h

        y_i         = solution([ x_i ])[0]
        y_i_plus_1  = solution([ x_i_plus_1 ])[0]
        y_i_minus_1 = solution([ x_i_minus_1 ])[0]

        f_i         = model(x_i        , y_i)[0]
        f_i_plus_1  = model(x_i_plus_1 , y_i_plus_1)[0]
        f_i_minus_1 = model(x_i_minus_1, y_i_minus_1)[0]

        this_hb = HB(
            x_i_minus_1, x_i, x_i_plus_1,
            y_i_minus_1, f_i_minus_1,
            y_i, f_i,
            y_i_plus_1, f_i_plus_1
        )

        max_error = float("-inf")
        for trial in [0.1, 0.3, 0.7, 1.2, 1.5, 1.7]:
            x_eval = x_i_minus_1 + trial * h
            error = this_hb.eval(x_eval) - solution([x_eval])[0]
            print(error)
            max_error = max(max_error, error)

        max_all_errors_at_this_h = max(max_all_errors_at_this_h, max_error)

    print("max all errors at this h", max_all_errors_at_this_h)
    log_10_max_all_errors_at_this_h = log10( max_all_errors_at_this_h )
    print("log10 max all errors at this h", log_10_max_all_errors_at_this_h)

    log_max_error_at_each_h.append( log_10_max_all_errors_at_this_h )
    log_hs.append( -log10(h) )


plt.figure()
plt.plot(log_hs, log_max_error_at_each_h)
plt.show()

print("we can see that the defect is also of order 6 no matter at which x we chose")