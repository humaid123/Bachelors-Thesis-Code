from math import exp, log10
from rk4 import one_step
import matplotlib.pyplot as plt

def model(t, y):
    return [(1/4)*y*(1-y/20)]

def solution(t):
    return [20/(1 + 19*exp(-x/4)) for x in t]
t_span = [0, 10]
y0 = [1]

log_error_at_each_h = []
log_hs = []

the_hs = [1e1, 1e0,1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10]
for h in the_hs:
    max_error = float("-inf")
    for xn in [0, 0.1, 1, 2, 5, 4.633, 10]:
        yn = solution([xn])[0]
        f_start = model(xn, yn)[0]
        (k, y_plus_1, y_plus_1_higher_order) = one_step(model, xn, yn, f_start, h)

        true_y_plus_1 = solution([ xn + h ])[0]
        error = abs(y_plus_1_higher_order - true_y_plus_1)
        max_error = max(max_error, error)
        print(error, max_error)
    print("max error", max_error)
    log_10_max_error = log10(max_error)
    print("log10 max error", log_10_max_error)
    log_error_at_each_h.append( log_10_max_error )
    log_hs.append( -log10(h) )
    print(log_error_at_each_h, log_hs)


plt.figure()
plt.plot(log_hs, log_error_at_each_h)
plt.show()

print("we can see that the order of the method is indeed O(h ^ 6) SO IT CAN TAKE BIG STEPS IF NEEDED")