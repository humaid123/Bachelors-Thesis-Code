from math import exp, log2
from case2_rk4 import one_step
import matplotlib.pyplot as plt

def model(t, y):
    return [(1/4)*y*(1-y/20)]

def solution(t):
    return [20/(1 + 19*exp(-x/4)) for x in t]
t_span = [0, 10]
y0 = [1]

log_error_at_each_h = []
log_hs = []

the_hs = [2, 1, 1/2, 1/(2**2), 1/(2**3), 1/(2**4), 1/(2**5), 1/(2**6), 1/(2**7), 1/(2**8), 1/(2**9), 1/(2**10), 1/(2**11)]
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
    log_max_error = log2(max_error)
    log_error_at_each_h.append( log_max_error )
    log_hs.append( -log2(h) )
    print(log_error_at_each_h, log_hs)


plt.figure()
plt.plot(log_hs, log_error_at_each_h)
plt.show()

print("we can see that the order of the method is indeed O(h ^ 6) SO IT CAN TAKE BIG STEPS IF NEEDED")

maximum = 0
for i in range(len(log_error_at_each_h) - 1):
    value = log_error_at_each_h[i] - log_error_at_each_h[i + 1]
    print("value", value)
    maximum = max(maximum, value)
print("macimum", maximum)