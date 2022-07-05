from math import exp, sqrt, pi
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from numpy import transpose

xa = -5
xb = 5
a = xa
b = xb

w1 = 0.0 # -1.0
w2 = 27.0 # 1.0
gamma = 9 # 0.2
xi1 = -4.0
xi2 =  2.0
eta = 0.2
l = 1.00 # 1


# def D_func(x):
#     D = (1.00/(exp(l * (w1 - x)) + 1.00)) - 1.00
#     D = (D + (1.00/(exp(l * (x - w2)) + 1.00))) * (1.00 - gamma)
#     D = D + gamma
#     return D

with_measures = 0.005
without_measures = 0.9
time_measures_implemented = 50
# steepness_of_change = 25 # smaller means that the change is smoother

# def time_beta_func(x):
#     if (x <= time_measures_implemented):
#         return without_measures # * exp(-0.01* (x - 27) ** 2)
#     else:
#         return (without_measures - with_measures) * exp(-steepness_of_change*(x - time_measures_implemented) ** 2) + with_measures

# def decrease(x, time_start_decrease, steepness_of_change):
#     return (without_measures - with_measures) * exp(-steepness_of_change*(x - time_start_decrease) ** 2) + with_measures

def sigmoid(x, time_start_increase, steepness_of_change):
    return (without_measures - with_measures) * ( 1 / (1 + exp(-steepness_of_change*(x - time_start_increase)))) + with_measures

def inverse_sigmoid(x, time_start_decrease, steepness_of_change):
    e = exp(-steepness_of_change*(x - time_start_decrease))
    return (without_measures - with_measures) * e / (1 + e) + with_measures


t_eval = [x/10 for x in range(1000)]
# decrease_eval = [inverse_sigmoid(t, time_measures_implemented) for t in t_eval]
# increase_eval =  [sigmoid(t, time_measures_implemented) for t in t_eval]

plt.figure()
for steepness_of_change in [0.1, 1, 10]:
    plt.plot(t_eval, [inverse_sigmoid(t, time_measures_implemented, steepness_of_change) for t in t_eval], label=f"a={steepness_of_change}")
# plt.title("inverse sigmoid with discontinuity at t=50 with steepness of change, a, at 0.1, 1, 10")
plt.xlabel("t")
plt.ylabel("beta(t) - inverse sigmoid")
plt.legend()
plt.show()


plt.figure()
for steepness_of_change in [0.1, 1, 10]:
    plt.plot(t_eval, [sigmoid(t, time_measures_implemented, steepness_of_change) for t in t_eval], label=f"a={steepness_of_change}")
# plt.title("sigmoid with discontinuity at t=50 with steepness of change, a, at 0.1, 1, 10")
plt.xlabel("t")
plt.ylabel("beta(t) - inverse sigmoid")
plt.legend()
plt.show()