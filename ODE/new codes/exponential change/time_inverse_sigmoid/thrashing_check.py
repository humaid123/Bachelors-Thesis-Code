from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt 
import numpy as np 
from math import exp
from timeit import default_timer as timer

with_measures = 0.005
without_measures = 0.9
time_measures_implemented = 27
steepness_of_change = 0.1 # smaller means that the change to 0.005 is slower

thrashing_nfev = 0
thrashings = []

def inverse_sigmoid(x, time_start_decrease):
    e = exp(-steepness_of_change*(x - time_start_decrease))
    return (without_measures - with_measures) * e / (1 + e) + with_measures


def model_inverse_sigmoid(t, y):
    (S, E, I, R) = y

    global thrashing_nfev
    thrashing_nfev += 1
    thrashings.append( (t, thrashing_nfev) )

    N = 37.741 * (10**6)
    alpha = 1.0/8.0
    beta = inverse_sigmoid(t, time_measures_implemented)
    gamma = 0.06
    mu = 0.01/365

    dSdt = mu*N - mu*S - (beta/N)*I*S
    dEdt = (beta/N)*I*S - alpha*E - mu*E
    dIdt = alpha*E - gamma*I - mu*I
    dRdt = gamma*I - mu*R

    return (dSdt, dEdt, dIdt, dRdt)

def experiment_with_if(method):
    tspan = [0, 95]
    E0 = 103
    I0 = 1
    R0 = 0
    N = 37.741 * (10**6)
    S0 = N - (E0 + I0 + R0)
    y0 = (S0, E0, I0, R0)
    t_eval = np.linspace(0, 95, 96)

    start = timer()
    sol = solve_ivp(model_inverse_sigmoid, tspan, y0, method=method, t_eval=t_eval)
    end = timer()
    return (t_eval, sol.y, sol.nfev, end-start)

thrashing_nfev = 0
thrashings = []
(times_lsoda, res_lsoda, nfev_lsoda, elapsed_lsoda) = experiment_with_if('LSODA')

# plt.plot(times_lsoda, res_lsoda[2], label="lsoda")
# plt.xlabel('time')
# plt.ylabel("I(t)")
# plt.title(f"with no disc hand - {steepness_of_change}")
# plt.legend()
# plt.show()

times_thrashing = [thrashing[0] for thrashing in thrashings]
nfev_thrashing = [thrashing[1] for thrashing in thrashings]

plt.plot(times_thrashing, nfev_thrashing, label="nfev")
plt.xlabel('time')
plt.ylabel("cumulative nfev")
plt.title(f"CUMULATIVE NFEV")
plt.legend()
plt.show()

print(f"lsoda & {nfev_lsoda} & {elapsed_lsoda} \\\\")
