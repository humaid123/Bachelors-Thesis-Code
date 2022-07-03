

from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt 
import numpy as np
from timeit import default_timer as timer
from math import exp

with_measures = 0.005
without_measures = 0.9
time_measures_implemented = 27
steepness_of_change = 25 # smaller means that the change to 0.005 is slower

def inverse_sigmoid(x, time_start_decrease):
    e = exp(-steepness_of_change*(x - time_start_decrease))
    return (without_measures - with_measures) * e / (1 + e) + with_measures


def model_inverse_sigmoid(t, y):
    (S, E, I, R) = y

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

def experiment(method):
    tspan = [0, time_measures_implemented]
    E0 = 103
    I0 = 1
    R0 = 0
    N = 37.741 * (10**6)
    S0 = N - (E0 + I0 + R0)
    y0 = (S0, E0, I0, R0)
    t_eval = np.linspace(0, time_measures_implemented, time_measures_implemented + 1)

    start1 = timer()
    sol1 = solve_ivp(model_inverse_sigmoid, tspan, y0, t_eval=t_eval, method=method)
    end1 = timer()

    tspan = [time_measures_implemented, 95]
    last_index = len(sol1.y[0]) - 1
    y0 = (sol1.y[0][last_index],
            sol1.y[1][last_index],
            sol1.y[2][last_index],
            sol1.y[3][last_index])
    t_eval = np.linspace(time_measures_implemented, 95, 95-time_measures_implemented+1)
    start2 = timer()
    sol2 = solve_ivp(model_inverse_sigmoid, tspan, y0, t_eval=t_eval, method=method)
    end2 = timer()
    res = (np.concatenate((sol1.t, sol2.t)), 
            np.concatenate((sol1.y, sol2.y), axis=1), 
            sol1.nfev + sol2.nfev, end1-start1 + end2-start2)
    return res

def high_accuracy():
    tspan = [0, time_measures_implemented]
    E0 = 103
    I0 = 1
    R0 = 0
    N = 37.741 * (10**6)
    S0 = N - (E0 + I0 + R0)
    y0 = (S0, E0, I0, R0)
    t_eval = np.linspace(0, time_measures_implemented, time_measures_implemented + 1)

    start1 = timer()
    sol1 = solve_ivp(model_inverse_sigmoid, tspan, y0, t_eval=t_eval, atol=1e-12, rtol=1e-12)
    end1 = timer()

    tspan = [time_measures_implemented, 95]
    last_index = len(sol1.y[0]) - 1
    y0 = (sol1.y[0][last_index],
            sol1.y[1][last_index],
            sol1.y[2][last_index],
            sol1.y[3][last_index])
    t_eval = np.linspace(time_measures_implemented, 95, 95-time_measures_implemented+1)
    start2 = timer()
    sol2 = solve_ivp(model_inverse_sigmoid, tspan, y0, t_eval=t_eval, atol=1e-12, rtol=1e-12)
    end2 = timer()
    res = (np.concatenate((sol1.t, sol2.t)), 
            np.concatenate((sol1.y, sol2.y), axis=1), 
            sol1.nfev + sol2.nfev, end1-start1 + end2-start2)
    return res

(times_lsoda, res_lsoda, nfev_lsoda, elapsed_lsoda) = experiment('LSODA')
(times_rk45, res_rk45, nfev_rk45, elapsed_rk45) = experiment('RK45')
(times_bdf, res_bdf, nfev_bdf, elapsed_bdf) = experiment('BDF')
(times_radau, res_radau, nfev_radau, elapsed_radau) = experiment('Radau')
(times_dop853, res_dop853, nfev_dop853, elapsed_dop853) = experiment('DOP853')
(times_rk23, res_rk23, nfev_rk23, elapsed_rk23) = experiment('RK23')
(times_high, res_high, nfev_high, elapsed_high) = high_accuracy()

plt.plot(times_bdf, res_bdf[2], label="bdf")
plt.plot(times_lsoda, res_lsoda[2], label="lsoda")
plt.plot(times_radau, res_radau[2], label="radau")
plt.plot(times_rk45, res_rk45[2], label="rk45")
plt.plot(times_dop853, res_dop853[2], label="dop853")
plt.plot(times_rk23, res_rk23[2], label="rk23")
plt.plot(times_high, res_high[2], label="high")
plt.xlabel('time')
plt.ylabel("I(t)")
plt.title(f"with disc hand - {steepness_of_change}")
plt.legend()
plt.show()

print(f"lsoda & {nfev_lsoda} & {elapsed_lsoda} \\\\")
print(f"rk45 & {nfev_rk45} & {elapsed_rk45} \\\\")
print(f"bdf & {nfev_bdf} & {elapsed_bdf} \\\\")
print(f"radau & {nfev_radau} & {elapsed_radau} \\\\")
print(f"dop853 & {nfev_dop853} & {elapsed_dop853} \\\\")
print(f"rk23 & {nfev_rk23} & {elapsed_rk23} \\\\")
print(f"high & {nfev_high} & {elapsed_high} \\\\")

