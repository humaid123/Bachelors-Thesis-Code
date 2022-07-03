from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt 
import numpy as np 
from math import exp
from timeit import default_timer as timer

with_measures = 0.005
without_measures = 0.9
time_measures_implemented = 27
steepness_of_change = 1 # smaller means that the change to 0.005 is slower

def beta_func(t):

    if (t <= time_measures_implemented):
        return without_measures # * exp(-0.01* (x - 27) ** 2)
    else:
        return (without_measures-with_measures) * exp(-steepness_of_change*(t - time_measures_implemented + 1) ** 2) + with_measures


def model_with_if(t, y):
    (S, E, I, R) = y

    N = 37.741 * (10**6)
    alpha = 1.0/8.0
    beta = beta_func(t)
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
    sol = solve_ivp(model_with_if, tspan, y0, method=method, t_eval=t_eval)
    end = timer()
    return (t_eval, sol.y, sol.nfev, end-start)


def model_before(t, y):
    def beta_before(t):
        return without_measures # * exp(-0.01* (x - 27) ** 2)

    (S, E, I, R) = y

    N = 37.741 * (10**6)
    alpha = 1.0/8.0
    beta = beta_before(t)
    gamma = 0.06
    mu = 0.01/365

    dSdt = mu*N - mu*S - (beta/N)*I*S
    dEdt = (beta/N)*I*S - alpha*E - mu*E
    dIdt = alpha*E - gamma*I - mu*I
    dRdt = gamma*I - mu*R

    return (dSdt, dEdt, dIdt, dRdt)

def model_after(t, y):
    def beta_after(t):
        return (without_measures-with_measures) * exp(-steepness_of_change*(t - time_measures_implemented + 1) ** 2) + with_measures

    (S, E, I, R) = y

    N = 37.741 * (10**6)
    alpha = 1.0/8.0
    beta = beta_after(t)
    gamma = 0.06
    mu = 0.01/365

    dSdt = mu*N - mu*S - (beta/N)*I*S
    dEdt = (beta/N)*I*S - alpha*E - mu*E
    dIdt = alpha*E - gamma*I - mu*I
    dRdt = gamma*I - mu*R

    return (dSdt, dEdt, dIdt, dRdt)


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
    sol1 = solve_ivp(model_before, tspan, y0, t_eval=t_eval, atol=1e-12, rtol=1e-12)
    end1 = timer()

    tspan = [time_measures_implemented, 95]
    last_index = len(sol1.y[0]) - 1
    y0 = (sol1.y[0][last_index],
            sol1.y[1][last_index],
            sol1.y[2][last_index],
            sol1.y[3][last_index])
    t_eval = np.linspace(time_measures_implemented, 95, 95-time_measures_implemented+1)
    start2 = timer()
    sol2 = solve_ivp(model_after, tspan, y0, t_eval=t_eval, atol=1e-12, rtol=1e-12)
    end2 = timer()
    res = (np.concatenate((sol1.t, sol2.t)), 
            np.concatenate((sol1.y, sol2.y), axis=1), 
            sol1.nfev + sol2.nfev, end1-start1 + end2-start2)
    return res

(times_lsoda, res_lsoda, nfev_lsoda, elapsed_lsoda) = experiment_with_if('LSODA')
(times_rk45, res_rk45, nfev_rk45, elapsed_rk45) = experiment_with_if('RK45')
(times_bdf, res_bdf, nfev_bdf, elapsed_bdf) = experiment_with_if('BDF')
(times_radau, res_radau, nfev_radau, elapsed_radau) = experiment_with_if('Radau')
(times_dop853, res_dop853, nfev_dop853, elapsed_dop853) = experiment_with_if('DOP853')
(times_rk23, res_rk23, nfev_rk23, elapsed_rk23) = experiment_with_if('RK23')
(times_high, res_high, nfev_high, elapsed_high) = high_accuracy()

plt.plot(times_bdf, res_bdf[2])
plt.plot(times_lsoda, res_lsoda[2])
plt.plot(times_radau, res_radau[2])
plt.plot(times_rk45, res_rk45[2])
plt.plot(times_dop853, res_dop853[2])
plt.plot(times_rk23, res_rk23[2])
plt.plot(times_high, res_high[2])
plt.xlabel('time')
plt.ylabel("I(t)")
plt.legend(['bdf', 'lsoda', 'radau', 'rk45', 'dop853', 'rk23', "high"], shadow=True)
plt.show()

print(f"lsoda & {nfev_lsoda} & {elapsed_lsoda} \\\\")
print(f"rk45 & {nfev_rk45} & {elapsed_rk45} \\\\")
print(f"bdf & {nfev_bdf} & {elapsed_bdf} \\\\")
print(f"radau & {nfev_radau} & {elapsed_radau} \\\\")
print(f"dop853 & {nfev_dop853} & {elapsed_dop853} \\\\")
print(f"rk23 & {nfev_rk23} & {elapsed_rk23} \\\\")
print(f"high & {nfev_high} & {elapsed_high} \\\\")

