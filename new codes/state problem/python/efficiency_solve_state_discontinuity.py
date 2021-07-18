from scipy.integrate import solve_ivp
import numpy as np 
import matplotlib.pyplot as plt 
from math import floor, ceil
from timeit import default_timer as timer

def model_no_measures(_, y):
    (S, E, I, R) = y

    N = 37.741 * 10**6
    alpha = 1.0/8.0
    beta = 0.9
    gamma = 0.06
    mu = 0.01/365

    dSdt = mu*N - mu*S - (beta/N)*I*S
    dEdt = (beta/N)*I*S - alpha*E - mu*E
    dIdt = alpha*E - gamma*I - mu*I
    dRdt = gamma*I - mu*R
    return (dSdt, dEdt, dIdt, dRdt)

def root_25000(_, y):
    E = y[1]
    return E - 25000
root_25000.terminal = True

def model_with_measures(_, y):
    (S, E, I, R) = y

    N = 37.741* (10**6)
    alpha = 1.0/8.0
    beta = 0.005
    gamma = 0.06
    mu = 0.01/365

    dSdt = mu*N - mu*S - (beta/N)*I*S
    dEdt = (beta/N)*I*S - alpha*E - mu*E
    dIdt = alpha*E - gamma*I - mu*I
    dRdt = gamma*I - mu*R

    return (dSdt, dEdt, dIdt, dRdt)

def root_10000(_, y):
    E = y[1]
    return E - 10000
root_10000.terminal = True

def experiment(method):
    N = 37.741 * (10**6)
    E0 = 103
    I0 = 1
    R0 = 0
    S0 = N - (E0 + I0 + R0)
    y0 = (S0, E0, I0, R0)
    t0 = 0

    res = np.array([[], [], [], []])
    times = np.array([])
    t_initial = t0
    y_initial = y0
    measures_implemented = False
    count_number_evaluations = 0
    time_elapsed = 0
    while t_initial < 180:
        tspan = [t_initial, 180]
        sol, start, end = None, None, None
        if (measures_implemented):
            start = timer()
            sol = solve_ivp(model_with_measures, tspan, y_initial, events=root_10000, dense_output=True, method=method)
            end = timer()
            measures_implemented = False
        else:
            start = timer()
            sol = solve_ivp(model_no_measures, tspan, y_initial, events=root_25000, dense_output=True, method=method)
            end = timer()
            measures_implemented = True
        t_event = 180 if len(sol.t_events[0]) == 0 else sol.t_events[0][0]
        t_calc = np.arange(t_initial, t_event, 1)
        y_cal = sol.sol(t_calc)
        res = np.concatenate((res, y_cal), axis=1)
        times = np.concatenate((times, t_calc))
        t_initial = t_event
        last_index = len(sol.y[0]) - 1
        y_initial = [sol.y[0][last_index], sol.y[1][last_index], sol.y[2][last_index], sol.y[3][last_index]]
        count_number_evaluations += sol.nfev
        time_elapsed += (end - start)
    return (times, res, count_number_evaluations, time_elapsed)

(times_lsoda, res_lsoda, nfev_lsoda, elapsed_lsoda) = experiment("LSODA")
(times_rk45, res_rk45, nfev_rk45, elapsed_rk45) = experiment('RK45')
(times_bdf, res_bdf, nfev_bdf, elapsed_bdf) = experiment('BDF')
(times_radau, res_radau, nfev_radau, elapsed_radau) = experiment('Radau')
(times_dop853, res_dop853, nfev_dop853, elapsed_dop853) = experiment('DOP853')
(times_rk23, res_rk23, nfev_rk23, elapsed_rk23) = experiment('RK23')

plt.plot(times_lsoda, res_lsoda[1])
plt.plot(times_bdf, res_bdf[1])
plt.plot(times_radau, res_radau[1])
plt.plot(times_rk45, res_rk45[1])
plt.plot(times_dop853, res_dop853[1])
plt.plot(times_rk23, res_rk23[1])
plt.xlabel('time')
plt.legend(['lsoda', 'bdf', 'radau', 'rk45', 'dop853', 'rk23'], shadow=True)
plt.show()


print(f"lsoda & {nfev_lsoda} & {elapsed_lsoda}")
print(f"bdf & {nfev_bdf} & {elapsed_bdf}")
print(f"radau & {nfev_radau} & {elapsed_radau}")
print(f"rk45 & {nfev_rk45} & {elapsed_rk45}")
print(f"dop853 & {nfev_dop853} & {elapsed_dop853}")
print(f"rk23 & {nfev_rk23} & {elapsed_rk23}")

