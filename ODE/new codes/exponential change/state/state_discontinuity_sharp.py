
from scipy.integrate import solve_ivp
import numpy as np 
import matplotlib.pyplot as plt 
from math import floor, ceil, exp
from timeit import default_timer as timer

measures_implemented = False
direction = "up"
time_last_change = 0

with_measures = 0.005
without_measures = 0.9
steepness_of_change = 10 # smaller means that the change is smoother

def sigmoid(x, time_start_increase):
    return (without_measures - with_measures) * ( 1 / (1 + exp(-steepness_of_change*(x - time_start_increase)))) + with_measures

def inverse_sigmoid(x, time_start_decrease):
    e = exp(-steepness_of_change*(x - time_start_decrease))
    return (without_measures - with_measures) * e / (1 + e) + with_measures


def model_with_if(t, y):
    (S, E, I, R) = y
    global measures_implemented, direction, time_last_change

    N = 37.741 * (10**6)
    alpha = 1.0/8.0
    if (direction == "up"):
        if (E > 25000):
            #print("switched to True")
            measures_implemented = True
            direction = "down"
            time_last_change = t
    else:
        if (E < 10000):
            #print("switched to false")
            measures_implemented = False
            direction = "up"
            time_last_change = t

    beta = None
    if measures_implemented:
        beta = inverse_sigmoid(t, time_last_change)
    else:
        beta = sigmoid(t, time_last_change)

    gamma = 0.06
    mu = 0.01/365

    #print(E, beta)

    dSdt = mu*N - mu*S - (beta/N)*I*S
    dEdt = (beta/N)*I*S - alpha*E - mu*E
    dIdt = alpha*E - gamma*I - mu*I
    dRdt = gamma*I - mu*R

    return (dSdt, dEdt, dIdt, dRdt)

def model_no_measures(t, y, time_start_increase):
    (S, E, I, R) = y

    N = 37.741 * 10**6
    alpha = 1.0/8.0
    beta = sigmoid(t, time_start_increase)
    gamma = 0.06
    mu = 0.01/365

    dSdt = mu*N - mu*S - (beta/N)*I*S
    dEdt = (beta/N)*I*S - alpha*E - mu*E
    dIdt = alpha*E - gamma*I - mu*I
    dRdt = gamma*I - mu*R
    return (dSdt, dEdt, dIdt, dRdt)

def root_25000(_, y, _time_changed):
    E = y[1]
    return E - 25000
root_25000.terminal = True

def model_with_measures(t, y, time_start_decrease):
    (S, E, I, R) = y

    N = 37.741* (10**6)
    alpha = 1.0/8.0
    beta = inverse_sigmoid(t, time_start_decrease)
    gamma = 0.06
    mu = 0.01/365

    dSdt = mu*N - mu*S - (beta/N)*I*S
    dEdt = (beta/N)*I*S - alpha*E - mu*E
    dIdt = alpha*E - gamma*I - mu*I
    dRdt = gamma*I - mu*R

    return (dSdt, dEdt, dIdt, dRdt)

def root_10000(_, y, _time_changed):
    E = y[1]
    return E - 10000
root_10000.terminal = True

def experiment_with_if(method):
    N = 37.741 * (10**6)
    E0 = 103
    I0 = 1
    R0 = 0
    S0 = N - (E0 + I0 + R0)
    y0 = (S0, E0, I0, R0)
    tspan = [0, 180]
    t_eval = np.linspace(0, 180, 362)
    global measures_implemented, direction, time_last_change

    measures_implemented = False
    direction = "up"
    time_last_change = 0
    start = timer()
    sol = solve_ivp(model_with_if, tspan, y0,
            method=method, t_eval=t_eval, atol=1e-12, rtol=1e-12)
    end = timer()
    return t_eval, sol.y, sol.nfev, end-start

def high_accuracy(method):
    N = 37.741 * (10**6)
    E0 = 103
    I0 = 1
    R0 = 0
    S0 = N - (E0 + I0 + R0)
    y0 = (S0, E0, I0, R0)
    t0 = 0

    time_changes = []

    res = np.array([[], [], [], []])
    times = np.array([])
    t_initial = t0
    y_initial = y0
    measures_implemented = False
    count_number_evaluations = 0
    time_elapsed = 0
    t_final = 180
    while t_initial < t_final:
        tspan = [t_initial, t_final]
        sol, start, end = None, None, None
        if (measures_implemented):
            start = timer()
            sol = solve_ivp(model_with_measures, tspan, y_initial, args={t_initial}, events=root_10000, dense_output=True, method=method, atol=1e-12, rtol=1e-12)
            end = timer()
            measures_implemented = False
        else:
            start = timer()
            sol = solve_ivp(model_no_measures, tspan, y_initial, events=root_25000,  args={t_initial}, dense_output=True, method=method, atol=1e-12, rtol=1e-12)
            end = timer()
            measures_implemented = True
        t_event = t_final if len(sol.t_events[0]) == 0 else sol.t_events[0][0]
        t_calc = np.arange(t_initial, t_event, 0.25)
        y_cal = sol.sol(t_calc)
        res = np.concatenate((res, y_cal), axis=1)
        times = np.concatenate((times, t_calc))
        t_initial = t_event

        time_changes.append(t_initial)
        
        last_index = len(sol.y[0]) - 1
        y_initial = [sol.y[0][last_index], sol.y[1][last_index], sol.y[2][last_index], sol.y[3][last_index]]
        count_number_evaluations += sol.nfev
        time_elapsed += (end - start)
    return (times, res, count_number_evaluations, time_elapsed, time_changes)


(times_bdf_if, res_bdf_if, nfev_bdf_if, time_bdf_if) = experiment_with_if("BDF")
(times_lsoda_if, res_lsoda_if, nfev_lsoda_if, time_lsoda_if) = experiment_with_if("LSODA")
(times_radau_if, res_radau_if, nfev_radau_if, time_radau_if) = experiment_with_if('Radau')
(times_rk45_if, res_rk45_if, nfev_rk45_if, time_rk45_if) = experiment_with_if('RK45')
(times_dop853_if, res_dop853_if, nfev_dop853_if, time_dop853_if) = experiment_with_if('DOP853')
(times_rk23_if, res_rk23_if, nfev_rk23_if, time_rk23_if) = experiment_with_if('RK23')

(times_high, res_high, nfev_high, elapsed_high, time_changes) = high_accuracy("LSODA")


plt.figure(1)
plt.plot(times_bdf_if, res_bdf_if[1], label="bdf")
plt.plot(times_lsoda_if, res_lsoda_if[1], label="lsoda")
plt.plot(times_rk45_if, res_rk45_if[1], label="rk45")
plt.plot(times_dop853_if, res_dop853_if[1], label="dop853")
plt.plot(times_rk23_if, res_rk23_if[1], label="rk23")
plt.plot(times_radau_if, res_radau_if[1], label="radau")
# plt.plot(times_radau_if, res_radau_if[1])
# plt.plot(times_high, res_high[1], label="high")
# for time in time_changes:
#     plt.axvline(x=time)
plt.xlabel('time')
plt.ylabel("E(t)")
plt.legend()
plt.show()

# plt.figure(2)
# plt.plot(times_radau_if, res_radau_if[1], label="radau")
# plt.xlabel('time')
# plt.ylabel("E(t)")
# plt.legend()
# plt.show()


print(f"lsoda & {nfev_lsoda_if} & {time_lsoda_if}")
print(f"bdf & {nfev_bdf_if} & {time_bdf_if}")
print(f"radau & {nfev_radau_if} & {time_radau_if}")
print(f"rk45 & {nfev_rk45_if} & {time_rk45_if}")
print(f"dop853 & {nfev_dop853_if} & {time_dop853_if}")
print(f"rk23 & {nfev_rk23_if} & {time_rk23_if}")