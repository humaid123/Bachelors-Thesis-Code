from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt 
import numpy as np 
from timeit import default_timer as timer
from math import log10, exp

tolerances = [1e-1, 1e-2, 1e-4, 1e-6, 1e-7, 1e-8, 1e-10, 1e-11]

### no event time experiment
measures_implemented = False
direction = "up"
time_last_change = 0

# method="LSODA"
method="RK45"

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
            measures_implemented = True
            direction = "down"
            time_last_change = t
    else:
        if (E < 10000):
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

def experiment_no_event(atol, rtol):
    N = 37.741 * (10**6)
    E0 = 103
    I0 = 1
    R0 = 0 
    S0 = N - (E0 + I0 + R0)
    y0 = (S0, E0, I0, R0)
    tspan = [0, 180]
    t_eval = np.linspace(0, 180, 362)

    global measures_implemented, direction, time_last_change
    
    sol = None
    sum_time = 0
    sum_nfev = 0
    sum_njev = 0
    sum_nlu = 0
    num_iter = 10
    for i in range(num_iter):
        measures_implemented = False
        direction = "up"
        time_last_change = 0

        start = timer()
        sol = solve_ivp(model_with_if, tspan, y0, method=method, 
                atol=atol, rtol=rtol, t_eval=t_eval)
        end = timer()
        sum_time += end-start
        sum_nfev += sol.nfev
        sum_njev += sol.njev
        sum_nlu += sol.nlu

    return (t_eval, sol.y, sum_nfev/num_iter, sum_time/num_iter, 
            sum_njev/num_iter, sum_nlu/num_iter)

### experiment WITH event
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

def root_25000(_, y, time_changed):
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

def root_10000(_, y, time_changed):
    E = y[1]
    return E - 10000
root_10000.terminal = True

def experiment_with_event(atol, rtol):
    N = 37.741 * (10**6)
    E0 = 103
    I0 = 1
    R0 = 0
    S0 = N - (E0 + I0 + R0)
    y0 = (S0, E0, I0, R0)
    t0 = 0

    sum_time = 0
    sum_nfev = 0
    sum_njev = 0
    sum_nlu = 0
    num_iter = 10
    times = None
    res = None

    for i in range(num_iter):
        res = np.array([[], [], [], []])
        times = np.array([])
        t_initial = t0
        y_initial = y0
        measures_implemented = False
        while t_initial < 180:
            tspan = [t_initial, 180]
            sol, start, end = None, None, None
            if (measures_implemented):
                start = timer()
                sol = solve_ivp(model_with_measures, tspan, y_initial, events=root_10000, args=[t_initial],
                        dense_output=True, atol=atol, rtol=rtol, method=method)
                end = timer()
                measures_implemented = False
            else:
                start = timer()
                sol = solve_ivp(model_no_measures, tspan, y_initial, events=root_25000, args=[t_initial],
                        dense_output=True, atol=atol, rtol=rtol, method=method)
                end = timer()
                measures_implemented = True
            
            t_event = 180 if len(sol.t_events[0]) == 0 else sol.t_events[0][0]
            t_calc = np.arange(t_initial, t_event, 0.25)
            y_cal = sol.sol(t_calc)
            
            res = np.concatenate((res, y_cal), axis=1)
            times = np.concatenate((times, t_calc))

            t_initial = t_event
            last_index = len(sol.y[0]) - 1
            y_initial = [sol.y[0][last_index], sol.y[1][last_index], 
                    sol.y[2][last_index], sol.y[3][last_index]]
            sum_nfev += sol.nfev
            sum_njev += sol.njev
            sum_nlu += sol.nlu
            sum_time += (end - start)
    return (times, res, sum_nfev/num_iter, sum_time/num_iter, sum_njev/num_iter, sum_nlu/num_iter)

### getting the data

ans_no_event = {}
for tolerance in tolerances:
    name = "no_event_" + str(abs(log10(tolerance)))
    ans_no_event[name] = experiment_no_event(tolerance, tolerance)

ans_with_event = {}
for tolerance in tolerances:
    name = "with_event_" + str(abs(log10(tolerance)))
    ans_with_event[name] = experiment_with_event(tolerance, tolerance)

### plotting the data
plt.figure(1)
for tolerance in tolerances:
    name = "no_event_" + str(abs(log10(tolerance)))
    times = ans_no_event[name][0]
    solution = ans_no_event[name][1]
    plt.plot(times, solution[1])
plt.ylabel("E(t)")
plt.xlabel('time')
plt.legend(tolerances, shadow=True)
plt.show()

plt.figure(2)
for tolerance in tolerances:
    name = "with_event_" + str(abs(log10(tolerance)))
    times = ans_with_event[name][0]
    solution = ans_with_event[name][1]
    plt.plot(times, solution[1])
plt.ylabel("E(t)")
plt.xlabel('time')
plt.legend(tolerances, shadow=True)
plt.show()


### printing the efficiency data
for tolerance in tolerances:
    name = "no_event_" + str(abs(log10(tolerance)))
    nfev_no_event = ans_no_event[name][2]
    time_no_event = ans_no_event[name][3]
    njev_no_event = ans_no_event[name][4]
    nlu_no_event = ans_no_event[name][5]

    name = "with_event_" + str(abs(log10(tolerance)))
    nfev_with_event = ans_with_event[name][2]
    time_with_event = ans_with_event[name][3]
    njev_with_event = ans_with_event[name][4]
    nlu_with_event = ans_with_event[name][5]

    print(f"{tolerance} & {int(nfev_no_event)} & {int(nfev_with_event)} \\\\")
