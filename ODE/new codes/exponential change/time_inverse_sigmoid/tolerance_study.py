from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt 
import numpy as np 
from timeit import default_timer as timer
from math import log10, exp

tolerances = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7]

with_measures = 0.005
without_measures = 0.9
time_measures_implemented = 27

steepness_of_change = 10 # smaller means that the change to 0.005 is slower
method = "RK45" 
# method = "LSODA"

def inverse_sigmoid(x, time_start_decrease):
    e = exp(-steepness_of_change*(x - time_start_decrease))
    return (without_measures - with_measures) * e / (1 + e) + with_measures


### no event time experiment
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

def experiment_no_event(atol, rtol):
    tspan = [0, 95]
    E0 = 103
    I0 = 1
    R0 = 0
    N = 37.741 * (10**6)
    S0 = N - (E0 + I0 + R0)
    y0 = (S0, E0, I0, R0)
    t_eval = np.linspace(0, 95, 96)

    sol = None
    sum_time = 0
    sum_nfev = 0
    sum_njev = 0
    sum_nlu = 0
    num_iter = 10
    for i in range(num_iter):
        start = timer()
        sol = solve_ivp(model_inverse_sigmoid, tspan, y0, method=method, 
                atol=atol, rtol=rtol, t_eval=t_eval)
        end = timer()
        sum_time += end-start
        sum_nfev += sol.nfev
        sum_njev += sol.njev
        sum_nlu += sol.nlu

    return (t_eval, sol.y, sum_nfev/num_iter, sum_time/num_iter, 
            sum_njev/num_iter, sum_nlu/num_iter)

### experiment WITH event
def experiment_with_event(atol, rtol):
    tspan_before = [0, 27]
    tspan_after = [27, 95]

    t_eval_before = np.linspace(0, 27, 28)
    t_eval_after = np.linspace(27, 95, 95-27+1)

    E0 = 103
    I0 = 1
    R0 = 0
    N = 37.741 * (10**6)
    S0 = N - (E0 + I0 + R0)
    y0_before = (S0, E0, I0, R0)

    times = None
    solution = None
    sum_time = 0
    sum_nfev = 0
    sum_njev = 0
    sum_nlu = 0
    num_iter = 10

    for i in range(num_iter):
        start1 = timer()
        sol1 = solve_ivp(model_inverse_sigmoid, tspan_before, y0_before, 
                t_eval=t_eval_before, method=method, atol=atol, rtol=rtol)
        end1 = timer()
    
        last_index = len(sol1.y[0]) - 1
        y0_after = (sol1.y[0][last_index], sol1.y[1][last_index],
            sol1.y[2][last_index], sol1.y[3][last_index])
   
        start2 = timer()
        sol2 = solve_ivp(model_inverse_sigmoid, tspan_after, y0_after, 
                t_eval=t_eval_after, method=method, atol=atol, rtol=rtol)
        end2 = timer()

        sum_time += end2-start2 + end1-start1
        sum_nfev += sol1.nfev + sol2.nfev
        sum_njev += sol1.njev + sol2.njev
        sum_nlu += sol1.nlu + sol2.nlu
        times = np.concatenate((sol1.t, sol2.t))
        solution = np.concatenate((sol1.y, sol2.y), axis=1)

    return (times, solution, sum_nfev/num_iter, sum_time/num_iter, 
                sum_njev/num_iter, sum_nlu/num_iter)

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
    plt.plot(times, solution[2])
plt.xlabel('time')
plt.ylabel("I(t)")
plt.legend(tolerances, shadow=True)
plt.show()

plt.figure(2)
for tolerance in tolerances:
    name = "with_event_" + str(abs(log10(tolerance)))
    times = ans_with_event[name][0]
    solution = ans_with_event[name][1]
    plt.plot(times, solution[2])
plt.xlabel('time')
plt.ylabel("I(t)")
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

    print(f"{tolerance} & {nfev_no_event} & {nfev_with_event} \\\\")