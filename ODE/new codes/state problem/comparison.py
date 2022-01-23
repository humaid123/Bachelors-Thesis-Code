from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt 
import numpy as np 

### no event time experiment
measures_implemented = False
direction = "up"
def model_with_if(_, y):
    (S, E, I, R) = y
    global measures_implemented, direction
    N = 37.741 * (10**6)
    alpha = 1.0/8.0
    if (direction == "up"):
        if (E > 25000):
            measures_implemented = True
            direction = "down"
    else:
        if (E < 10000):
            measures_implemented = False
            direction = "up"

    beta = 0.005 if measures_implemented else 0.9
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
    t_eval = np.linspace(0, 180, 181)

    global measures_implemented, direction
    measures_implemented = False
    direction = "up"
    sol = solve_ivp(model_with_if, tspan, y0, method="LSODA", 
            atol=atol, rtol=rtol, t_eval=t_eval)

    return (t_eval, sol.y, sol.nfev)

### experiment WITH event
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

def experiment_with_event(atol, rtol):
    N = 37.741 * (10**6)
    E0 = 103
    I0 = 1
    R0 = 0
    S0 = N - (E0 + I0 + R0)
    y0 = (S0, E0, I0, R0)
    t0 = 0

    sum_nfev = 0
    times = None
    res = None

    res = np.array([[], [], [], []])
    times = np.array([])
    t_initial = t0
    y_initial = y0
    measures_implemented = False
    while t_initial < 180:
        tspan = [t_initial, 180]
        sol = None
        if (measures_implemented):
            sol = solve_ivp(model_with_measures, tspan, y_initial, events=root_10000, 
                    dense_output=True, atol=atol, rtol=rtol, method="LSODA")
            measures_implemented = False
        else:
            sol = solve_ivp(model_no_measures, tspan, y_initial, events=root_25000, 
                    dense_output=True, atol=atol, rtol=rtol, method="LSODA")
            measures_implemented = True
        
        t_event = 180 if len(sol.t_events[0]) == 0 else sol.t_events[0][0]
        t_calc = np.arange(t_initial, t_event, 1)
        y_cal = sol.sol(t_calc)
        
        res = np.concatenate((res, y_cal), axis=1)
        times = np.concatenate((times, t_calc))

        t_initial = t_event
        last_index = len(sol.y[0]) - 1
        y_initial = [sol.y[0][last_index], sol.y[1][last_index], 
                sol.y[2][last_index], sol.y[3][last_index]]
        sum_nfev += sol.nfev

    return (times, res, sum_nfev)

### getting the data
# default
(times_default, res_default, nfev_default) = experiment_no_event(1e-6, 1e-3)
(times_sharp, res_sharp, nfev_sharp) = experiment_no_event(1e-12, 1e-12)
(times_solved, res_solved, nfev_solved) = experiment_with_event(1e-6, 1e-3)

plt.plot(times_default, res_default[1])
plt.plot(times_sharp, res_sharp[1])
plt.plot(times_solved, res_solved[1])
plt.ylabel("E(t)")
plt.xlabel('time')
plt.legend(['default', 'sharpest', 'solved'], shadow=True)
plt.show()

print(f"default = {nfev_default}")
print(f"sharpest = {nfev_sharp}")
print(f"solved = {nfev_solved}")
