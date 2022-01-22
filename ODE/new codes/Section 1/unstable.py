from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt 
import numpy as np 
from timeit import default_timer as timer
from math import log10

E0s = [70, 80, 90, 100, 110, 120]

### experiment WITH event
def model_before(t, y):
    (S, E, I, R) = y

    N = 37.741 * (10**6)
    alpha = 1.0/8.0
    beta = 0.9
    gamma = 0.06
    mu = 0.01/365

    dSdt = mu*N - mu*S - (beta/N)*I*S
    dEdt = (beta/N)*I*S - alpha*E - mu*E
    dIdt = alpha*E - gamma*I - mu*I
    dRdt = gamma*I - mu*R

    return (dSdt, dEdt, dIdt, dRdt)

def model_after(t, y):
    (S, E, I, R) = y

    N = 37.741 * (10**6)
    alpha = 1.0/8.0
    beta = 0.005
    gamma = 0.06
    mu = 0.01/365

    dSdt = mu*N - mu*S - (beta/N)*I*S
    dEdt = (beta/N)*I*S - alpha*E - mu*E
    dIdt = alpha*E - gamma*I - mu*I
    dRdt = gamma*I - mu*R

    return (dSdt, dEdt, dIdt, dRdt)

def stability_experiment(E0, t_final):
    I0 = 1
    R0 = 0
    N = 37.741 * (10**6)
    S0 = N - (E0 + I0 + R0)
    y0_before = (S0, E0, I0, R0)


    if (t_final < 40):
        tspan_before = [0, t_final]
        t_eval_before = np.linspace(0, t_final, t_final + 1)
        sol1 = solve_ivp(model_before, tspan_before, y0_before, 
            t_eval=t_eval_before, method="LSODA")
        return (sol1.t, sol1.y)
    else:
        tspan_before = [0, 40]
        tspan_after = [40, t_final]
        t_eval_before = np.linspace(0, 40, 41)
        t_eval_after = np.linspace(40, t_final, t_final-40+1)

        sol1 = solve_ivp(model_before, tspan_before, y0_before, 
            t_eval=t_eval_before, method="LSODA")
        last_index = len(sol1.y[0]) - 1
        y0_after = (sol1.y[0][last_index], sol1.y[1][last_index],
            sol1.y[2][last_index], sol1.y[3][last_index])
        sol2 = solve_ivp(model_after, tspan_after, y0_after, 
            t_eval=t_eval_after, method="LSODA")
        times = np.concatenate((sol1.t, sol2.t))
        solution = np.concatenate((sol1.y, sol2.y), axis=1)
        return (times, solution)

### getting the data

unstable = {}
for E0 in E0s:
    name = "E0=" + str(E0)
    unstable[name] = stability_experiment(E0, 39)

regain_stability = {}
for E0 in E0s:
    name = "E0=" + str(E0)
    regain_stability[name] = stability_experiment(E0, 120)

### plotting the data
plt.figure(1)
for E0 in E0s:
    name = "E0=" + str(E0)
    times = unstable[name][0]
    solution = unstable[name][1]
    plt.plot(times, solution[2])
plt.ylabel("I(t)")
plt.xlabel('time')
plt.legend(E0s, shadow=True)
plt.show()

plt.figure(1)
for E0 in E0s:
    name = "E0=" + str(E0)
    times = regain_stability[name][0]
    solution = regain_stability[name][1]
    plt.plot(times, solution[2])
plt.ylabel("I(t)")
plt.xlabel('time')
plt.legend(E0s, shadow=True)
plt.show()
