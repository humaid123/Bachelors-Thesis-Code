from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt 
import numpy as np 
from timeit import default_timer as timer

my_nfev = 0
count = []

def model_with_if(t, y):
    (S, E, I, R) = y
    global my_nfev, count
    count.append((t, my_nfev))
    my_nfev += 1

    N = 37.741 * (10**6)
    alpha = 1.0/8.0
    beta = 0.005
    if t < 27:
        beta = 0.9
    gamma = 0.06
    mu = 0.01/365

    dSdt = mu*N - mu*S - (beta/N)*I*S
    dEdt = (beta/N)*I*S - alpha*E - mu*E
    dIdt = alpha*E - gamma*I - mu*I
    dRdt = gamma*I - mu*R

    return (dSdt, dEdt, dIdt, dRdt)

def discontinuity_with(method, i):
    tspan = [0, 50]
    E0 = 103
    I0 = 1
    R0 = 0
    N = 37.741 * (10**6)
    S0 = N - (E0 + I0 + R0)
    y0 = (S0, E0, I0, R0)

    global my_nfev, count
    my_nfev = 0
    count = []
    sol = solve_ivp(model_with_if, tspan, y0, method=method)

    times_nfev = []
    cumulative_nfev = []
    for (t, my_nfev) in count:
        times_nfev.append(t)
        cumulative_nfev.append(my_nfev)
    plt.figure(i)
    plt.axvline(x=27, c="red")
    plt.plot(times_nfev, cumulative_nfev)
    plt.legend(["line at 27", f"{method} v/s discontinuity"])
    plt.show()
    return


discontinuity_with('LSODA', 0)
#discontinuity_with('DOP853', 1)
#discontinuity_with('RK45', 2)
#discontinuity_with('BDF', 3)
discontinuity_with('Radau', 4)
#discontinuity_with('RK23', 5)
