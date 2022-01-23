

from scipy.integrate import solve_ivp
import numpy as np 
import matplotlib.pyplot as plt 
from math import floor, ceil
from timeit import default_timer as timer

measures_implemented = False
direction = "up"
def model_with_if(_, y):
    (S, E, I, R) = y
    global measures_implemented, direction
    N = 37.741 * (10**6)
    alpha = 1.0/8.0
    if (direction == "up"):
        if (E > 25000):
            #print("switched to True")
            measures_implemented = True
            direction = "down"
    else:
        if (E < 10000):
            #print("switched to false")
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

def experiment_with_if(method):
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
    start = timer()
    sol = solve_ivp(model_with_if, tspan, y0,
            method=method, t_eval=t_eval, atol=1e-12, rtol=1e-12)
    end = timer()
    return t_eval, sol.y, sol.nfev, end-start

(times_bdf_if, res_bdf_if, nfev_bdf_if, time_bdf_if) = experiment_with_if("BDF")
(times_lsoda_if, res_lsoda_if, nfev_lsoda_if, time_lsoda_if) = experiment_with_if("LSODA")
(times_radau_if, res_radau_if, nfev_radau_if, time_radau_if) = experiment_with_if('Radau')
(times_rk45_if, res_rk45_if, nfev_rk45_if, time_rk45_if) = experiment_with_if('RK45')
(times_dop853_if, res_dop853_if, nfev_dop853_if, time_dop853_if) = experiment_with_if('DOP853')
(times_rk23_if, res_rk23_if, nfev_rk23_if, time_rk23_if) = experiment_with_if('RK23')

plt.figure(1)
plt.plot(times_bdf_if, res_bdf_if[1])
plt.plot(times_lsoda_if, res_lsoda_if[1])
plt.plot(times_rk45_if, res_rk45_if[1])
plt.plot(times_dop853_if, res_dop853_if[1])
plt.plot(times_rk23_if, res_rk23_if[1])
plt.xlabel('time')
plt.ylabel("E(t)")
plt.legend(['bdf', 'lsoda', 'rk45', 'dop853', 'rk23'], shadow=True)
plt.show()

plt.figure(2)
plt.plot(times_radau_if, res_radau_if[1])
plt.xlabel('time')
plt.ylabel("E(t)")
plt.legend(['radau'], shadow=True)
plt.show()

print(f"lsoda & {nfev_lsoda_if} & {time_lsoda_if}")
print(f"bdf & {nfev_bdf_if} & {time_bdf_if}")
print(f"radau & {nfev_radau_if} & {time_radau_if}")
print(f"rk45 & {nfev_rk45_if} & {time_rk45_if}")
print(f"dop853 & {nfev_dop853_if} & {time_dop853_if}")
print(f"rk23 & {nfev_rk23_if} & {time_rk23_if}")

