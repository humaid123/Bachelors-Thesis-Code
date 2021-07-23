from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt 
import numpy as np 
from timeit import default_timer as timer

def model_with_if(t, y):
    (S, E, I, R) = y

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

(times_lsoda, res_lsoda, nfev_lsoda, elapsed_lsoda) = experiment_with_if('LSODA')
(times_rk45, res_rk45, nfev_rk45, elapsed_rk45) = experiment_with_if('RK45')
(times_bdf, res_bdf, nfev_bdf, elapsed_bdf) = experiment_with_if('BDF')
(times_radau, res_radau, nfev_radau, elapsed_radau) = experiment_with_if('Radau')
#(times_dop853, res_dop853, nfev_dop853, elapsed_dop853) = experiment_with_if('DOP853')
(times_rk23, res_rk23, nfev_rk23, elapsed_rk23) = experiment_with_if('RK23')

plt.plot(times_bdf, res_bdf[2])
plt.plot(times_lsoda, res_lsoda[2])
plt.plot(times_radau, res_radau[2])
plt.plot(times_rk45, res_rk45[2])
#plt.plot(times_dop853, res_dop853[2])
plt.plot(times_rk23, res_rk23[2])
plt.xlabel('time')
plt.legend(['bdf', 'lsoda', 'radau', 'rk45', 'dop853', 'rk23'], shadow=True)
plt.show()

print(f"lsoda & {nfev_lsoda} & {elapsed_lsoda} \\\\")
print(f"rk45 & {nfev_rk45} & {elapsed_rk45} \\\\")
print(f"bdf & {nfev_bdf} & {elapsed_bdf} \\\\")
print(f"radau & {nfev_radau} & {elapsed_radau} \\\\")
print(f"dop853 & {nfev_dop853} & {elapsed_dop853} \\\\")
print(f"rk23 & {nfev_rk23} & {elapsed_rk23} \\\\")


