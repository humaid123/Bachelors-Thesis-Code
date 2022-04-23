from cProfile import label
from re import I
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt 
import numpy as np 
from timeit import default_timer as timer

def model(t, y):
    (S, E, I, R) = y

    N = 37.741 * (10**6)
    alpha = 1.0/8.0
    beta = 0.005
    # if t < 27:
    #     beta = 0.9
    gamma = 0.06
    mu = 0.01/365

    dSdt = mu*N - mu*S - (beta/N)*I*S
    dEdt = (beta/N)*I*S - alpha*E - mu*E
    dIdt = alpha*E - gamma*I - mu*I
    dRdt = gamma*I - mu*R

    return (dSdt, dEdt, dIdt, dRdt)

def experiment(method):
    tspan = [0, 95]
    E0 = 103
    I0 = 1
    R0 = 0
    N = 37.741 * (10**6)
    S0 = N - (E0 + I0 + R0)
    y0 = (S0, E0, I0, R0)
    t_eval = np.linspace(0, 95, 96)

    sol = solve_ivp(model, tspan, y0, method=method, t_eval=t_eval)

    plt.figure()
    # plt.plot(t_eval, sol.y[0], label="S(t)")
    plt.plot(t_eval, sol.y[1], label="E(t)")
    plt.plot(t_eval, sol.y[2], label="I(t)")
    plt.plot(t_eval, sol.y[3], label="R(t)")
    plt.legend()
    plt.xlabel("t")
    plt.ylabel("solution")
    plt.show()


experiment('LSODA')
# experiment('RK45')
# experiment('BDF')
# experiment('Radau')
# experiment('DOP853')
# experiment('RK23')



