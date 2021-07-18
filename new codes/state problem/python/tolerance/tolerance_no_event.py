
# writing a dictionary as a json
# import json

# # as requested in comment
# exDict = {'exDict': exDict}

# with open('file.txt', 'w') as file:
#       file.write(json.dumps(exDict)) # use `json.loads` to do the reverse

from scipy.integrate import solve_ivp
import numpy as np 
import matplotlib.pyplot as plt 
from math import floor, ceil, log10
from timeit import default_timer as timer
import json

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

def experiment_with_if(method, atol, rtol):
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
    sol = solve_ivp(model_with_if, tspan, y0, method=method, 
            t_eval=t_eval, atol=atol, rtol=rtol)
    end = timer()
    return t_eval.tolist(), sol.y.tolist(), sol.nfev, end-start

ans = {}
tolerances = [1e-1, 1e-2, 1e-4, 1e-6, 1e-8, 1e-10, 1e-12, 1e-14, 1e-15]
for tolerance in tolerances:
    atol, rtol = tolerance, tolerance
    ans["lsoda_" + str(log10(tolerance))] = experiment_with_if("LSODA", atol, rtol)
    ans["bdf_" + str(log10(tolerance))] = experiment_with_if("BDF", atol, rtol)
    ans["radau_" + str(log10(tolerance))] = experiment_with_if('Radau', atol, rtol)
    ans["rk45_" + str(log10(tolerance))] = experiment_with_if('RK45', atol, rtol)
    ans["dop853_" + str(log10(tolerance))] = experiment_with_if('DOP853', atol, rtol)
    ans["rk23_" + str(log10(tolerance))] = experiment_with_if('RK23', atol, rtol)

exDict = {'state-tolerance': ans}

with open('file.txt', 'w') as file:
    file.write(json.dumps(exDict)) # use `json.loads` to do the reverse
