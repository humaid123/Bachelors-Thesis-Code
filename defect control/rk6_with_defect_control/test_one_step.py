from math import sqrt, exp, log10


C = [0, 1/2, 1/2, 1]
B = [1/6, 1/3, 1/3, 1/6]

A = [ 
	  [0, 0, 0, 0],
	  [1/2, 0,   0, 0],
      [0,   1/2, 0, 0],
      [0,   0,   1, 0]
	]
n_stages = 4


## test problem
def model(t, y):
    return [(1/4)*y*(1-y/20)]

def solution(t):
    return [20/(1 + 19*exp(-x/4)) for x in t]
t_span = [0, 1000]
y0 = [1]
## end test problem

def sigma_prod(arr1, arr2, start, end):
    res = 0
    for i in range(start, end):
        res += (  arr1[i] * arr2[i]  )
    return res

def one_step(func, xn, yn, f_start, h):
    # theory:
        # y_n_plus_1 = y_n + h * sigma( b[i] * k[i] ) # where i goes from 1 to n_stages
        # k[i] = f(tn + c[i] * h, yn + h * sigma( A[i][j]) * k[j] ) where i is the current stage, j goes from 0 to current_stage-1

    k = [0] * n_stages
    k[0] = f_start     # we assume it is precomputed => k[0] = f(xn, yn) and that k[-1] = f(xn_plus_1, yn_plus_1)

    for i in range(1, n_stages):
        k[i] = func(   
            xn + C[i], 
            yn + h * sigma_prod(A[i], k, 0, i) # sigma sums the product of the two array from start up to but excluding end
        )[0] # f returnes an array for each component, for now everytme we can f, we will just extract the first component

    yn_plus_1 = yn + h * sigma_prod(B, k, 0, n_stages) # sigma sums the product of the two array from start up to but excluding end

    # yn_plus_1_higher_order = yn + h * sigma_prod(B_HAT, k, 0, n_stages) 

    return (k, yn_plus_1) # , yn_plus_1_higher_order)


def perfect_rk4_step(f, xn, yn, f_start, h):
    # copied from https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods
    k1 = f_start
    k2 = f(xn + h/2, yn + h * k1/2)[0]
    k3 = f(xn + h/2, yn + h * k2/2)[0]
    k4 = f(xn + h  , yn + h * k3)[0]

    yn_plus_1 = yn + (1/6) * h * (k1 + 2*k2 + 2*k3 + k4)

    return ([k1, k2, k3, k4], yn_plus_1)


def rk4_comparison(fun, t_span, y0, tol):
    xn = t_span[0]
    xend = t_span[1]

    yn = y0
    f_start = fun(xn, yn)[0] # each time we call, the function 'fun', we will have to extract the first value as solve_ivp wants vectorised model functions
    
    res = [(xn, yn)]
    fn_s = [f_start]

    # we take 100 steps
    # h = (xend-xn) / 100
    # I force it to take huge steps to check how the function changes with time
    h = (xend-xn) / 100

    while xn < xend:
        (k, yn_plus_1) = one_step(fun, xn, yn, f_start, h)
        (k_actual, yn_plus_1_actual) = perfect_rk4_step(fun, xn, yn, f_start, h)

        # testing if everything is correct
        print("testing if the step is correct")
        print(yn_plus_1 - yn_plus_1_actual)
        for i, _ in enumerate(k):
            print(k[i] - k_actual[i])
        print("end of testing if the step is correct\n\n")

        # accept the step, by moving the x and the y
        xn = xn + h
        yn = yn_plus_1
        res.append( (xn, yn) )

        # I use to think:
            # we note that the pair I got from Jim Verner does not seem to follow k[-1] = f(xn_plus_1, yn_plus_1)
            # so for now, I do an additional function evaluation
            # f_eval = fun(xn, yn)[0]
            # print("k[-1] and actual f_eval", k[-1] - f_eval)
            # f_start = f_eval # k[-1]
        # but it turns out that even python makes an additional function evaluation
        # so it seems f_start = fun(x_i_plus_1, y_i_plus_1) is the way to go
        f_start = fun(xn, yn)[0]
        fn_s.append(f_start)    
    return res



rk4_comparison(model, t_span, y0[0], 1e-6)

