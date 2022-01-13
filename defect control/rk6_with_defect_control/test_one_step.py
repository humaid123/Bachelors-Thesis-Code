from math import sqrt, exp, log10
import matplotlib.pyplot as plt


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
t_span = [0, 10]
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
            xn + C[i] * h, 
            yn + h * sigma_prod(A[i], k, 0, i) # sigma sums the product of the two array from start up to but excluding end
        )[0] # f returnes an array for each component, for now everytme we can f, we will just extract the first component

    yn_plus_1 = yn + h * sigma_prod(B, k, 0, n_stages) # sigma sums the product of the two array from start up to but excluding end

    # yn_plus_1_higher_order = yn + h * sigma_prod(B_HAT, k, 0, n_stages) 

    return (k, yn_plus_1) # , yn_plus_1_higher_order)


def perfect_rk4_step(f, xn, yn, f_start, h):
    # copied from https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods
    # checked versus https://www.codesansar.com/numerical-methods/runge-kutta-fourth-order-rk4-python-program.htm
    k1 = f_start
    k2 = f(xn + h/2, yn + h * k1/2)[0]
    k3 = f(xn + h/2, yn + h * k2/2)[0]
    k4 = f(xn + h  , yn + h * k3)[0]

    k = ( k1+ 2*k2 + 2*k3 +k4 ) / 6
    yn_plus_1 = yn + h * k

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

# testing if order of one_step is indeed 4
log_error_at_each_h = []
log_hs = []

the_hs = [1e0, 1e-1, 1e-2, 1e-3]
for h in the_hs:
    max_error = float('-inf')
    xn = 0.0
    yn = solution([xn])[0]
    f_start = model(xn, yn)[0]
    (k, y_plus_1_higher_order) = one_step(model, xn, yn, f_start, h)
    (k, y_plus_1_perfect) = perfect_rk4_step(model, xn, yn, f_start, h)

    print("TEST", y_plus_1_higher_order - y_plus_1_perfect)

    true_y_plus_1 = solution([ xn + h ])[0]
    # one_step => error = abs(y_plus_1_higher_order - true_y_plus_1)
    error = abs(y_plus_1_perfect - true_y_plus_1)
    max_error = max(max_error, error)

    print("h", h)
    print(error, max_error)
    
    print("max error", max_error)
    log_10_max_error = log10(max_error)
    print("log10 max error", log_10_max_error)
    log_error_at_each_h.append( log_10_max_error )
    log_hs.append( -log10(h) )
    print(log_error_at_each_h, log_hs)


plt.figure()
plt.plot(log_hs, log_error_at_each_h)
plt.show()