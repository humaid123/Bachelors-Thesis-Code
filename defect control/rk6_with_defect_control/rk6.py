from math import sqrt
from HB import HB

A = [
    [0, 0, 0, 0, 0, 0, 0, 0, 0], 
    [0.06, 0, 0, 0, 0, 0, 0, 0, 0], 
    [0.019239962962962962, 0.07669337037037037, 0, 0, 0, 0, 0, 0, 0], 
    [0.035975, 0, 0.107925, 0, 0, 0, 0, 0, 0], 
    [1.3186834152331484, 0, -5.042058063628562, 4.220674648395414, 0, 0, 0, 0, 0], 
    [-41.872591664327516, 0, 159.4325621631375, -122.11921356501003, 5.531743066200054, 0, 0, 0, 0], 
    [-54.430156935316504, 0, 207.06725136501848, -158.61081378459, 6.991816585950242, -0.018597231062203234, 0, 0, 0], 
    [-54.66374178728198, 0, 207.95280625538936, -159.2889574744995, 7.018743740796944, -0.018338785905045722, -0.0005119484997882099, 0, 0], 
    [0.03438957868357036, 0, 0, 0.2582624555633503, 0.4209371189673537, 4.40539646966931, -176.48311902429865, 172.36413340141507, 0]
]
C = [0, 0.06, 0.09593333333333333, 0.1439, 0.4973, 0.9725, 0.9995, 1, 1]
B = [0.03438957868357036, 0, 0, 0.2582624555633503, 0.4209371189673537, 4.40539646966931, -176.48311902429865, 172.36413340141507, 0]
B_HAT = [0.04301298296577121, 0, 0, 0.23882842561019763, 0.4493871915553917, 2.2956854086040193, -73.02457612433467, 70.96432878226597, 0.03333333333333333]
n_stages = 9

def sigma_prod(arr1, arr2, start, end):
    res = 0
    for i in range(len(start, end)):
        res += (  arr1[i] * arr2[i]  )

def one_step(func, xn, yn, f_start, h):
    # theory:
        # y_n_plus_1 = y_n + h * sigma( b[i] * k[i] ) # where i goes from 1 to n_stages
        # k[i] = f(tn + c[i] * h, yn + h * sigma( A[i][j]) * k[j] ) where i is the current stage, j goes from 0 to current_stage-1

    k = [0] * n_stages
    k[0] = f_start     # we assume it is precomputed => k[0] = f(xn, yn) and that k[-1] = f(xn_plus_1, yn_plus_1)

    for i in range(n_stages):
        k[i] = func(   
            xn + C[i], 
            yn + h * sigma_prod(A[i], k, 0, i) # sigma sums the product of the two array from start up to but excluding end
        )[0] # f returnes an array for each component, for now everytme we can f, we will just extract the first component

    yn_plus_1 = yn + h * sigma_prod(B, k, 0, n_stages) # sigma sums the product of the two array from start up to but excluding end

    yn_plus_1_higher_order = yn + h * sigma_prod(B_HAT, k, 0, n_stages) 

    return (k, yn_plus_1, yn_plus_1_higher_order)

def rk_error_control(fun, t_span, y0, tol):
    # theory:
        # we take step by step from the start in t_span[0] to the end in t_span[1]
        # each step is a rk_step as per one_step(), the size of the step, h depends on how the step taken satisfied the tolerance
        # for now: 
        #       if the error on the step is less than the tolerance, we accept, else we reject and half the step size to retake the step
        #       if the error satisfied tol/10 as well, we double the step size for the next step
    xn = t_span[0]
    xend = t_span[1]

    yn = y0
    f_start = fun(xn, yn)[0] # each time we call, the function 'fun', we will have to extract the first value as solve_ivp wants vectorised model functions
    
    res = [(xn, yn)]
    fn_s = [f_start]

    h = sqrt(tol)
    while xn < xend:
        (k, yn_plus_1, yn_plus_1_higher_order) = one_step(fun, xn, yn, f_start, h)

        error = abs(yn_plus_1_higher_order - yn_plus_1)

        if error < tol:
            # accept the step, by moving the x and the y
            xn = xn + h
            yn = yn_plus_1
            res.append( (xn, yn) )

            # we note that the pair I got from Jim Verner does not seem to follow k[-1] = f(xn_plus_1, yn_plus_1)
            # so for now, I do an additional function evaluation
            f_eval = fun(xn, yn)[0]
            print(k[-1] - f_eval)
            f_start = f_eval # k[-1]
            fn_s.append(f_start)

            if error < (tol / 10):
                h *= 2
        else:
            h /= 2
    
    return res

def first_step(fun, xn, yn, f_start, tol):
    h = sqrt(tol)
    error = float('inf')

    stricter_tol = tol / 100

    while error > stricter_tol:
        (k, yn_plus_1, yn_plus_1_higher_order) = one_step(fun, xn, yn, f_start, h)

        # we do error control on the first few steps until we get enough function evaluations to be able to do defect control
        error = abs(yn_plus_1_higher_order - yn_plus_1)
        
        # we control the error at the end of the step with a stricter tolerance for the first step
        if error < stricter_tol:
            # accept the step, by moving the x and the y
            xn = xn + h
            yn = yn_plus_1_higher_order

            # we note that the pair I got from Jim Verner does not seem to follow k[-1] = f(xn_plus_1, yn_plus_1)
            # so for now, I do an additional function evaluation
            f_eval = fun(xn, yn)[0]
            print(k[-1] - f_eval)
            f_start = f_eval # k[-1]
            
            if error < (tol / 10):
                h *= 2
        else:
            h /= 2
    return (xn, yn, f_start, h)

def rk_defect_control(fun, t_span, y0, tol):
    # theory:
        # we take step by step from the start in t_span[0] to the end in t_span[1]
        # each step is a rk_step as per one_step(), the size of the step, h depends on how the step taken satisfied the tolerance
        # for now: 
        #       if the error on the step is less than the tolerance, we accept, else we reject and half the step size to retake the step
        #       if the error satisfied tol/10 as well, we double the step size for the next step
    xn = t_span[0]
    xend = t_span[1]

    yn = y0
    f_start = fun(xn, yn)[0] # each time we call, the function 'fun', we will have to extract the first value as solve_ivp wants vectorised model functions
    
    res = [(xn, yn)]
    fn_s = [f_start]
    interps = [None]

    # we do strict error control on the first step
    # we update all the values of xn, yn, fstart based on what happened on the first step
    (xn, yn, f_start, h) = first_step(fun, xn, yn, f_start, tol)
    res.append( (xn, yn) )
    fn_s.append(f_start)

    while xn < xend:
        (k, yn_plus_1, yn_plus_1_higher_order) = one_step(fun, xn, yn, f_start, h)

        x_i, y_i = res[-1]
        x_i_minus_1, y_i_minus_1 = res[-2]
        x_i_plus_1, y_i_plus_1 = x_i + h, yn_plus_1_higher_order

        f_i = fn_s[-1]
        f_i_minus_1 = fn_s[-2]
        # we note that the pair I got from Jim Verner does not seem to follow k[-1] = f(xn_plus_1, yn_plus_1)
        # so for now, I do an additional function evaluation
        f_eval = fun(xn, yn)[0]
        print(k[-1] - f_eval)
        f_i_plus_1 = f_eval # k[-1]

        this_interp = HB(
            x_i_minus_1, x_i, x_i_plus_1,
            y_i_minus_1, f_i_minus_1,
            y_i, f_i,
            y_i_plus_1, f_i_plus_1 
        )

        # we test the interpolant to check if the Hermite Birkhoff conditions are met as intended
        print(this_interp.eval(x_i_minus_1) - y_i_minus_1)
        print(this_interp.eval(x_i) - y_i)
        print(this_interp.eval(x_i_plus_1) - y_i_plus_1)
        print(this_interp.prime(x_i_minus_1) - f_i_minus_1)
        print(this_interp.prime(x_i) - f_i)
        print(this_interp.prime(x_i_plus_1) - f_i_plus_1)

        h_i = x_i_plus_1 - x_i
        x_sample_1 = x_i + 0.35 * h_i
        defect_sample_1 = abs( 
            this_interp.prime(x_sample_1) - fun( x_sample_1, this_interp.eval(x_sample_1) ) 
        )

        x_sample_2 = x_i + 0.75 * h_i
        defect_sample_2 = abs(
            this_interp.prime(x_sample_2) - fun( x_sample_2, this_interp.eval(x_sample_2) )
        )

        max_defect = max(defect_sample_1, defect_sample_2)

        if max_defect < tol:
            # accept the step, by moving the x and the y
            xn = x_i_plus_1
            yn = y_i_plus_1
            res.append( (xn, yn) )

            fn_s.append(f_i_plus_1)

            if max_defect < (tol / 10):
                h *= 2
        else:
            h /= 2
    
    return res