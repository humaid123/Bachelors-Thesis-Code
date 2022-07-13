from math import sqrt
from CRK6Interp import CRK6ContinuousSolution, CRK6Interp
from CRK5Interp import CRK5ContinuousSolution, CRK5Interp

# http://people.math.sfu.ca/~jverner/RKV65.IIIXb.Efficient.00000144617.081204.RATOnWeb
A = [
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
    [0.06, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
    [0.019239962962962962, 0.07669337037037037, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
    [0.035975, 0, 0.107925, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
    [1.3186834152331484, 0, -5.042058063628562, 4.220674648395414, 0, 0, 0, 0, 0, 0, 0, 0], 
    [-41.872591664327516, 0, 159.4325621631375, -122.11921356501003, 5.531743066200054, 0, 0, 0, 0, 0, 0, 0], 
    [-54.430156935316504, 0, 207.06725136501848, -158.61081378459, 6.991816585950242, -0.018597231062203234, 0, 0, 0, 0, 0, 0], 
    [-54.66374178728198, 0, 207.95280625538936, -159.2889574744995, 7.018743740796944, -0.018338785905045722, -0.0005119484997882099, 0, 0, 0, 0, 0], 
    [0.03438957868357036, 0, 0, 0.2582624555633503, 0.4209371189673537, 4.40539646966931, -176.48311902429865, 172.36413340141507, 0, 0, 0, 0], 
    [0.016524159013572806, 0, 0, 0.3053128187514179, 0.2071200938201979, -1.293879140655123, 57.11988411588149, -55.87979207510932, 0.024830028297766014, 0, 0, 0], 
    [0.038150081818627744, 0, 0, 0.2502358252513705, 0.3249441447817608, 1.8224606658327962, -67.7137233269262, 66.03587911808127, -0.0363881087495127, 0.106441599909888, 0, 0], 
    [0.11178168039666012, 0, 0, 0.025757505109345213, 3.785140856363646, 92.34088993695727, -3819.461508432344, 3732.492711530704, -1.0756940209963033, -3.231539970732086, -4.7075390854586345, 0]
]
C = [0, 0.06, 0.09593333333333333, 0.1439, 0.4973, 0.9725, 0.9995, 1, 1, 0.5, 0.828, 0.28] 
B = [0.03438957868357036, 0, 0, 0.2582624555633503, 0.4209371189673537, 4.40539646966931, -176.48311902429865, 172.36413340141507, 0, 0, 0, 0]
B_HAT = [0.04301298296577121, 0, 0, 0.23882842561019763, 0.4493871915553917, 2.2956854086040193, -73.02457612433467, 70.96432878226597, 0.03333333333333333, 0, 0, 0]
n_stages = 12

def sigma_prod(arr1, arr2, start, end):
    res = 0
    for i in range(start, end):
        res += (  arr1[i] * arr2[i]  )
    return res

# ===================================================================================================================================
# taking one rk step with a Runge Kutta pair
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

    # B and B_HAT contains 0s at the positions for the extra stages
    yn_plus_1 = yn + h * sigma_prod(B, k, 0, n_stages) # sigma sums the product of the two array from start up to but excluding end
    yn_plus_1_higher_order = yn + h * sigma_prod(B_HAT, k, 0, n_stages) 

    return (k, yn_plus_1, yn_plus_1_higher_order)

# ===============================================================================================================================
# start of error control solver

def rk_error_control(fun, t_span, y0, tol):
    # theory:
        # we take step by step from the start in t_span[0] to the end in t_span[1]
        # each step is a rk_step as per one_step(), the size of the step, h depends on how the step taken satisfied the tolerance
        # for now: 
        #       if the error on the step is less than the tolerance, we accept, else we reject and half the step size to retake the step
        #       if the error satisfied tol/10 as well, we double the step size for the next step
    xn, xend = t_span

    yn = y0
    f_start = fun(xn, yn)[0] # each time we call, the function 'fun', we will have to extract the first value as solve_ivp wants vectorised model functions
    
    res = [(xn, yn)]
    fn_s = [f_start]

    continous_sol = CRK6ContinuousSolution()
    lower_continous_sol = CRK5ContinuousSolution()

    nsteps = 0
    n_successful_steps = 0

    h = sqrt(tol)
    while xn < xend:
        (k, yn_plus_1, yn_plus_1_higher_order) = one_step(fun, xn, yn, f_start, h)

        x_i = xn
        x_i_plus_1 = xn + h
        y_i = yn
        
        crk6 = CRK6Interp(k, y_i, x_i, x_i_plus_1)
        crk5 = CRK5Interp(k, y_i, x_i, x_i_plus_1)

        h_i = x_i_plus_1 - x_i
        x_sample_1 = x_i + 0.5 * h_i
        error_sample1 = abs(
            crk5.eval(x_sample_1) - crk6.eval(x_sample_1)
        )

        x_sample_2 = x_i + 0.8 * h_i
        error_sample2 = abs(
            crk5.eval(x_sample_2) - crk6.eval(x_sample_2)
        )

        max_error_estimate = max(error_sample1, error_sample2)

        nsteps += 1

        if max_error_estimate < tol:
            # accept the step, by moving the x and the y
            xn = x_i_plus_1
            yn = yn_plus_1
            res.append( (xn, yn) )

            # so for now, I do an additional function evaluation
            f_start = fun(xn, yn)[0]
            fn_s.append(f_start)

            n_successful_steps += 1

            # create interps and append to continuous solutionhere
            continous_sol.append(crk6)
            lower_continous_sol.append(crk5)
            if max_error_estimate < (tol / 10):
                h *= 2
        else:
            h /= 2
    print("nsteps =", nsteps)
    print("nsuccessful_steps =", n_successful_steps)
    print("=========================================================")
    return (
        res, 
        continous_sol,
        lower_continous_sol,
        continous_sol.create_error_sampling(),
        lower_continous_sol.create_error_sampling()
    )

# =================================================================================================================
# start of fixed step-size solver

# def rk_fixed_step(fun, t_span, y0, nsteps=100):
#     # theory:
#         # we take step by step from the start in t_span[0] to the end in t_span[1]
#         # each step is a rk_step as per one_step(), the size of the step, h, the step-size, is fixed based on the nubmer of steps the user wants to take
#     xn = t_span[0]
#     xend = t_span[1]

#     yn = y0
#     f_start = fun(xn, yn)[0] # each time we call, the function 'fun', we will have to extract the first value as solve_ivp wants vectorised model functions
    
#     res = [(xn, yn)]
#     fn_s = [f_start]

#     h = (xend - xn) / nsteps
#     while xn < xend:
#         (k, yn_plus_1, yn_plus_1_higher_order) = one_step(fun, xn, yn, f_start, h)

#         # error = abs(yn_plus_1_higher_order - yn_plus_1)
        
#         # accept the step, by moving the x and the y
#         xn = xn + h
#         yn = yn_plus_1_higher_order
#         res.append( (xn, yn) )
        
#         f_start = fun(xn, yn)[0] # we make a final function evalution at the current step
#         fn_s.append(f_start)

#     continous_sol = create_continuous_sol_from_results(res, fn_s)
#     return (
#         res, 
#         continous_sol.eval,
#         continous_sol.prime,
#         create_defect_samplings(res, fn_s)
#     )
