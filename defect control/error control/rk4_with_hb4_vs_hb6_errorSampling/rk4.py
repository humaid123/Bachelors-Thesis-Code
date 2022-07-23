from math import sqrt
from HB6 import HB6, ContinuousSolution, create_defect_samplings
from HB4 import HB4

# http://people.math.sfu.ca/~jverner/RKV65.IIIXb.Efficient.00000144617.081204.RATOnWeb
C = [0, 1/2, 1/2, 1]
B = [1/6, 1/3, 1/3, 1/6]

A = [ 
	  [0, 0, 0, 0],
	  [1/2, 0,   0, 0],
      [0,   1/2, 0, 0],
      [0,   0,   1, 0]
	]
n_stages = 4

# ===================================================================================================================================
# taking one rk step with a Runge Kutta pair

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
            # BUG => I you need to multiply by h => xn + C[i], 
            xn + C[i] * h, 
            yn + h * sigma_prod(A[i], k, 0, i) # sigma sums the product of the two array from start up to but excluding end
        )[0] # f returnes an array for each component, for now everytme we can f, we will just extract the first component

    yn_plus_1 = yn + h * sigma_prod(B, k, 0, n_stages) # sigma sums the product of the two array from start up to but excluding end

    yn_plus_1_higher_order = yn_plus_1 # yn + h * sigma_prod(B_HAT, k, 0, n_stages) 

    return (k, yn_plus_1, yn_plus_1_higher_order)

class Monitor:
    def __init__(self) -> None:
        self.different_values_alpha = set()
        self.n_steps=0
        self.n_successful_steps=0
    def print(self):
        print("alpha values", list(self.different_values_alpha))
        print("n_steps", self.n_steps)
        print("n_successful_steps", self.n_successful_steps)

def rk_error_control_perfect_first_step(fun, t_span, y0, tol, solution, gauss_rule=5):
    xn, xend = t_span
    yn = y0
    f_start = fun(xn, yn)[0] # each time we call, the function 'fun', we will have to extract the first value as solve_ivp wants vectorised model functions
    
    res = [(xn, yn)]
    fn_s = [f_start]
    interps = []
    lower_order_interps = []

    # we do a perfect step for the one_step
    h = 1e-2 # HB6 V is at 1e-2. sqrt(tol)
    xn = xn + h
    yn = solution([xn])[0]
    f_start = fun(xn, yn)[0]
    res.append( (xn, yn) )
    fn_s.append(f_start)

    monitor6 = Monitor()
    while xn < xend:
        (k, yn_plus_1, yn_plus_1_higher_order) = one_step(fun, xn, yn, f_start, h)

        x_i, y_i = res[-1]
        x_i_minus_1, y_i_minus_1 = res[-2]
        x_i_plus_1, y_i_plus_1 = x_i + h, yn_plus_1_higher_order

        f_i = fn_s[-1]
        f_i_minus_1 = fn_s[-2]
        f_i_plus_1 = fun(x_i_plus_1, y_i_plus_1)[0]

        this_interp_hb6 = HB6(
            x_i_minus_1, x_i, x_i_plus_1,
            y_i_minus_1, f_i_minus_1,
            y_i, f_i,
            y_i_plus_1, f_i_plus_1,
            monitor6
        )

        this_interp_hb4 = HB4(
                x_i, x_i_plus_1,
                y_i, f_i,
                y_i_plus_1, f_i_plus_1 
        )

        monitor6.n_steps += 1

        """
        # we test the interpolant to check if the Hermite Birkhoff conditions are met as intended
        if (abs(this_interp.eval(x_i_minus_1) - y_i_minus_1))  > 1e-12: print("wrong y_i_minus_1", abs(this_interp.eval(x_i_minus_1) - y_i_minus_1))
        if (abs(this_interp.eval(x_i)         - y_i))          > 1e-12: print("wrong y_i",         abs(this_interp.eval(x_i)         - y_i))
        if (abs(this_interp.eval(x_i_plus_1)  - y_i_plus_1))   > 1e-12: print("wrong y_i_plus_1",  abs(this_interp.eval(x_i_plus_1)  - y_i_plus_1))

        if (abs(this_interp.prime(x_i_minus_1) - f_i_minus_1)) > 1e-12: print("wrong f_i_minus_1", abs(this_interp.prime(x_i_minus_1) - f_i_minus_1))
        if (abs(this_interp.prime(x_i)         - f_i))         > 1e-12: print("wrong f_i",         abs(this_interp.prime(x_i)         - f_i))
        if (abs(this_interp.prime(x_i_plus_1)  - f_i_plus_1))  > 1e-12: print("wrong f_i_plus_1",  abs(this_interp.prime(x_i_plus_1)  - f_i_plus_1))
        """

        # error control on [x_i to x_i_plus_1]
        h_i = x_i_plus_1 - x_i
        x_sample_1 = x_i + 0.6 * h_i
        error_sample1 = abs(
            this_interp_hb4.eval(x_sample_1) - this_interp_hb6.eval(x_sample_1)
        )

        x_sample_2 = x_i + 0.8 * h_i
        error_sample2 = abs(
            this_interp_hb4.eval(x_sample_2) - this_interp_hb6.eval(x_sample_2)
        )

        max_error_estimate = max(error_sample1, error_sample2)
        
        if max_error_estimate < tol:
            # accept the step, by moving the x and the y
            xn = x_i_plus_1
            yn = y_i_plus_1
            res.append( (xn, yn) )

            f_start = f_i_plus_1
            fn_s.append(f_start)

            monitor6.n_successful_steps += 1


            interps.append(this_interp_hb6)
            lower_order_interps.append(this_interp_hb4)

            if max_error_estimate < (tol / 10):
                h *= 2
        else:
            h /= 2

    print("tolerance=", tol)
    monitor6.print()
    print("================================\n")
    
    continuous_sol = ContinuousSolution()
    continuous_sol.extend(interps)

    lower_continuous_sol = ContinuousSolution()
    lower_continuous_sol.extend(lower_order_interps)
    return (
        res, 
        continuous_sol.eval,
        continuous_sol.prime,
        continuous_sol.create_error_samplings(),
        lower_continuous_sol,
        lower_continuous_sol.create_error_samplings(),
    )

# =================================================================================
# the following attempt is when the solver is to keep alpha at 1 throughout the integration
# def rk_error_control_static_alpha_beta(fun, t_span, y0, tol, solution, gauss_rule=5):
#     xn, xend = t_span
#     yn = y0
#     f_start = fun(xn, yn)[0] # each time we call, the function 'fun', we will have to extract the first value as solve_ivp wants vectorised model functions
    
#     res = [(xn, yn)]
#     fn_s = [f_start]
#     interps = []

#     # we do a perfect step for the one_step
#     h = 1e-2 # HB6 V is at 1e-2. sqrt(tol)
#     xn = xn + h
#     yn = solution([xn])[0]
#     f_start = fun(xn, yn)[0]
#     res.append( (xn, yn) )
#     fn_s.append(f_start)

#     xn = xn + h
#     yn = solution([xn])[0]
#     f_start = fun(xn, yn)[0]
#     res.append( (xn, yn) )
#     fn_s.append(f_start)

#     x_i_plus_1, y_i_plus_1   = res[-1]
#     x_i, y_i                 = res[-2]
#     x_i_minus_1, y_i_minus_1 = res[-3]

#     f_i_plus_1  = fn_s[-1]
#     f_i         = fn_s[-2]
#     f_i_minus_1 = fn_s[-3]

#     monitor6 = Monitor()

#     this_interp = HB6(
#             x_i_minus_1, x_i, x_i_plus_1,
#             y_i_minus_1, f_i_minus_1,
#             y_i, f_i,
#             y_i_plus_1, f_i_plus_1,
#             monitor6
#     )

#     continuous_sol = ContinuousSolution()
#     continuous_sol.append(this_interp)

#     while xn < xend:
#         (k, yn_plus_1, yn_plus_1_higher_order) = one_step(fun, xn, yn, f_start, h)

#         x_i, y_i = res[-1]
#         x_i_plus_1, y_i_plus_1 = x_i + h, yn_plus_1_higher_order

#         f_i = fn_s[-1]
#         f_i_plus_1 = fun(x_i_plus_1, y_i_plus_1)[0]


#         print("evaluating x_i_minus_1")
#         x_i_minus_1 = x_i - h
#         y_i_minus_1 = continuous_sol.eval(x_i_minus_1)
#         f_i_minus_1 = continuous_sol.prime(x_i_minus_1)

#         print("creating hb6")
#         this_interp_hb6 = HB6(
#             x_i_minus_1, x_i, x_i_plus_1,
#             y_i_minus_1, f_i_minus_1,
#             y_i, f_i,
#             y_i_plus_1, f_i_plus_1,
#             monitor6
#         )

#         print("creating hb4")
#         this_interp_hb4 = HB4(
#                 x_i, x_i_plus_1,
#                 y_i, f_i,
#                 y_i_plus_1, f_i_plus_1 
#         )

#         monitor6.n_steps += 1

#         """
#         # we test the interpolant to check if the Hermite Birkhoff conditions are met as intended
#         if (abs(this_interp.eval(x_i_minus_1) - y_i_minus_1))  > 1e-12: print("wrong y_i_minus_1", abs(this_interp.eval(x_i_minus_1) - y_i_minus_1))
#         if (abs(this_interp.eval(x_i)         - y_i))          > 1e-12: print("wrong y_i",         abs(this_interp.eval(x_i)         - y_i))
#         if (abs(this_interp.eval(x_i_plus_1)  - y_i_plus_1))   > 1e-12: print("wrong y_i_plus_1",  abs(this_interp.eval(x_i_plus_1)  - y_i_plus_1))

#         if (abs(this_interp.prime(x_i_minus_1) - f_i_minus_1)) > 1e-12: print("wrong f_i_minus_1", abs(this_interp.prime(x_i_minus_1) - f_i_minus_1))
#         if (abs(this_interp.prime(x_i)         - f_i))         > 1e-12: print("wrong f_i",         abs(this_interp.prime(x_i)         - f_i))
#         if (abs(this_interp.prime(x_i_plus_1)  - f_i_plus_1))  > 1e-12: print("wrong f_i_plus_1",  abs(this_interp.prime(x_i_plus_1)  - f_i_plus_1))
#         """

#         # error control on [x_i to x_i_plus_1]
#         print("error sample 1")
#         h_i = x_i_plus_1 - x_i
#         x_sample_1 = x_i + 0.1 * h_i
#         error_sample1 = abs(
#             this_interp_hb4.eval(x_sample_1) - this_interp_hb6.eval(x_sample_1)
#         )

#         print("error sample 2")
#         x_sample_2 = x_i + 0.9 * h_i
#         error_sample2 = abs(
#             this_interp_hb4.eval(x_sample_2) - this_interp_hb6.eval(x_sample_2)
#         )

#         max_error_estimate = max(error_sample1, error_sample2)
        
#         if max_error_estimate < tol:
#             # accept the step, by moving the x and the y
#             xn = x_i_plus_1
#             yn = y_i_plus_1
#             res.append( (xn, yn) )

#             f_start = f_i_plus_1
#             fn_s.append(f_start)

#             monitor6.n_successful_steps += 1


#             continuous_sol.append(this_interp_hb6)
#             if max_error_estimate < (tol / 10):
#                 h *= 2
#         else:
#             h /= 2

#     print("tolerance=", tol)
#     monitor6.print()
#     print("================================\n")
    
#     return (
#         res, 
#         continuous_sol.eval,
#         continuous_sol.prime,
#         continuous_sol.create_error_samplings()
#     )

