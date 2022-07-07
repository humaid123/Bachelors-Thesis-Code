from math import sqrt
from HB10 import HB, ContinuousSolution, create_defect_samplings, create_continuous_sol_from_results


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
        self.different_values_beta = set()
        self.n_steps=0
        self.n_successful_steps=0
    def print(self):
        print("alpha values", list(self.different_values_alpha))
        print("beta values", list(self.different_values_beta))
        print("n_steps", self.n_steps)
        print("n_successful_steps", self.n_successful_steps)
# ===============================================================================================
# RK defect control that allows the first step be perfect
# this allows us to do a successful first step...

def rk_defect_control_perfect_first_step(fun, t_span, y0, tol, solution):
    xn, xend = t_span
    yn = y0
    f_start = fun(xn, yn)[0] # each time we call, the function 'fun', we will have to extract the first value as solve_ivp wants vectorised model functions
    
    res = [(xn, yn)]
    fn_s = [f_start]
    interps = []

    # we do a perfect step for the one_step
    h = 1e-1 # as HB10 V is at 1e-1 sqrt(tol)
    for _ in range(4):
        xn = xn + h
        yn = solution([ xn ])[0]
        f_start = fun(xn, yn)[0]
        res.append( (xn, yn) )
        fn_s.append(f_start)

    x_i_plus_1, y_i_plus_1       = res[-1]
    x_i, y_i                     = res[-2]
    x_i_minus_0_5, y_i_minus_0_5 = res[-3]
    x_i_minus_1, y_i_minus_1     = res[-4]
    x_i_minus_2, y_i_minus_2     = res[-5]

    f_i_plus_1    = fn_s[-1]
    f_i           = fn_s[-2]
    f_i_minus_0_5 = fn_s[-3]
    f_i_minus_1   = fn_s[-4]
    f_i_minus_2   = fn_s[-5]

    monitor = Monitor()
    this_interp = HB(
            x_i_minus_2, x_i_minus_1, x_i_minus_0_5, x_i, x_i_plus_1,
            y_i_minus_2, f_i_minus_2,
            y_i_minus_1, f_i_minus_1,
            y_i_minus_0_5, f_i_minus_0_5,
            y_i, f_i,
            y_i_plus_1, f_i_plus_1,
            monitor
    )
    interps.append(this_interp)

    while xn < xend:
        (k, yn_plus_1, yn_plus_1_higher_order) = one_step(fun, xn, yn, f_start, h)

        x_i, y_i = res[-1]
        x_i_minus_1, y_i_minus_1 = res[-2]
        x_i_minus_2, y_i_minus_2 = res[-3]
        x_i_plus_1, y_i_plus_1 = x_i + h, yn_plus_1_higher_order

        f_i = fn_s[-1]
        f_i_minus_1 = fn_s[-2]
        f_i_minus_2 = fn_s[-3]
        f_i_plus_1 = fun(x_i_plus_1, y_i_plus_1)[0]

        # to eliminate one parameter, we define an x value between x_i_minus_1 and x_i
        # this allows us to have a more resilient interpolant
        # we are forced to make a function evaluation so that the data is of the correct order and not affected by interpolant order
        prev_interp = interps[-1]
        x_i_minus_0_5 = x_i - ((x_i - x_i_minus_1) / 2)
        y_i_minus_0_5 = solution([ x_i_minus_0_5 ])[0]
        y_i_minus_0_5_eval = prev_interp.eval(x_i_minus_0_5)
        print("difference between sol and eval", abs(y_i_minus_0_5 - y_i_minus_0_5_eval))
        f_i_minus_0_5 = fun(x_i_minus_0_5, y_i_minus_0_5)[0]

        this_interp = HB(
            x_i_minus_2, x_i_minus_1, x_i_minus_0_5, x_i, x_i_plus_1,
            y_i_minus_2, f_i_minus_2,
            y_i_minus_1, f_i_minus_1,
            # BUG => need to use the actual solution y_i_minus_0_5 instead of the eval
            # y_i_minus_0_5 is too sensitive.... even an error between the eval and the actual solution of 1e-9 
            # makes an interpolant that is bad....
            y_i_minus_0_5_eval, f_i_minus_0_5,
            # y_i_minus_0_5, f_i_minus_0_5,
            y_i, f_i,
            y_i_plus_1, f_i_plus_1,
            monitor
        )

        monitor.n_steps += 1

        """
        # we test the interpolant to check if the Hermite Birkhoff conditions are met as intended
        if (abs(this_interp.eval(x_i_minus_2) - y_i_minus_2))  > 1e-9: print("wrong y_i_minus_2", abs(this_interp.eval(x_i_minus_2) - y_i_minus_2))
        if (abs(this_interp.eval(x_i_minus_1) - y_i_minus_1))  > 1e-9: print("wrong y_i_minus_1", abs(this_interp.eval(x_i_minus_1) - y_i_minus_1))
        if (abs(this_interp.eval(x_i_minus_0_5) - y_i_minus_0_5))  > 1e-9: print("wrong y_i_minus_0_5", abs(this_interp.eval(x_i_minus_0_5) - y_i_minus_0_5))
        if (abs(this_interp.eval(x_i)         - y_i))          > 1e-9: print("wrong y_i",         abs(this_interp.eval(x_i)         - y_i))
        if (abs(this_interp.eval(x_i_plus_1)  - y_i_plus_1))   > 1e-9: print("wrong y_i_plus_1",  abs(this_interp.eval(x_i_plus_1)  - y_i_plus_1))

        if (abs(this_interp.prime(x_i_minus_2) - f_i_minus_2)) > 1e-9: print("wrong f_i_minus_2", abs(this_interp.prime(x_i_minus_2) - f_i_minus_2))
        if (abs(this_interp.prime(x_i_minus_1) - f_i_minus_1)) > 1e-9: print("wrong f_i_minus_1", abs(this_interp.prime(x_i_minus_1) - f_i_minus_1))
        if (abs(this_interp.prime(x_i_minus_0_5) - f_i_minus_0_5)) > 1e-9: print("wrong f_i_minus_1", abs(this_interp.prime(x_i_minus_0_5) - f_i_minus_0_5))
        if (abs(this_interp.prime(x_i)         - f_i))         > 1e-9: print("wrong f_i",         abs(this_interp.prime(x_i)         - f_i))
        if (abs(this_interp.prime(x_i_plus_1)  - f_i_plus_1))  > 1e-9: print("wrong f_i_plus_1",  abs(this_interp.prime(x_i_plus_1)  - f_i_plus_1))
        """        

        # defect control on [x_i to x_i_plus_1]
        h_i = x_i_plus_1 - x_i
        x_sample_1 = x_i + 0.4 * h_i
        defect_sample_1 = abs( 
            this_interp.prime(x_sample_1) - fun( x_sample_1, this_interp.eval(x_sample_1) )[0] 
        )

        x_sample_2 = x_i + 0.8 * h_i
        defect_sample_2 = abs(
            this_interp.prime(x_sample_2) - fun( x_sample_2, this_interp.eval(x_sample_2) )[0]
        )
        max_defect = max(defect_sample_1, defect_sample_2)

        if max_defect < tol:
            # accept the step, by moving the x and the y
            xn = x_i_plus_1
            yn = y_i_plus_1

            # push the intermediate value back into the res array
            res.pop() # pops (x_i, y_i)
            res.append( (x_i_minus_0_5, y_i_minus_0_5) )
            res.append( (x_i, y_i) )

            # push the new step into the array
            res.append( (xn, yn) )

            f_start = f_i_plus_1

            # push the intermediate function evaluations into the res array
            fn_s.pop()
            fn_s.append( f_i_minus_0_5 )
            fn_s.append( f_i )

            # push the function evaluation at x_i_plus_1
            fn_s.append(f_start)

            monitor.n_successful_steps += 1
            print("================== accept step", x_i, "with step-size", h)
            interps.append(this_interp)
            if max_defect < (tol / 10):
                h *= 2
        else:
            h /= 2
            print("failed step at", x_i, "new_h", h)

    print("tolerance=", tol)
    monitor.print()
    print("================================\n")
    continuous_sol = ContinuousSolution()
    continuous_sol.extend(interps)
    return (
        res, 
        continuous_sol.eval,
        continuous_sol.prime,
        create_defect_samplings(res, fn_s, monitor)
    )

# =================================================================================
# the following attempt is when the solver is to keep alpha at 1 throughout the integration

# will also have solution for the first step as a proof of concept
def rk_defect_control_static_alpha_beta_func_call(fun, t_span, y0, tol, solution):
    xn, xend = t_span
    yn = y0
    f_start = fun(xn, yn)[0] 
    
    res = [ (xn, yn) ]
    fn_s = [f_start]

    h = 1e-1 # as HB10 V is at 1e-1 sqrt(tol)

    for _ in range(4): # take 4 steps to build the first interpolant
        xn = xn + h
        yn = solution([ xn ])[0]
        f_start = fun(xn, yn)[0]
        res.append( (xn, yn) )
        fn_s.append(f_start)

    x_i_plus_1, y_i_plus_1       = res[-1]
    x_i, y_i                     = res[-2]
    x_i_minus_0_5, y_i_minus_0_5 = res[-3]
    x_i_minus_1, y_i_minus_1     = res[-4]
    x_i_minus_2, y_i_minus_2     = res[-5]

    f_i_plus_1    = fn_s[-1]
    f_i           = fn_s[-2]
    f_i_minus_0_5 = fn_s[-3]
    f_i_minus_1   = fn_s[-4]
    f_i_minus_2   = fn_s[-5]

    monitor = Monitor()
    this_interp = HB(
        x_i_minus_2, x_i_minus_1, x_i_minus_0_5, x_i, x_i_plus_1,
        y_i_minus_2, f_i_minus_2,
        y_i_minus_1, f_i_minus_1,
        y_i_minus_0_5, f_i_minus_0_5,
        y_i, f_i,
        y_i_plus_1, f_i_plus_1,
        monitor
    )
    continous_sol = ContinuousSolution()
    continous_sol.append(this_interp)

    while xn < xend:
        (k, yn_plus_1, yn_plus_1_higher_order) = one_step(fun, xn, yn, f_start, h)

        x_i, y_i = res[-1]
        x_i_plus_1, y_i_plus_1 = x_i + h, yn_plus_1_higher_order

        f_i = fn_s[-1]
        f_i_plus_1 = fun(x_i_plus_1, y_i_plus_1)[0]

        x_i_minus_1 = x_i - h
        y_i_minus_1 = continous_sol.eval(x_i_minus_1)
        f_i_minus_1 = fun(x_i_minus_1, y_i_minus_1)[0]

        x_i_minus_2 = x_i_minus_1 - h
        y_i_minus_2 = continous_sol.eval(x_i_minus_2)
        f_i_minus_2 = fun(x_i_minus_2, y_i_minus_2)[0]

        # we define a point BETWEEN x_i_minus_1 to x_i so that we can use one less param
        x_i_minus_0_5 = (x_i + x_i_minus_1) / 2
        y_i_minus_0_5 = continous_sol.eval(x_i_minus_0_5)
        f_i_minus_0_5 = fun(x_i_minus_0_5, y_i_minus_0_5)[0]

        this_interp = HB(
            x_i_minus_2, x_i_minus_1, x_i_minus_0_5, x_i, x_i_plus_1,
            y_i_minus_2, f_i_minus_2,
            y_i_minus_1, f_i_minus_1,
            y_i_minus_0_5, f_i_minus_0_5,
            y_i, f_i,
            y_i_plus_1, f_i_plus_1,
            monitor
        )

        """
        # we test the interpolant to check if the Hermite Birkhoff conditions are met as intended
        if (abs(this_interp.eval(x_i_minus_1) - y_i_minus_1))  > 1e-12: print("wrong y_i_minus_1", abs(this_interp.eval(x_i_minus_1) - y_i_minus_1))
        if (abs(this_interp.eval(x_i)         - y_i))          > 1e-12: print("wrong y_i",         abs(this_interp.eval(x_i)         - y_i))
        if (abs(this_interp.eval(x_i_plus_1)  - y_i_plus_1))   > 1e-12: print("wrong y_i_plus_1",  abs(this_interp.eval(x_i_plus_1)  - y_i_plus_1))

        if (abs(this_interp.prime(x_i_minus_1) - f_i_minus_1)) > 1e-12: print("wrong f_i_minus_1", abs(this_interp.prime(x_i_minus_1) - f_i_minus_1))
        if (abs(this_interp.prime(x_i)         - f_i))         > 1e-12: print("wrong f_i",         abs(this_interp.prime(x_i)         - f_i))
        if (abs(this_interp.prime(x_i_plus_1)  - f_i_plus_1))  > 1e-12: print("wrong f_i_plus_1",  abs(this_interp.prime(x_i_plus_1)  - f_i_plus_1))
        """        

        # defect control on [x_i to x_i_plus_1]
        h_i = x_i_plus_1 - x_i
        x_sample_1 = x_i + 0.4 * h_i
        defect_sample_1 = abs( 
            this_interp.prime(x_sample_1) - fun( x_sample_1, this_interp.eval(x_sample_1) )[0] 
        )

        x_sample_2 = x_i + 0.8 * h_i
        defect_sample_2 = abs(
            this_interp.prime(x_sample_2) - fun( x_sample_2, this_interp.eval(x_sample_2) )[0]
        )
        max_defect = max(defect_sample_1, defect_sample_2)

        monitor.n_steps += 1

        # print("max_defect", max_defect)
        if True: # max_defect < tol:
            # accept the step, by moving the x and the y
            xn = x_i_plus_1
            yn = y_i_plus_1
            res.append( (xn, yn) )

            f_start = f_i_plus_1
            fn_s.append(f_start)

            continous_sol.append(this_interp)

            monitor.n_successful_steps += 1

            if max_defect < (tol / 10):
                pass
                # h *= 2
        else:
            # print("tolerance not satisfied")
            h /= 2

    print("tolerance=", tol)
    monitor.print()
    print("================================\n")

    return (
        res, 
        continous_sol.eval,
        continous_sol.prime,
        create_defect_samplings(res, fn_s, monitor)
    )

# fixed step size solver to debug interpolant
def rk_fixed_step(fun, t_span, y0, nsteps=100):
    # theory:
        # we take step by step from the start in t_span[0] to the end in t_span[1]
        # each step is a rk_step as per one_step(), the size of the step, h, the step-size, is fixed based on the nubmer of steps the user wants to take
    xn = t_span[0]
    xend = t_span[1]

    yn = y0
    f_start = fun(xn, yn)[0] # each time we call, the function 'fun', we will have to extract the first value as solve_ivp wants vectorised model functions
    
    res = [(xn, yn)]
    fn_s = [f_start]

    h = (xend - xn) / nsteps
    while xn < xend:
        (k, yn_plus_1, yn_plus_1_higher_order) = one_step(fun, xn, yn, f_start, h)

        # error = abs(yn_plus_1_higher_order - yn_plus_1)
        
        # accept the step, by moving the x and the y
        xn = xn + h
        yn = yn_plus_1_higher_order
        res.append( (xn, yn) )
        
        f_start = fun(xn, yn)[0] # we make a final function evalution at the current step
        fn_s.append(f_start)

    monitor = Monitor()
    continous_sol = create_continuous_sol_from_results(res, fn_s, monitor)
    return (
        res, 
        continous_sol.eval,
        continous_sol.prime,
        create_defect_samplings(res, fn_s, monitor)
    )
