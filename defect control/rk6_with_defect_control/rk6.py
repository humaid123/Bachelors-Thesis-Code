from math import sqrt
from HB import HB

# http://people.math.sfu.ca/~jverner/RKV65.IIIXb.Efficient.00000144617.081204.RATOnWeb
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

    yn_plus_1_higher_order = yn + h * sigma_prod(B_HAT, k, 0, n_stages) 

    return (k, yn_plus_1, yn_plus_1_higher_order)

# ===============================================================================================================
# Helper functions

def create_continuous_sol_from_interpolants(interps):
    def sol(x) -> float:
        for hb in interps:
            # if (hb.x_i_minus_1 <= x <= hb.x_i):
            if (hb.x_i <= x <= hb.x_i_plus_1):
                return hb.eval(x)
        # last_hb = interps[-1]
        # if (last_hb.x_i <= x <= last_hb.x_i_plus_1):
        #     return last_hb.eval(x)
        
        first_hb = interps[0]
        if (first_hb.x_i_minus_1 <= x <= first_hb.x_i):
            return first_hb.eval(x)
        print(f"ERROR: {x} is outside of the solution range: {interps[0].x_i_minus_1} <= x <= {interps[-1].x_i_plus_1}")
        return -1

    return sol


def create_continuous_sol_from_results(res, fn_s):
    interps = []
    
    for i in range(len(res) - 2):
        x_i_minus_1, y_i_minus_1 = res[i]    
        x_i, y_i                 = res[i + 1]    
        x_i_plus_1, y_i_plus_1   = res[i + 2]

        f_i_minus_1 = fn_s[i]    
        f_i         = fn_s[i + 1]    
        f_i_plus_1  = fn_s[i + 2]
        
        interps.append(
            HB (
                x_i_minus_1, x_i, x_i_plus_1,
                y_i_minus_1, f_i_minus_1,
                y_i, f_i,
                y_i_plus_1, f_i_plus_1 
            )
        )
    return create_continuous_sol_from_interpolants(interps)

def create_continuous_first_derivatives_from_interpolants(interps):
    def sol(x) -> float:
        for hb in interps:
            # if (hb.x_i_minus_1 <= x <= hb.x_i):
            if (hb.x_i <= x <= hb.x_i_plus_1):
                return hb.prime(x)
        # last_hb = interps[-1]
        # if (last_hb.x_i <= x <= last_hb.x_i_plus_1):
        #     return last_hb.prime(x)

        first_hb = interps[0]
        if (first_hb.x_i_minus_1 <= x <= first_hb.x_i):
            return first_hb.prime(x)
        
        print(f"ERROR: {x} is outside of the solution range: {interps[0].x_i_minus_1} <= x <= {interps[-1].x_i_plus_1}")
        return -1

    return sol


def create_continuous_first_derivatives_from_results(res, fn_s):
    interps = []
    
    for i in range(len(res) - 2):
        x_i_minus_1, y_i_minus_1 = res[i]    
        x_i, y_i                 = res[i + 1]    
        x_i_plus_1, y_i_plus_1   = res[i + 2]

        f_i_minus_1 = fn_s[i]    
        f_i         = fn_s[i + 1]    
        f_i_plus_1  = fn_s[i + 2]
        
        interps.append(
            HB (
                x_i_minus_1, x_i, x_i_plus_1,
                y_i_minus_1, f_i_minus_1,
                y_i, f_i,
                y_i_plus_1, f_i_plus_1 
            )
        )
    return create_continuous_first_derivatives_from_interpolants(interps)

def create_defect_samplings(res, fn_s):
    result = []
    for i in range(len(res) - 2):
        x_i_minus_1, y_i_minus_1 = res[i]    
        x_i, y_i                 = res[i + 1]    
        x_i_plus_1, y_i_plus_1   = res[i + 2]

        f_i_minus_1 = fn_s[i]    
        f_i         = fn_s[i + 1]    
        f_i_plus_1  = fn_s[i + 2]
        
        interp = HB (
                x_i_minus_1, x_i, x_i_plus_1,
                y_i_minus_1, f_i_minus_1,
                y_i, f_i,
                y_i_plus_1, f_i_plus_1 
        )
        result.append( (x_i_minus_1, x_i, x_i_plus_1, interp) )
    return result
        

# =================================================================================================================
# start of fixed step-size solver

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


    return (
        res, 
        create_continuous_sol_from_results(res, fn_s),
        create_continuous_first_derivatives_from_results(res, fn_s),
        create_defect_samplings(res, fn_s)
    )

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

    h = sqrt(tol)
    while xn < xend:
        (k, yn_plus_1, yn_plus_1_higher_order) = one_step(fun, xn, yn, f_start, h)

        error = abs(yn_plus_1_higher_order - yn_plus_1)

        if error < tol:
            # accept the step, by moving the x and the y
            xn = xn + h
            yn = yn_plus_1
            res.append( (xn, yn) )

            # so for now, I do an additional function evaluation
            f_start = fun(xn, yn)[0]
            fn_s.append(f_start)

            if error < (tol / 10):
                h *= 2
        else:
            h /= 2
    
    return (
        res, 
        create_continuous_sol_from_results(res, fn_s),
        create_continuous_first_derivatives_from_results(res, fn_s),
        create_defect_samplings(res, fn_s)
    )


# =============================================================================================================================
# start of defect control solver

def first_step(fun, xn, yn, f_start, tol):
    """
    For the first step, we take a step using a stronger tolerance.
    We use this stronger tolerance to also decide on the initial h value.
    """
    h = sqrt(tol)
    stricter_tol = tol / 100
    error = float('inf')
    while error > stricter_tol:
        (k, yn_plus_1, yn_plus_1_higher_order) = one_step(fun, xn, yn, f_start, h)

        # we do error control on the first few steps until we get enough function evaluations to be able to do defect control
        error = abs(yn_plus_1_higher_order - yn_plus_1)
        
        # we control the error at the end of the step with a stricter tolerance for the first step
        if error < stricter_tol:
            # accept the step, by moving the x and the y
            xn = xn + h
            yn = yn_plus_1_higher_order
            h *= 2
        else:
            h /= 2

    f_start = fun(xn, yn)[0]

    # print("first step", h, sqrt(tol))
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
    interps = []

    # we do strict error control on the first step
    # we update all the values of xn, yn, fstart based on what happened on the first step
    (xn, yn, f_start, h) = first_step(fun, xn, yn, f_start, tol)
    res.append( (xn, yn) )
    fn_s.append(f_start)

    n_steps = 0
    n_successful_steps = 0

    while xn < xend:
        (k, yn_plus_1, yn_plus_1_higher_order) = one_step(fun, xn, yn, f_start, h)

        x_i, y_i = res[-1]
        x_i_minus_1, y_i_minus_1 = res[-2]
        x_i_plus_1, y_i_plus_1 = x_i + h, yn_plus_1_higher_order

        f_i = fn_s[-1]
        f_i_minus_1 = fn_s[-2]
        # BUG => this is x_i, y_i NOT x_i_plus_1, y_i_plus_1 => f_i_plus_1 = fun(xn, yn)[0]
        f_i_plus_1 = fun(x_i_plus_1, y_i_plus_1)[0]

        this_interp = HB(
            x_i_minus_1, x_i, x_i_plus_1,
            y_i_minus_1, f_i_minus_1,
            y_i, f_i,
            y_i_plus_1, f_i_plus_1 
        )

        n_steps += 1

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

        if max_defect < tol:
            # accept the step, by moving the x and the y
            xn = x_i_plus_1
            yn = y_i_plus_1
            res.append( (xn, yn) )

            f_start = f_i_plus_1
            fn_s.append(f_start)

            n_successful_steps += 1

            interps.append(this_interp)
            if max_defect < (tol / 10):
                h *= 2
        else:
            h /= 2

    print("tolerance=", tol)
    print("number of successful steps=", n_successful_steps)
    print("number of steps=", n_steps)
    print("================================\n")
    return (
        res, 
        create_continuous_sol_from_interpolants(interps),
        create_continuous_first_derivatives_from_interpolants(interps),
        create_defect_samplings(res, fn_s)
    )


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
    h = sqrt(tol)
    xn = xn + h
    yn = solution([xn])[0]
    f_start = fun(xn, yn)[0]
    res.append( (xn, yn) )
    fn_s.append(f_start)

    n_steps = 0
    n_successful_steps = 0

    while xn < xend:
        (k, yn_plus_1, yn_plus_1_higher_order) = one_step(fun, xn, yn, f_start, h)

        x_i, y_i = res[-1]
        x_i_minus_1, y_i_minus_1 = res[-2]
        x_i_plus_1, y_i_plus_1 = x_i + h, yn_plus_1_higher_order

        f_i = fn_s[-1]
        f_i_minus_1 = fn_s[-2]
        f_i_plus_1 = fun(x_i_plus_1, y_i_plus_1)[0]

        this_interp = HB(
            x_i_minus_1, x_i, x_i_plus_1,
            y_i_minus_1, f_i_minus_1,
            y_i, f_i,
            y_i_plus_1, f_i_plus_1 
        )

        n_steps += 1

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

        if max_defect < tol:
            # accept the step, by moving the x and the y
            xn = x_i_plus_1
            yn = y_i_plus_1
            res.append( (xn, yn) )

            f_start = f_i_plus_1
            fn_s.append(f_start)

            n_successful_steps += 1

            interps.append(this_interp)
            if max_defect < (tol / 10):
                h *= 2
        else:
            h /= 2

    print("tolerance=", tol)
    print("number of successful steps=", n_successful_steps)
    print("number of steps=", n_steps)
    print("================================\n")
    return (
        res, 
        create_continuous_sol_from_interpolants(interps),
        create_continuous_first_derivatives_from_interpolants(interps),
        create_defect_samplings(res, fn_s)
    )

# =================================================================================
# the following attempt is when the solver is to keep alpha at 1 throughout the integration
class ContinuousSolution:
    def __init__(self) -> None:
        self.interps = []
    
    def eval(self, x) -> float:
        for hb in self.interps:
            if (hb.x_i <= x <= hb.x_i_plus_1):
                return hb.eval(x)

        first_hb = self.interps[0]
        if (first_hb.x_i_minus_1 <= x <= first_hb.x_i):
            return first_hb.eval(x)
        print(f"ERROR: {x} is outside of the solution range: {first_hb.x_i_minus_1} <= x <= {self.interps[-1].x_i_plus_1}")
        return -1

    def prime(self, x) -> float:
        for hb in self.interps:
            if (hb.x_i <= x <= hb.x_i_plus_1):
                return hb.prime(x)

        first_hb = self.interps[0]
        if (first_hb.x_i_minus_1 <= x <= first_hb.x_i):
            return first_hb.prime(x)
        
        print(f"ERROR: {x} is outside of the solution range: {first_hb.x_i_minus_1} <= x <= {self.interps[-1].x_i_plus_1}")
        return -1
    
    def append(self, interp) -> None:
        self.interps.append(interp)

# will also have solution for the first step as a proof of concept
def rk_defect_control_static_alpha(fun, t_span, y0, tol, solution):
    xn, xend = t_span
    yn = y0
    f_start = fun(xn, yn)[0] 
    
    res = [ (xn, yn) ]
    fn_s = [f_start]

    # first solution
    h = sqrt(tol)
    xn = xn + h
    yn = solution([ xn ])[0]
    f_start = fun(xn, yn)[0]
    res.append( (xn, yn) )
    fn_s.append(f_start)

    xn = xn + h
    yn = solution([ xn ])[0]
    f_start = fun(xn, yn)[0]
    res.append( (xn, yn) )
    fn_s.append(f_start)

    x_i_plus_1, y_i_plus_1   = res[-1]
    x_i, y_i                 = res[-2]
    x_i_minus_1, y_i_minus_1 = res[-3]
    f_i_plus_1  = fn_s[-1]
    f_i         = fn_s[-2]
    f_i_minus_1 = fn_s[-3]

    this_interp = HB(
        x_i_minus_1, x_i, x_i_plus_1,
        y_i_minus_1, f_i_minus_1,
        y_i, f_i,
        y_i_plus_1, f_i_plus_1 
    )

    continous_sol = ContinuousSolution()
    continous_sol.append(this_interp)

    n_steps = 0
    n_successful_steps = 0

    index = 0
    while xn < xend:
        (k, yn_plus_1, yn_plus_1_higher_order) = one_step(fun, xn, yn, f_start, h)

        x_i, y_i = res[-1]
        x_i_plus_1, y_i_plus_1 = x_i + h, yn_plus_1_higher_order

        f_i = fn_s[-1]
        f_i_plus_1 = fun(x_i_plus_1, y_i_plus_1)[0]

        x_i_minus_1 = x_i - h

        y_i_minus_1 = continous_sol.eval(x_i_minus_1)
        f_i_minus_1 = continous_sol.prime(x_i_minus_1)

        this_interp = HB(
            x_i_minus_1, x_i, x_i_plus_1,
            y_i_minus_1, f_i_minus_1,
            y_i, f_i,
            y_i_plus_1, f_i_plus_1 
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

        n_steps += 1

        # print("max_defect", max_defect)
        if max_defect < tol:
            # accept the step, by moving the x and the y
            xn = x_i_plus_1
            yn = y_i_plus_1
            res.append( (xn, yn) )

            f_start = f_i_plus_1
            fn_s.append(f_start)

            continous_sol.append(this_interp)

            index += 1
            n_successful_steps += 1

            if max_defect < (tol / 10):
                h *= 2
        else:
            # print("tolerance not satisfied")
            h /= 2

    print("tolerance=", tol)
    print("number of successful steps=", n_successful_steps)
    print("number of steps=", n_steps)
    print("================================\n")

    return (
        res, 
        continous_sol.eval,
        continous_sol.prime,
        create_defect_samplings(res, fn_s)
    )