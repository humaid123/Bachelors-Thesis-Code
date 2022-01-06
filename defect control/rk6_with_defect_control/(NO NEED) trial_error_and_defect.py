
THIS DOES NOT IMPROVE THE SITUATION BY A LOT...

# ============================================================================================================
# here we have a solver that does defect control ONLY after the error has been satisfied
# we thus sometimes save some function evaluations when doing the defect control

def rk_error_and_defect_control(fun, t_span, y0, tol):
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

    # from my experiments we can see that the solver is too strict
    # we relax the tolerance a little bit to allow faster execution
    # tol *= 5
    

    while xn < xend:
        (k, yn_plus_1, yn_plus_1_higher_order) = one_step(fun, xn, yn, f_start, h)

        error = abs(yn_plus_1_higher_order - yn_plus_1)
        if (error > tol):
            print("error not satisfied, so we reject and don't even do dfect control")
            h /= 2
            continue

        x_i, y_i = res[-1]
        x_i_minus_1, y_i_minus_1 = res[-2]
        x_i_plus_1, y_i_plus_1 = x_i + h, yn_plus_1_higher_order

        f_i = fn_s[-1]
        f_i_minus_1 = fn_s[-2]
        f_i_plus_1 = fun(xn, yn)[0]
        # print("test k[-1]", k[-1] - f_i_plus_1)

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
        
        # print("h", h, initial_h)
        # print("max_defect", max_defect, max_defect < tol)

        # print( "num succesful steps taken so far = ", len(interps) + 1 )
        # print( "x_i_minus_1", x_i_minus_1, "h",h)

        if max_defect < tol:
            # accept the step, by moving the x and the y
            xn = x_i_plus_1
            yn = y_i_plus_1
            res.append( (xn, yn) )

            f_start = f_i_plus_1
            fn_s.append(f_start)


            interps.append(this_interp)
            if max_defect < (tol / 1.5):
                h *= 3
            
            """
            The code below seems to WORSEN the situation not improve it
            # THROUGHOUT MY EXPERIMENTS, we have seen that an alpha closer to 1
            # and a h value close to 1e-2 lead to better result.
            # I try to force it to take bigger steps if it is taking similar sized steps
            if len(interps) >= 2:
                if (interps[-2].alpha == 1) and (interps[-1].alpha == 1) and (h < 1e-2):
                    # from what we know from the perfect convergence experiment
                    # when alpha is 1, the v-shaped graph's v is around 1e-2, 1e-3
                    # so there is a chance that if we keep decreasing, the step => the solver starts getting INNACURATE instead
                    # of more accurate...
                    # print("making it take a bigger step")
                    h *= 3
            """
            
        else:
            h /= 2

    return (
        res, 
        create_continuous_sol_from_interpolants(interps),
        create_continuous_first_derivatives_from_interpolants(interps),
        create_defect_samplings(res, fn_s)
    )

