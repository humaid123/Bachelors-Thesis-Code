from math import sqrt
from CRK8Interp import CRK8ContinuousSolution, CRK8Interp
from CRK7Interp import CRK7ContinuousSolution, CRK7Interp

# http://people.math.sfu.ca/~jverner/RKV65.IIIXb.Efficient.00000144617.081204.RATOnWeb
A = [
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
    [0.05, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
    [-0.0069931640625, 0.1135556640625, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
    [0.0399609375, 0.0, 0.1198828125, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
    [0.36139756280045754, 0.0, -1.3415240667004928, 1.3701265039000352, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
    [0.049047202797202795, 0.0, 0.0, 0.23509720422144048, 0.18085559298135673, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
    [0.06169289044289044, 0.0, 0.0, 0.11236568314640277, -0.03885046071451367, 0.01979188712522046, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
    [-1.767630240222327, 0.0, 0.0, -62.5, -6.061889377376669, 5.6508231982227635, 65.62169641937624, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
    [-1.1809450665549708, 0.0, 0.0, -41.50473441114321, -4.434438319103725, 4.260408188586133, 43.75364022446172, 0.00787142548991231, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
    [-1.2814059994414884, 0.0, 0.0, -45.047139960139866, -4.731362069449577, 4.514967016593808, 47.44909557172985, 0.010592282971116612, -0.0057468422638446166, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
    [-1.7244701342624853, 0.0, 0.0, -60.92349008483054, -5.951518376222393, 5.556523730698456, 63.98301198033305, 0.014642028250414961, 0.06460408772358203, -0.0793032316900888, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
    [-3.301622667747079, 0.0, 0.0, -118.01127235975251, -10.141422388456112, 9.139311332232058, 123.37594282840426, 4.62324437887458, -3.3832777380682018, 4.527592100324618, -5.828495485811623, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
    [-3.039515033766309, 0.0, 0.0, -109.26086808941763, -9.290642497400293, 8.43050498176491, 114.20100103783314, -0.9637271342145479, -5.0348840888021895, 5.958130824002923, 0.0, 0.0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
    [0.04427989419007951, 0.0, 0.0, 0.0, 0.0, 0.3541049391724449, 0.2479692154956438, -15.694202038838084, 25.084064965558564, -31.738367786260277, 22.938283273988784, -0.2361324633071542, 0.0, 0, 0, 0, 0, 0, 0, 0, 0], 
    [0.04620700646754963, 0.0, 0.0, 0.0, 0.0, 0.045039041608424805, 0.23368166977134244, 37.83901368421068, -15.949113289454246, 23.028368351816102, -44.85578507769412, -0.06379858768647444, 0.0, -0.012595035543861663, 0, 0, 0, 0, 0, 0, 0], 
    [0.05037946855482041, 0.0, 0.0, 0.0, 0.0, 0.041098361310460796, 0.17180541533481958, 4.614105319981519, -1.7916678830853965, 2.531658930485041, -5.324977860205731, -0.03065532595385635, 0.0, -0.005254479979429613, -0.08399194644224793, 0, 0, 0, 0, 0, 0], 
    [0.0408289713299708, 0.0, 0.0, 0.0, 0.0, 0.4244479514247632, 0.23260915312752345, 2.677982520711806, 0.7420826657338945, 0.1460377847941461, -3.579344509890565, 0.11388443896001738, 0.0, 0.012677906510331901, -0.07443436349946675, 0.047827480797578516, 0, 0, 0, 0, 0], 
    [0.052126823936684136, 0.0, 0.0, 0.0, 0.0, 0.053925083967447975, 0.01660758097434641, -4.45448575792678, 6.835218278632146, -8.711334822181994, 6.491635839232917, -0.07072551809844346, 0.0, -0.018540314919932164, 0.023504021054353848, 0.2344795103407822, -0.08241072501152899, 0, 0, 0, 0], 
    [0.05020102870355714, 0.0, 0.0, 0.0, 0.0, 0.1552209034795498, 0.1264268424089235, -5.149206303539847, 8.46834099903693, -10.662130681081495, 7.541833224959729, -0.07436968113832143, 0.0, -0.020558876866183826, 0.07753795264710298, 0.10462592203525443, -0.11792133064519794, 0.0, 0, 0, 0], 
    [0.03737341446457826, 0.0, 0.0, 0.0, 0.0, 0.35049307053383166, 0.49226528193730257, 8.553695439359313, -10.353172990305913, 13.83320427252915, -12.280924330784618, 0.17191515956565098, 0.0, 0.036415831143144964, 0.02961920580288763, -0.2651793938627067, 0.09429503961738067, 0.0, 0.0, 0, 0], 
    [0.039390583455282506, 0.0, 0.0, 0.0, 0.0, 0.3558516141234424, 0.419738222595261, 0.8720449778071941, 0.8989520834876595, -0.6305806161059884, -1.1218872205954835, 0.04298219512400197, 0.0, 0.013325575668739157, 0.018762270539641482, -0.18594111329221055, 0.17736142719246029, 0.0, 0.0, 0.0, 0]
]
C = [0.0, 0.05, 0.1065625, 0.15984375, 0.39, 0.465, 0.155, 0.943, 0.901802041735857, 0.909, 0.94, 1.0, 1.0, 1.0, 0.3110177634953864, 0.1725, 0.7846, 0.37, 0.5, 0.7, 0.9]
B = [0.04427989419007951, 0.0, 0.0, 0.0, 0.0, 0.3541049391724449, 0.2479692154956438, -15.694202038838084, 25.084064965558564, -31.738367786260277, 22.938283273988784, -0.2361324633071542, 0.0, 0, 0, 0, 0, 0, 0, 0, 0]
B_HAT = [0.044312615229089795, 0.0, 0.0, 0.0, 0.0, 0.35460956423432266, 0.2478480431366653, 4.4481347324757845, 19.846886366118735, -23.58162337746562, 0.0, 0.0, 
-0.36016794372897754, 0, 0, 0, 0, 0, 0, 0, 0]
n_stages =  21

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

    continous_sol = CRK8ContinuousSolution()
    lower_continous_sol = CRK7ContinuousSolution()

    nsteps = 0
    n_successful_steps = 0

    h = sqrt(tol)
    while xn < xend:
        (k, yn_plus_1, yn_plus_1_higher_order) = one_step(fun, xn, yn, f_start, h)

        x_i = xn
        x_i_plus_1 = xn + h
        y_i = yn
        
        crk8 = CRK8Interp(k, y_i, x_i, x_i_plus_1)
        crk7 = CRK7Interp(k, y_i, x_i, x_i_plus_1)

        h_i = x_i_plus_1 - x_i
        x_sample_1 = x_i + 0.5 * h_i
        error_sample1 = abs(
            crk7.eval(x_sample_1) - crk8.eval(x_sample_1)
        )

        x_sample_2 = x_i + 0.8 * h_i
        error_sample2 = abs(
            crk7.eval(x_sample_2) - crk8.eval(x_sample_2)
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
            continous_sol.append(crk8)
            lower_continous_sol.append(crk7)     
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
