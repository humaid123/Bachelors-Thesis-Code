from math import sqrt
from math import sqrt
from HB8_second_scheme import HB8, HB8ContinuousSolution
from scipy.integrate import fixed_quad
from HB10_fourth_scheme import HB10, HB10ContinuousSolution


# http://people.math.sfu.ca/~jverner/RKV65.IIIXb.Efficient.00000144617.081204.RATOnWeb
# A = [
#     [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
#     [0.0556, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
#     [0.007953676685071909, 0.09462409627853294, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
#     [0.03846666486135181, 0.0, 0.11539999458405548, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
#     [0.3843917952499957, 0.0, -1.4413754967533892, 1.4415837015033934, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
#     [0.04616799272252459, 0.0, 0.0, 0.23076667147858013, 0.1845653357988953, 0, 0, 0, 0, 0, 0, 0, 0], 
#     [0.05983406569816849, 0.0, 0.0, 0.11107098836580696, -0.034214310915191934, 0.017109256851216476, 0, 0, 0, 0, 0, 0, 0], 
#     [-0.5379500775278769, 0.0, 0.0, -6.937648213098321, -4.662453820973333, 3.995152111599525, 9.000000000000005, 0, 0, 0, 0, 0, 0], 
#     [-1.6324274408037667, 0.0, 0.0, -10.827155649216506, -12.412770216566443, 9.727368979607842, 16.199350914620496, -0.10384430809172447, 0, 0, 0, 0, 0], 
#     [0.4379695061761497, 0.0, 0.0, 3.9395318803002466, 2.860770346811585, -1.7743107088338461, -4.895390517637671, 0.21302485881866956, -0.0593953656351389, 0, 0, 0, 0], 
#     [-1.4741971530123092, 0.0, 0.0, -10.99400456884099, -11.347103595578742, 8.95698732807952, 15.893778872755856, -0.0987525742060384, 0.0048885045849518735, -0.0040968137822505165, 0, 0, 0], 
#     [-2.630059332774936, 0.0, 0.0, -9.174218051381494, -19.181392627716956, 14.64255869375726, 17.529319464433122, -0.37191756017857125, -0.7009961538281585, 0.05101601661228638, 0.8356895510774236, 0, 0], 
#     [0.2157603225684609, 0.0, 0.0, 8.345147326820761, 2.1856623855150725, -1.6872364803128412, -8.7118979011482, 0.024441459344689206, 0.08463787994499963, 0.5434850072670662, 0.0, 0.0, 0]
# ]
# C = [0.0, 0.0556, 0.10257777296360485, 0.15386665944540728, 0.3846,
#      0.4615, 0.1538, 0.8571, 0.9505222795498989, 0.7222, 0.9375, 1.0, 1.0]
# B = [0.04391770364440195, 0.0, 0.0, 0.0, 0.0, 0.35102462530126355, 0.24614282635491153, 0.900324493055835,
#      4.549418727273593, 0.004802501518691216, -4.741054352139637, -0.3545765250090589, 0.0]
# B_HAT = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
# n_stages = 13
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
n_stages =  13 # 21
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
    lower_order_interps = []

    # we do a perfect step for the one_step
    h = 1e-1 # as HB8 Vs at 1e-1 to 1e-2 sqrt(tol)
    xn = xn + h
    yn = solution([xn])[0]
    f_start = fun(xn, yn)[0]
    res.append( (xn, yn) )
    fn_s.append(f_start)

    xn = xn + h
    yn = solution([xn])[0]
    f_start = fun(xn, yn)[0]
    res.append( (xn, yn) )
    fn_s.append(f_start)

    monitor8 = Monitor()
    monitor10 = Monitor()
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

        this_interp_hb8 = HB8(
            x_i_minus_2, x_i_minus_1, x_i, x_i_plus_1,
            y_i_minus_2, f_i_minus_2,
            y_i_minus_1, f_i_minus_1,
            y_i, f_i,
            y_i_plus_1, f_i_plus_1,
            monitor8
        )


        # to eliminate one parameter, we define an x value between x_i_minus_1 and x_i
        # this allows us to have a more resilient interpolant
        # we are forced to make a function evaluation so that the data is of the correct order and not affected by interpolant order
        prev_interp = this_interp_hb8 # interps[-1] if len(interps) >= 1 else this_interp_hb8
        x_i_minus_0_5 = x_i - ((x_i - x_i_minus_1) / 2)
        y_i_minus_0_5 = solution([ x_i_minus_0_5 ])[0]
        y_i_minus_0_5_eval = prev_interp.eval_bary(x_i_minus_0_5)
        # print("difference between sol and eval", abs(y_i_minus_0_5 - y_i_minus_0_5_eval))
        f_i_minus_0_5 = fun(x_i_minus_0_5, y_i_minus_0_5_eval)[0]

        this_interp_hb10 = HB10(
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
            monitor10
        )

        monitor8.n_steps += 1

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
            this_interp_hb10.prime(x_sample_1) - fun( x_sample_1, this_interp_hb10.eval(x_sample_1) )[0] 
        )

        x_sample_2 = x_i + 0.8 * h_i
        defect_sample_2 = abs(
            this_interp_hb10.prime(x_sample_2) - fun( x_sample_2, this_interp_hb10.eval(x_sample_2) )[0]
        )
        max_defect = max(defect_sample_1, defect_sample_2)

        if max_defect < tol:
            # accept the step, by moving the x and the y
            xn = x_i_plus_1
            yn = y_i_plus_1
            res.append( (xn, yn) )

            f_start = f_i_plus_1
            fn_s.append(f_start)

            monitor8.n_successful_steps += 1

            interps.append(this_interp_hb10)
            lower_order_interps.append(this_interp_hb8)

            if max_defect < (tol / 10):
                h *= 2
        else:
            h /= 2

    print("tolerance=", tol)
    print("Monitor8\n===================================")
    monitor8.print()
    print("\n\nMonitor10\n===============================")
    monitor10.print()

    print("================================\n")
    continuous_sol = HB10ContinuousSolution()
    continuous_sol.extend(interps)
    
    lower_continuous_sol = HB8ContinuousSolution()
    lower_continuous_sol.extend(lower_order_interps)
    return (
        res, 
        continuous_sol.eval,
        continuous_sol.prime,
        continuous_sol.create_error_samplings(),
    )

# # =================================================================================
# # the following attempt is when the solver is to keep alpha at 1 throughout the integration

# # will also have solution for the first step as a proof of concept
def rk_defect_control_static_alpha_beta(fun, t_span, y0, tol, solution):
    xn, xend = t_span
    yn = y0
    f_start = fun(xn, yn)[0] # each time we call, the function 'fun', we will have to extract the first value as solve_ivp wants vectorised model functions
    
    res = [(xn, yn)]
    fn_s = [f_start]
    interps = []
    lower_order_interps = []

    # we do a perfect step for the one_step
    h = 1e-1 # as HB8 Vs at 1e-1 to 1e-2 sqrt(tol)
    for _ in range(3):
        xn = xn + h
        yn = solution([xn])[0]
        f_start = fun(xn, yn)[0]
        res.append( (xn, yn) )
        fn_s.append(f_start)

    x_i_plus_1, y_i_plus_1   = res[-1]
    x_i, y_i                 = res[-2]
    x_i_minus_1, y_i_minus_1 = res[-3]
    x_i_minus_2, y_i_minus_2 = res[-4]

    f_i_plus_1  = fn_s[-1]
    f_i         = fn_s[-2]
    f_i_minus_1 = fn_s[-3]
    f_i_minus_2 = fn_s[-4]

    continuous_sol = HB10ContinuousSolution()
    lower_continuous_sol = HB8ContinuousSolution()
    monitor8 = Monitor()
    monitor10 = Monitor()

    this_interp_hb8 = HB8(
        x_i_minus_2, x_i_minus_1, x_i, x_i_plus_1,
        y_i_minus_2, f_i_minus_2,
        y_i_minus_1, f_i_minus_1,
        y_i, f_i,
        y_i_plus_1, f_i_plus_1,
        monitor8
    )

    prev_interp = this_interp_hb8 # interps[-1] if len(interps) >= 1 else this_interp_hb8
    x_i_minus_0_5 = x_i - ((x_i - x_i_minus_1) / 2)
    y_i_minus_0_5 = solution([ x_i_minus_0_5 ])[0]
    y_i_minus_0_5_eval = prev_interp.eval_bary(x_i_minus_0_5)
    # print("difference between sol and eval", abs(y_i_minus_0_5 - y_i_minus_0_5_eval))
    f_i_minus_0_5 = fun(x_i_minus_0_5, y_i_minus_0_5_eval)[0]

    this_interp_hb10 = HB10(
        x_i_minus_2, x_i_minus_1, x_i_minus_0_5, x_i, x_i_plus_1,
        y_i_minus_2, f_i_minus_2,
        y_i_minus_1, f_i_minus_1,
        y_i_minus_0_5_eval, f_i_minus_0_5,
        y_i, f_i,
        y_i_plus_1, f_i_plus_1,
        monitor10
    )
    continuous_sol.append(this_interp_hb10)
    lower_continuous_sol.append(this_interp_hb8)

    while xn < xend:
        (k, yn_plus_1, yn_plus_1_higher_order) = one_step(fun, xn, yn, f_start, h)

        x_i, y_i = res[-1]
        x_i_plus_1, y_i_plus_1 = x_i + h, yn_plus_1_higher_order
        f_i = fn_s[-1]
        f_i_plus_1 = fun(x_i_plus_1, y_i_plus_1)[0]

        # the following is not good enough, we need the hb8...
        # x_i_minus_0_5 = x_i - h # go back 2*h as we break this step for the hb10 to have alpha and beta at 1 
        # y_i_minus_0_5 = continuous_sol.eval(x_i_minus_0_5) # might have to use lower_continuous_sol
        # f_i_minus_0_5 = fun(x_i_minus_0_5, y_i_minus_0_5)[0] # continuous_sol.prime(x_i_minus_1)
        
        x_i_minus_1 = x_i - 2*h # go back 2*h as we break this step for the hb10 to have alpha and beta at 1 
        y_i_minus_1 = continuous_sol.eval(x_i_minus_1) # might have to use lower_continuous_sol
        f_i_minus_1 = fun(x_i_minus_1, y_i_minus_1)[0] # continuous_sol.prime(x_i_minus_1)

        x_i_minus_2 = x_i_minus_1 - h
        y_i_minus_2 = continuous_sol.eval(x_i_minus_2)
        f_i_minus_2 = fun(x_i_minus_2, y_i_minus_2)[0] # continuous_sol.prime(x_i_minus_2)

        this_interp_hb8 = HB8(
            x_i_minus_2, x_i_minus_1, x_i, x_i_plus_1,
            y_i_minus_2, f_i_minus_2,
            y_i_minus_1, f_i_minus_1,
            y_i, f_i,
            y_i_plus_1, f_i_plus_1,
            monitor8
        )


        prev_interp = this_interp_hb8 # interps[-1] if len(interps) >= 1 else this_interp_hb8
        x_i_minus_0_5 = x_i - ((x_i - x_i_minus_1) / 2)
        y_i_minus_0_5 = solution([ x_i_minus_0_5 ])[0]
        y_i_minus_0_5_eval = prev_interp.eval_bary(x_i_minus_0_5)
        # print("difference between sol and eval", abs(y_i_minus_0_5 - y_i_minus_0_5_eval))
        f_i_minus_0_5 = fun(x_i_minus_0_5, y_i_minus_0_5_eval)[0]

        this_interp_hb10 = HB10(
            x_i_minus_2, x_i_minus_1, x_i_minus_0_5, x_i, x_i_plus_1,
            y_i_minus_2, f_i_minus_2,
            y_i_minus_1, f_i_minus_1,
            y_i_minus_0_5_eval, f_i_minus_0_5,
            y_i, f_i,
            y_i_plus_1, f_i_plus_1,
            monitor10
        )

        monitor8.n_steps += 1

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
            this_interp_hb10.prime(x_sample_1) - fun( x_sample_1, this_interp_hb10.eval(x_sample_1) )[0] 
        )

        x_sample_2 = x_i + 0.8 * h_i
        defect_sample_2 = abs(
            this_interp_hb10.prime(x_sample_2) - fun( x_sample_2, this_interp_hb10.eval(x_sample_2) )[0]
        )
        max_defect = max(defect_sample_1, defect_sample_2)

        if max_defect < tol:
            # accept the step, by moving the x and the y
            xn = x_i_plus_1
            yn = y_i_plus_1
            res.append( (xn, yn) )

            f_start = f_i_plus_1
            fn_s.append(f_start)

            monitor8.n_successful_steps += 1

            continuous_sol.append(this_interp_hb10)
            lower_continuous_sol.append(this_interp_hb8)

            if max_defect < (tol / 10):
                h *= 2
        else:
            h /= 2

    print("tolerance=", tol)
    print("Monitor8\n===================================")
    monitor8.print()
    print("\n\nMonitor10\n===============================")
    monitor10.print()

    print("================================\n")
    return (
        res, 
        continuous_sol.eval,
        continuous_sol.prime,
        continuous_sol.create_error_samplings(),
    )
