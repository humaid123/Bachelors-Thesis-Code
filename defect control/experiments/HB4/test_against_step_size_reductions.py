from HB4 import HB
from math import exp, log10, sqrt, sin, cos
import matplotlib.pyplot as plt

def create_pts(a, b):
    num_points = 100
    curr = a
    res = [curr]
    h = (b - a) / num_points
    x = 0
    xs = [x]
    for i in range(1, num_points + 1):
        curr += h
        res.append(curr)
        x += h
        xs.append(x)
    return res, xs

def perfect_convergence(model, solution, model_num):
    # the_xs = [0, 0.1, 1, 5, 9]
    the_xs = [5]
    for x0 in the_xs:
        convergences = []
        for h in [1e0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7]:
            (x_i, y_i)                = x0, solution([x0])[0]
            (x_i_plus_1, y_i_plus_1)  = x_i + h, solution([x_i + h])[0]

            f_i         =   model(x_i, y_i)[0]   
            f_i_plus_1  =   model(x_i_plus_1, y_i_plus_1)[0]

            this_interp = HB(
                x_i, x_i_plus_1, 
                y_i, f_i, 
                y_i_plus_1, f_i_plus_1
            )

            pts_to_sample, indices = create_pts(x_i, x_i_plus_1)
            defects = []
            defects_horner = []
            defects_bary = []
            for pt in pts_to_sample:
                y = solution([pt])[0]
                # y = this_hb(pt)
                f_eval  = model(pt, y)[0]
                hb_prime_eval = this_interp.prime(pt)
                hb_prime_horner_eval = this_interp.prime_horner(pt)
                hb_prime_bary_eval = this_interp.prime_bary(pt)
                defects.append( abs(hb_prime_eval - f_eval) )
                defects_horner.append( abs( hb_prime_horner_eval - f_eval ) )
                defects_bary.append( abs( hb_prime_bary_eval - f_eval ) )
            # plt.figure()
            # plt.plot(indices, defects, label=f"defect_{str(x_i_minus_1)}_{str(x_i_plus_1)}")
            # plt.title(f"h = {h}")
            # plt.show()

            convergences.append( 
                ( h, max(defects), max(defects_horner), max(defects_bary) )
            )

        hs = [-log10(convergence[0]) for convergence in convergences]
        max_defects = [log10(abs(convergence[1])) for convergence in convergences]
        max_defects_horner = [log10(abs(convergence[2])) for convergence in convergences]
        max_defects_bary = [log10(abs(convergence[3])) for convergence in convergences]
        plt.figure()
        plt.plot(hs, max_defects, label="h vs max_defect")
        plt.plot(hs, max_defects_horner, label="h vs max_defect_horner")
        plt.plot(hs, max_defects_bary, label="h vs max_defect_bary")
        plt.ylabel("log of defect")
        plt.xlabel("log of hs")
        plt.title(f"h vs max defect at x0={x0}, model={model_num}")
        plt.legend()
        plt.show()
    return (hs, max_defects, max_defects_horner, max_defects_bary)

def model1(t, y):
    return [(-1/2) * y**3]
def solution1(t):
    return [1/sqrt(1+x) for x in t]

def model2(t, y):
    return [-2*t*y**2]
def solution2(t):
    return [1/(1+x**2) for x in t]

def model3(t, y):
    return [(1/4)*y*(1-y/20)]
def solution3(t):
    return [20/(1 + 19*exp(-x/4)) for x in t]

def model6(t, y):
    return [-y/(t+1)]
def solution6(t):
    return [1/(x+1) for x in t]

def model7(t, y):
    alpha = 0.1
    return [ -alpha*y - exp(-alpha*t)*sin(t)]
def solution7(t):
    alpha = 0.1
    return [exp(-alpha*x)*cos(x) for x in t]

def model11(t, y):
    return [-2*y + t]
def solution11(t):
    return [1/4 * (-1 + 5 * exp(-2 * x) + 2 * x) for x in t]

def experiments():
    (hs_1, max_defects_1, max_defects_horner_1, max_defects_bary_1) = perfect_convergence(model1, solution1, 1)
    (hs_2, max_defects_2, max_defects_horner_2, max_defects_bary_2) = perfect_convergence(model2, solution2, 2)
    (hs_3, max_defects_3, max_defects_horner_3, max_defects_bary_3) = perfect_convergence(model3, solution3, 3)
    (hs_6, max_defects_6, max_defects_horner_6, max_defects_bary_6) = perfect_convergence(model6, solution6, 6)
    (hs_7, max_defects_7, max_defects_horner_7, max_defects_bary_7) = perfect_convergence(model7, solution7, 7)
    (hs_11, max_defects_11, max_defects_horner_11, max_defects_bary_11) = perfect_convergence(model11, solution11, 11)

    plt.figure()

    plt.plot(hs_1, max_defects_1, label="h vs max_defect, model=1")
    plt.plot(hs_1, max_defects_horner_1, label="h vs max_defect_horner, model=1")
    plt.plot(hs_1, max_defects_bary_1, label="h vs max_defect_bary, model=1")
    
    plt.plot(hs_2, max_defects_2, label="h vs max_defect, model=2")
    plt.plot(hs_2, max_defects_horner_2, label="h vs max_defect_horner, model=2")
    plt.plot(hs_2, max_defects_bary_2, label="h vs max_defect_bary, model=2")

    plt.plot(hs_3, max_defects_3, label="h vs max_defect, model=3")
    plt.plot(hs_3, max_defects_horner_3, label="h vs max_defect_horner, model=3")
    plt.plot(hs_3, max_defects_bary_3, label="h vs max_defect_bary, model=3")

    plt.plot(hs_6, max_defects_6, label="h vs max_defect, model=6")
    plt.plot(hs_6, max_defects_horner_6, label="h vs max_defect_horner, model=6")
    plt.plot(hs_6, max_defects_bary_6, label="h vs max_defect_bary, model=6")
    
    plt.plot(hs_7, max_defects_7, label="h vs max_defect, model=7")
    plt.plot(hs_7, max_defects_horner_7, label="h vs max_defect_horner, model=7")
    plt.plot(hs_7, max_defects_bary_7, label="h vs max_defect_bary, model=7")

    plt.plot(hs_11, max_defects_11, label="h vs max_defect, model=11")
    plt.plot(hs_11, max_defects_horner_11, label="h vs max_defect_horner, model=11")
    plt.plot(hs_11, max_defects_bary_11, label="h vs max_defect_bary, model=11")

    plt.ylabel("log of defect")
    plt.xlabel("log of hs")
    plt.title(f"h vs max defect")
    plt.legend(loc="upper right")
    plt.show()

experiments()