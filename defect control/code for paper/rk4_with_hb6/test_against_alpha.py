
from HB6 import HB
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

class Monitor:
    def __init__(self) -> None:
        self.different_values_alpha = set()
        self.n_steps=0
        self.n_successful_steps=0
    def print(self):
        print("alpha values", list(self.different_values_alpha))
        print("n_steps", self.n_steps)
        print("n_successful_steps", self.n_successful_steps)

def perfect_convergence(alpha, model, solution, model_num):
    monitor = Monitor()
    # the_xs = [0, 0.1, 1, 5, 9]
    the_xs = [5]
    for x0 in the_xs:
        convergences = []
        for h in [1e0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7]:
            (x_i_minus_1, y_i_minus_1) =    x0, solution([x0])[0]
            (x_i, y_i)                 =    x_i_minus_1 + alpha*h, solution([x_i_minus_1 + alpha*h])[0]
            (x_i_plus_1, y_i_plus_1)   =    x_i + h, solution([x_i + h])[0]

            f_i_minus_1 =   model(x_i_minus_1, y_i_minus_1)[0]
            f_i         =   model(x_i, y_i)[0]   
            f_i_plus_1  =   model(x_i_plus_1, y_i_plus_1)[0]

            this_interp = HB(
                x_i_minus_1, x_i, x_i_plus_1,
                y_i_minus_1, f_i_minus_1,
                y_i, f_i,
                y_i_plus_1, f_i_plus_1,
                monitor
            )

            pts_to_sample, indices = create_pts(x_i_minus_1, x_i_plus_1)
            defects = []
            defects_horner = []
            for pt in pts_to_sample:
                y = solution([pt])[0]
                # y = this_hb(pt)
                f_eval  = model(pt, y)[0]
                hb_prime_eval = this_interp.prime(pt)
                hb_prime_horner_eval = this_interp.prime_horner(pt)
                defects.append( abs(hb_prime_eval - f_eval) )
                defects_horner.append( abs( hb_prime_horner_eval - f_eval ) )
            # plt.figure()
            # plt.plot(indices, defects, label=f"defect_{str(x_i_minus_1)}_{str(x_i_plus_1)}")
            # plt.title(f"h = {h}")
            # plt.show()

            convergences.append( 
                ( h, max(defects), max(defects_horner) )
            )

        hs = [-log10(convergence[0]) for convergence in convergences]
        max_defects = [log10(abs(convergence[1])) for convergence in convergences]
        # max_defects_horner = [log10(abs(convergence[2])) for convergence in convergences]
        # plt.figure()
        # plt.plot(hs, max_defects, label="h vs max_defect")
        # plt.plot(hs, max_defects_horner, label="h vs max_defect_horner")
        # plt.ylabel("log of defect")
        # plt.xlabel("log of hs")
        # plt.title(f"h vs max defect at x0={x0}, alpha={alpha}, model={model_num}")
        # plt.legend()
        # plt.show()
        # monitor.print()

        # print("h", "\t\t", "max_defect", "\t\t", "log(h)", "\t\t", "log(defect)")
        # for (h, max_defect) in convergences:
        #    print(h, "\t\t", max_defect, "\t\t", log10(h), "\t\t", log10(max_defect)) 
        return (hs, max_defects)

def experiments(model, solution, model_num):
    (hs_1, max_defects_1) = perfect_convergence(1, model, solution, model_num)
    (hs_2, max_defects_2) = perfect_convergence(2, model, solution, model_num)
    (hs_1_2, max_defects_1_2) = perfect_convergence(1/2, model, solution, model_num)
    (hs_4, max_defects_4) = perfect_convergence(4, model, solution, model_num)
    (hs_1_4, max_defects_1_4) = perfect_convergence(1/4, model, solution, model_num)
    # ALPHA IS NOT ALLOWED TO BE TOO LARGE => the closer alpha is to 1, the better the codes...
    (hs_256, max_defects_256) = perfect_convergence(2**8, model, solution, model_num)
    (hs_1_256, max_defects_1_256) = perfect_convergence(1/(2**8), model, solution, model_num)

    plt.figure()

    plt.plot(hs_1, max_defects_1, label="alpha=1")
    
    plt.plot(hs_2, max_defects_2, label="alpha=2")
    plt.plot(hs_1_2, max_defects_1_2, label="alpha=1/2")

    plt.plot(hs_4, max_defects_4, label="alpha=4")
    plt.plot(hs_4, max_defects_1_4, label="alpha=1/4")

    plt.plot(hs_1, max_defects_256, label="alpha=256")
    plt.plot(hs_1, max_defects_1_256, label="alpha=1/256")

    plt.ylabel("log of defect")
    plt.xlabel("- log of hs")
    plt.title(f"h vs max defect")
    plt.legend(loc="lower right")
    plt.show()


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

# experiments(model1, solution1, 1)
experiments(model2, solution2, 2)
# experiments(model3, solution3, 3)
# experiments(model6, solution6, 6)
# experiments(model7, solution7, 7)
# experiments(model11, solution11, 11)
