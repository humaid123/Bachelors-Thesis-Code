
from math import cos, pi
from scipy.interpolate import BarycentricInterpolator

def get_Chebyshev_nodes(a, b, n):
    res = []
    for k in range(1, n+1):
        res.append(
            (a+b)/2 + (b-a)/2 * cos( (2*k - 1) / (2*n) * pi)
        )
    return res

class HB8ContinuousSolution:
    def __init__(self) -> None:
        self.interps = []
    
    def eval(self, x) -> float:
        for hb in self.interps:
            if (hb.x_i <= x <= hb.x_i_plus_1):
                return hb.eval(x)

        first_hb = self.interps[0]
        if (first_hb.x_i_minus_2 <= x <= first_hb.x_i):
            return first_hb.eval(x)
        print(f"ERROR: {x} is outside of the solution range: {first_hb.x_i_minus_2} <= x <= {self.interps[-1].x_i_plus_1}")
        return -1

    def prime(self, x) -> float:
        for hb in self.interps:
            if (hb.x_i <= x <= hb.x_i_plus_1):
                return hb.prime(x)

        first_hb = self.interps[0]
        if (first_hb.x_i_minus_2 <= x <= first_hb.x_i):
            return first_hb.prime(x)
        
        print(f"ERROR: {x} is outside of the solution range: {first_hb.x_i_minus_2} <= x <= {self.interps[-1].x_i_plus_1}")
        return -1
    
    def append(self, interp) -> None:
        self.interps.append(interp)

    def extend(self, newInterps) -> None:
        self.interps.extend(newInterps)

    def create_error_samplings(self):
        res = []
        for interp in self.interps:
            interp: HB8 = interp 
            res.append(
                (interp.x_i, interp.x_i_plus_1, interp)
            )
        return res


def create_continuous_sol_from_results(res, fn_s, monitor):
    interps = []
    for i in range(len(res) - 3):
        x_i_minus_2, y_i_minus_2 = res[i]    
        x_i_minus_1, y_i_minus_1 = res[i + 1]    
        x_i, y_i                 = res[i + 2]    
        x_i_plus_1, y_i_plus_1   = res[i + 3]

        f_i_minus_2 = fn_s[i]    
        f_i_minus_1 = fn_s[i + 1]    
        f_i         = fn_s[i + 2]    
        f_i_plus_1  = fn_s[i + 3]
        
        interps.append(
            HB8 (
                x_i_minus_2, x_i_minus_1, x_i, x_i_plus_1,
                y_i_minus_2, f_i_minus_2,
                y_i_minus_1, f_i_minus_1,
                y_i, f_i,
                y_i_plus_1, f_i_plus_1,
                monitor 
            )
        )
    continuous_sol = HB8ContinuousSolution()
    continuous_sol.extend(interps)
    return continuous_sol

def create_defect_samplings(res, fn_s, monitor):
    result = []
    for i in range(len(res) - 3):
        x_i_minus_2, y_i_minus_2 = res[i]    
        x_i_minus_1, y_i_minus_1 = res[i + 1]    
        x_i, y_i                 = res[i + 2]    
        x_i_plus_1, y_i_plus_1   = res[i + 3]

        f_i_minus_2 = fn_s[i]    
        f_i_minus_1 = fn_s[i + 1]    
        f_i         = fn_s[i + 2]    
        f_i_plus_1  = fn_s[i + 3]
        
        interp = HB8 (
                x_i_minus_2, x_i_minus_1, x_i, x_i_plus_1,
                y_i_minus_2, f_i_minus_2,
                y_i_minus_1, f_i_minus_1,
                y_i, f_i,
                y_i_plus_1, f_i_plus_1,
                monitor
        )
        result.append( (x_i_minus_1, x_i, x_i_plus_1, interp) )
    return result
        
def d0(x, alpha, beta):  
    return ( (6*alpha**2 + 4*alpha*beta + 8*alpha + 2*beta + 2)/(alpha**9 + 3*alpha**8*beta + 6*alpha**8 + 3*alpha**7*beta**2 + 15*alpha**7*beta + 15*alpha**7 + alpha**6*beta**3 + 12*alpha**6*beta**2 + 30*alpha**6*beta + 20*alpha**6 + 3*alpha**5*beta**3 + 18*alpha**5*beta**2 + 30*alpha**5*beta + 15*alpha**5 + 3*alpha**4*beta**3 + 12*alpha**4*beta**2 + 15*alpha**4*beta + 6*alpha**4 + alpha**3*beta**3 + 3*alpha**3*beta**2 + 3*alpha**3*beta + alpha**3)*(x**7) + (7*alpha**3 - 7*alpha**2*beta + 28*alpha**2 - 8*alpha*beta**2 - alpha*beta + 27*alpha - 4*beta**2 + 2*beta + 6)/(alpha**9 + 3*alpha**8*beta + 6*alpha**8 + 3*alpha**7*beta**2 + 15*alpha**7*beta + 15*alpha**7 + alpha**6*beta**3 + 12*alpha**6*beta**2 + 30*alpha**6*beta + 20*alpha**6 + 3*alpha**5*beta**3 + 18*alpha**5*beta**2 + 30*alpha**5*beta + 15*alpha**5 + 3*alpha**4*beta**3 + 12*alpha**4*beta**2 + 15*alpha**4*beta + 6*alpha**4 + alpha**3*beta**3 + 3*alpha**3*beta**2 + 3*alpha**3*beta + alpha**3)*(x**6) + (-14*alpha**3*beta + 14*alpha**3 - 4*alpha**2*beta**2 - 46*alpha**2*beta + 38*alpha**2 + 4*alpha*beta**3 - 22*alpha*beta**2 - 36*alpha*beta + 30*alpha + 2*beta**3 - 10*beta**2 - 6*beta + 6)/(alpha**9 + 3*alpha**8*beta + 6*alpha**8 + 3*alpha**7*beta**2 + 15*alpha**7*beta + 15*alpha**7 + alpha**6*beta**3 + 12*alpha**6*beta**2 + 30*alpha**6*beta + 20*alpha**6 + 3*alpha**5*beta**3 + 18*alpha**5*beta**2 + 30*alpha**5*beta + 15*alpha**5 + 3*alpha**4*beta**3 + 12*alpha**4*beta**2 + 15*alpha**4*beta + 6*alpha**4 + alpha**3*beta**3 + 3*alpha**3*beta**2 + 3*alpha**3*beta + alpha**3)*(x**5) + (7*alpha**3*beta**2 - 28*alpha**3*beta + 7*alpha**3 + 5*alpha**2*beta**3 + 
8*alpha**2*beta**2 - 71*alpha**2*beta + 16*alpha**2 + 15*alpha*beta**3 - 9*alpha*beta**2 - 53*alpha*beta + 11*alpha + 6*beta**3 - 6*beta**2 - 10*beta + 2)/(alpha**9 + 3*alpha**8*beta + 6*alpha**8 + 3*alpha**7*beta**2 + 15*alpha**7*beta + 15*alpha**7 + alpha**6*beta**3 + 12*alpha**6*beta**2 + 30*alpha**6*beta + 20*alpha**6 + 3*alpha**5*beta**3 + 18*alpha**5*beta**2 + 30*alpha**5*beta + 15*alpha**5 + 3*alpha**4*beta**3 + 12*alpha**4*beta**2 + 15*alpha**4*beta + 6*alpha**4 + alpha**3*beta**3 + 3*alpha**3*beta**2 + 3*alpha**3*beta + alpha**3)*(x**4) + (14*alpha**3*beta**2 - 14*alpha**3*beta + 10*alpha**2*beta**3 + 28*alpha**2*beta**2 - 32*alpha**2*beta + 18*alpha*beta**3 + 16*alpha*beta**2 - 22*alpha*beta + 6*beta**3 + 2*beta**2 - 4*beta)/(alpha**9 + 3*alpha**8*beta + 6*alpha**8 + 3*alpha**7*beta**2 + 15*alpha**7*beta + 15*alpha**7 + alpha**6*beta**3 + 12*alpha**6*beta**2 + 30*alpha**6*beta + 20*alpha**6 + 3*alpha**5*beta**3 + 18*alpha**5*beta**2 + 30*alpha**5*beta + 15*alpha**5 + 3*alpha**4*beta**3 + 12*alpha**4*beta**2 + 15*alpha**4*beta + 6*alpha**4 + alpha**3*beta**3 + 3*alpha**3*beta**2 + 3*alpha**3*beta + alpha**3)*(x**3) + (7*alpha**2*beta**2 + 5*alpha*beta**3 + 9*alpha*beta**2 + 2*beta**3 + 2*beta**2)/(alpha**8 + 3*alpha**7*beta + 5*alpha**7 + 3*alpha**6*beta**2 + 12*alpha**6*beta + 10*alpha**6 + alpha**5*beta**3 + 9*alpha**5*beta**2 + 18*alpha**5*beta + 10*alpha**5 + 2*alpha**4*beta**3 + 9*alpha**4*beta**2 + 12*alpha**4*beta + 
5*alpha**4 + alpha**3*beta**3 + 3*alpha**3*beta**2 + 3*alpha**3*beta + alpha**3)*(x**2) + 0*(x) + 0 )
def d1(x, alpha, beta):
    return ( 1/(alpha**6 + 2*alpha**5*beta + 4*alpha**5 + alpha**4*beta**2 + 6*alpha**4*beta + 6*alpha**4 + 2*alpha**3*beta**2 + 6*alpha**3*beta + 4*alpha**3 + alpha**2*beta**2 + 2*alpha**2*beta + alpha**2)*(x**7) + (alpha - 2*beta + 3)/(alpha**6 + 2*alpha**5*beta + 4*alpha**5 + alpha**4*beta**2 + 6*alpha**4*beta + 6*alpha**4 + 2*alpha**3*beta**2 + 6*alpha**3*beta + 4*alpha**3 + alpha**2*beta**2 + 2*alpha**2*beta + alpha**2)*(x**6) + (-2*alpha*beta + 2*alpha + beta**2 - 6*beta + 3)/(alpha**6 + 2*alpha**5*beta + 4*alpha**5 + alpha**4*beta**2 + 6*alpha**4*beta + 6*alpha**4 + 2*alpha**3*beta**2 + 6*alpha**3*beta + 4*alpha**3 + alpha**2*beta**2 + 2*alpha**2*beta + alpha**2)*(x**5) + (alpha*beta**2 - 4*alpha*beta + alpha + 
3*beta**2 - 6*beta + 1)/(alpha**6 + 2*alpha**5*beta + 4*alpha**5 + alpha**4*beta**2 + 6*alpha**4*beta + 6*alpha**4 + 2*alpha**3*beta**2 + 6*alpha**3*beta + 4*alpha**3 + alpha**2*beta**2 + 2*alpha**2*beta + alpha**2)*(x**4) + (2*alpha*beta**2 - 2*alpha*beta + 3*beta**2 - 2*beta)/(alpha**6 + 2*alpha**5*beta + 4*alpha**5 + alpha**4*beta**2 + 6*alpha**4*beta + 6*alpha**4 + 2*alpha**3*beta**2 + 6*alpha**3*beta + 4*alpha**3 + alpha**2*beta**2 + 2*alpha**2*beta + alpha**2)*(x**3) + beta**2/(alpha**5 + 2*alpha**4*beta + 3*alpha**4 + alpha**3*beta**2 + 4*alpha**3*beta + 3*alpha**3 + alpha**2*beta**2 + 2*alpha**2*beta + alpha**2)*(x**2) + 0*(x) + 0 )
def d2(x, alpha, beta):
    return ( (2*alpha*beta + 4*alpha - 2*beta - 2)/(alpha**3*beta**3 + 3*alpha**3*beta**2 + 3*alpha**3*beta + alpha**3)*(x**7) + (4*alpha**2*beta + 8*alpha**2 - 4*alpha*beta**2 - 5*alpha*beta + 9*alpha + 4*beta**2 - 2*beta - 6)/(alpha**3*beta**3 + 3*alpha**3*beta**2 + 3*alpha**3*beta + alpha**3)*(x**6) + (2*alpha**3*beta + 4*alpha**3 - 8*alpha**2*beta**2 - 8*alpha**2*beta + 16*alpha**2 + 2*alpha*beta**3 - 2*alpha*beta**2 - 18*alpha*beta + 6*alpha - 2*beta**3 + 10*beta**2 + 6*beta - 6)/(alpha**3*beta**3 + 3*alpha**3*beta**2 + 3*alpha**3*beta + alpha**3)*(x**5) + (-4*alpha**3*beta**2 - 5*alpha**3*beta + 5*alpha**3 + 4*alpha**2*beta**3 - 8*alpha**2*beta**2 - 28*alpha**2*beta + 8*alpha**2 + 3*alpha*beta**3 + 9*alpha*beta**2 - 13*alpha*beta + alpha - 6*beta**3 + 6*beta**2 + 10*beta - 2)/(alpha**3*beta**3 + 3*alpha**3*beta**2 + 3*alpha**3*beta + alpha**3)*(x**4) + (2*alpha**3*beta**3 - 2*alpha**3*beta**2 - 10*alpha**3*beta + 8*alpha**2*beta**3 + 8*alpha**2*beta**2 - 16*alpha**2*beta + 8*alpha*beta**2 - 2*alpha*beta - 6*beta**3 - 2*beta**2 + 4*beta)/(alpha**3*beta**3 + 3*alpha**3*beta**2 + 3*alpha**3*beta + alpha**3)*(x**3) + (3*alpha**3*beta**3 + 5*alpha**3*beta**2 + 4*alpha**2*beta**3 + 8*alpha**2*beta**2 - alpha*beta**3 + alpha*beta**2 - 2*beta**3 - 2*beta**2)/(alpha**3*beta**3 + 3*alpha**3*beta**2 + 3*alpha**3*beta + alpha**3)*(x**2) + 0*(x) + 0 )
def d3(x, alpha, beta):
    return ( 1/(alpha**2*beta**2 + 2*alpha**2*beta + alpha**2)*(x**7) + (2*alpha - 2*beta + 3)/(alpha**2*beta**2 + 2*alpha**2*beta + alpha**2)*(x**6) + (alpha**2 - 4*alpha*beta + 4*alpha + beta**2 - 6*beta + 3)/(alpha**2*beta**2 + 2*alpha**2*beta + alpha**2)*(x**5) + (-2*alpha**2*beta + alpha**2 + 2*alpha*beta**2 - 8*alpha*beta + 2*alpha + 3*beta**2 - 6*beta + 1)/(alpha**2*beta**2 + 2*alpha**2*beta + alpha**2)*(x**4) + (alpha**2*beta**2 - 2*alpha**2*beta + 4*alpha*beta**2 - 4*alpha*beta + 3*beta**2 - 2*beta)/(alpha**2*beta**2 + 2*alpha**2*beta + alpha**2)*(x**3) + (alpha**2*beta**2 + 2*alpha*beta**2 + beta**2)/(alpha**2*beta**2 + 2*alpha**2*beta + alpha**2)*(x**2) + 0*(x) + 0 )
def d4(x, alpha, beta):
    return ( (-2*alpha*beta + 2*alpha - 4*beta + 2)/(alpha**3*beta**3 + 3*alpha**2*beta**3 + 3*alpha*beta**3 + beta**3)*(x**7) + (-4*alpha**2*beta + 4*alpha**2 + 4*alpha*beta**2 - 19*alpha*beta + 12*alpha + 8*beta**2 - 19*beta + 8)/(alpha**3*beta**3 + 3*alpha**2*beta**3 + 3*alpha*beta**3 + beta**3)*(x**6) + (-2*alpha**3*beta + 2*alpha**3 + 8*alpha**2*beta**2 - 22*alpha**2*beta + 14*alpha**2 - 2*alpha*beta**3 + 32*alpha*beta**2 - 54*alpha*beta + 24*alpha - 4*beta**3 + 32*beta**2 - 36*beta + 12)/(alpha**3*beta**3 + 3*alpha**2*beta**3 + 3*alpha*beta**3 + beta**3)*(x**5) + (4*alpha**3*beta**2 - 7*alpha**3*beta + 4*alpha**3 - 4*alpha**2*beta**3 + 32*alpha**2*beta**2 - 41*alpha**2*beta + 16*alpha**2 - 15*alpha*beta**3 + 72*alpha*beta**2 - 68*alpha*beta + 20*alpha - 15*beta**3 + 48*beta**2 - 34*beta + 8)/(alpha**3*beta**3 + 3*alpha**2*beta**3 + 3*alpha*beta**3 + beta**3)*(x**4) + (-2*alpha**3*beta**3 + 8*alpha**3*beta**2 - 8*alpha**3*beta + 2*alpha**3 - 14*alpha**2*beta**3 + 40*alpha**2*beta**2 - 32*alpha**2*beta + 6*alpha**2 - 30*alpha*beta**3 + 64*alpha*beta**2 - 40*alpha*beta + 6*alpha - 20*beta**3 + 32*beta**2 - 16*beta + 2)/(alpha**3*beta**3 + 3*alpha**2*beta**3 + 3*alpha*beta**3 + beta**3)*(x**3) + (-3*alpha**2*beta**2 + 4*alpha**2*beta - 3*alpha**2 - 10*alpha*beta**2 + 12*alpha*beta - 6*alpha - 10*beta**2 + 8*beta - 3)/(alpha**2*beta**2 + 2*alpha*beta**2 + beta**2)*(x**2) + 0*(x) + 1 )
def d5(x, alpha, beta):
    return ( 1/(alpha**2*beta**2 + 2*alpha*beta**2 + beta**2)*(x**7) + (2*alpha - 2*beta + 4)/(alpha**2*beta**2 + 2*alpha*beta**2 + beta**2)*(x**6) + (alpha**2 - 4*alpha*beta + 6*alpha + beta**2 - 8*beta + 6)/(alpha**2*beta**2 + 2*alpha*beta**2 + beta**2)*(x**5) + (-2*alpha**2*beta + 2*alpha**2 + 2*alpha*beta**2 - 12*alpha*beta + 6*alpha + 4*beta**2 - 12*beta + 
4)/(alpha**2*beta**2 + 2*alpha*beta**2 + beta**2)*(x**4) + (alpha**2*beta**2 - 4*alpha**2*beta + alpha**2 + 6*alpha*beta**2 - 12*alpha*beta + 2*alpha + 6*beta**2 - 8*beta + 1)/(alpha**2*beta**2 + 2*alpha*beta**2 + beta**2)*(x**3) + (2*alpha*beta - 2*alpha + 4*beta - 2)/(alpha*beta + beta)*(x**2) + 1*(x) + 0 )
def d6(x, alpha, beta):
    return ( (-4*alpha*beta - 2*alpha - 6*beta**2 - 8*beta - 2)/(alpha**3*beta**6 + 3*alpha**3*beta**5 + 3*alpha**3*beta**4 + alpha**3*beta**3 + 3*alpha**2*beta**7 + 12*alpha**2*beta**6 + 18*alpha**2*beta**5 + 12*alpha**2*beta**4 + 3*alpha**2*beta**3 + 3*alpha*beta**8 + 15*alpha*beta**7 + 30*alpha*beta**6 + 30*alpha*beta**5 + 15*alpha*beta**4 + 3*alpha*beta**3 + beta**9 + 6*beta**8 + 15*beta**7 + 20*beta**6 + 15*beta**5 + 6*beta**4 + beta**3)*(x**7) + (-8*alpha**2*beta - 4*alpha**2 - 7*alpha*beta**2 - 29*alpha*beta - 12*alpha + 7*beta**3 - 
14*beta**2 - 29*beta - 8)/(alpha**3*beta**6 + 3*alpha**3*beta**5 + 3*alpha**3*beta**4 + alpha**3*beta**3 + 3*alpha**2*beta**7 + 12*alpha**2*beta**6 + 18*alpha**2*beta**5 + 12*alpha**2*beta**4 + 3*alpha**2*beta**3 + 3*alpha*beta**8 + 15*alpha*beta**7 + 30*alpha*beta**6 + 30*alpha*beta**5 + 15*alpha*beta**4 + 3*alpha*beta**3 + beta**9 + 6*beta**8 + 15*beta**7 + 20*beta**6 + 15*beta**5 + 6*beta**4 + beta**3)*(x**6) + (-4*alpha**3*beta - 2*alpha**3 + 4*alpha**2*beta**2 - 26*alpha**2*beta - 14*alpha**2 + 14*alpha*beta**3 + 4*alpha*beta**2 - 54*alpha*beta - 24*alpha + 28*beta**3 + 4*beta**2 - 36*beta - 12)/(alpha**3*beta**6 + 3*alpha**3*beta**5 + 3*alpha**3*beta**4 + alpha**3*beta**3 + 3*alpha**2*beta**7 + 12*alpha**2*beta**6 + 18*alpha**2*beta**5 + 12*alpha**2*beta**4 + 3*alpha**2*beta**3 + 3*alpha*beta**8 + 15*alpha*beta**7 + 30*alpha*beta**6 + 30*alpha*beta**5 + 15*alpha*beta**4 + 3*alpha*beta**3 + beta**9 + 6*beta**8 + 15*beta**7 + 20*beta**6 + 15*beta**5 + 6*beta**4 + beta**3)*(x**5) + (5*alpha**3*beta**2 - 5*alpha**3*beta - 4*alpha**3 + 7*alpha**2*beta**3 + 28*alpha**2*beta**2 - 19*alpha**2*beta - 16*alpha**2 + 42*alpha*beta**3 + 54*alpha*beta**2 - 28*alpha*beta - 20*alpha + 42*beta**3 + 36*beta**2 - 14*beta - 8)/(alpha**3*beta**6 + 3*alpha**3*beta**5 + 
3*alpha**3*beta**4 + alpha**3*beta**3 + 3*alpha**2*beta**7 + 12*alpha**2*beta**6 + 18*alpha**2*beta**5 + 12*alpha**2*beta**4 + 3*alpha**2*beta**3 + 3*alpha*beta**8 + 15*alpha*beta**7 + 30*alpha*beta**6 + 30*alpha*beta**5 + 15*alpha*beta**4 + 3*alpha*beta**3 + beta**9 + 6*beta**8 + 15*beta**7 + 20*beta**6 + 15*beta**5 + 6*beta**4 + beta**3)*(x**4) + (10*alpha**3*beta**2 + 2*alpha**3*beta - 2*alpha**3 + 14*alpha**2*beta**3 + 44*alpha**2*beta**2 + 8*alpha**2*beta - 6*alpha**2 + 42*alpha*beta**3 + 68*alpha*beta**2 + 10*alpha*beta - 6*alpha + 28*beta**3 + 34*beta**2 + 4*beta - 2)/(alpha**3*beta**6 + 3*alpha**3*beta**5 + 3*alpha**3*beta**4 + alpha**3*beta**3 + 3*alpha**2*beta**7 + 12*alpha**2*beta**6 + 18*alpha**2*beta**5 + 
12*alpha**2*beta**4 + 3*alpha**2*beta**3 + 3*alpha*beta**8 + 15*alpha*beta**7 + 30*alpha*beta**6 + 30*alpha*beta**5 + 15*alpha*beta**4 + 3*alpha*beta**3 + beta**9 + 6*beta**8 + 15*beta**7 + 20*beta**6 + 15*beta**5 + 6*beta**4 + beta**3)*(x**3) + (5*alpha**3*beta + 3*alpha**3 + 7*alpha**2*beta**2 + 20*alpha**2*beta + 9*alpha**2 + 14*alpha*beta**2 + 25*alpha*beta 
+ 9*alpha + 7*beta**2 + 10*beta + 3)/(alpha**3*beta**5 + 3*alpha**3*beta**4 + 3*alpha**3*beta**3 + alpha**3*beta**2 + 3*alpha**2*beta**6 + 12*alpha**2*beta**5 + 18*alpha**2*beta**4 + 12*alpha**2*beta**3 + 3*alpha**2*beta**2 + 3*alpha*beta**7 + 15*alpha*beta**6 + 30*alpha*beta**5 + 30*alpha*beta**4 + 15*alpha*beta**3 + 3*alpha*beta**2 + beta**8 + 6*beta**7 + 15*beta**6 + 20*beta**5 + 15*beta**4 + 6*beta**3 + beta**2)*(x**2) + 0*(x) + 0 )
def d7(x, alpha, beta):
    return ( 1/(alpha**2*beta**4 + 2*alpha**2*beta**3 + alpha**2*beta**2 + 2*alpha*beta**5 + 6*alpha*beta**4 + 6*alpha*beta**3 + 2*alpha*beta**2 + beta**6 + 4*beta**5 + 6*beta**4 + 4*beta**3 + beta**2)*(x**7) + (2*alpha - beta + 4)/(alpha**2*beta**4 + 2*alpha**2*beta**3 + alpha**2*beta**2 + 2*alpha*beta**5 + 6*alpha*beta**4 + 6*alpha*beta**3 + 2*alpha*beta**2 + 
beta**6 + 4*beta**5 + 6*beta**4 + 4*beta**3 + beta**2)*(x**6) + (alpha**2 - 2*alpha*beta + 6*alpha - 4*beta + 6)/(alpha**2*beta**4 + 2*alpha**2*beta**3 + alpha**2*beta**2 + 2*alpha*beta**5 + 6*alpha*beta**4 + 6*alpha*beta**3 + 2*alpha*beta**2 + beta**6 + 4*beta**5 + 6*beta**4 + 4*beta**3 + beta**2)*(x**5) + (-alpha**2*beta + 2*alpha**2 - 6*alpha*beta + 6*alpha - 6*beta + 4)/(alpha**2*beta**4 + 2*alpha**2*beta**3 + alpha**2*beta**2 + 2*alpha*beta**5 + 6*alpha*beta**4 + 6*alpha*beta**3 + 2*alpha*beta**2 + beta**6 + 4*beta**5 + 6*beta**4 + 4*beta**3 + beta**2)*(x**4) + (-2*alpha**2*beta + alpha**2 - 6*alpha*beta + 2*alpha - 4*beta + 1)/(alpha**2*beta**4 + 2*alpha**2*beta**3 + alpha**2*beta**2 + 2*alpha*beta**5 + 6*alpha*beta**4 + 6*alpha*beta**3 + 2*alpha*beta**2 + beta**6 + 4*beta**5 + 6*beta**4 + 4*beta**3 + beta**2)*(x**3) + (-alpha**2 - 2*alpha - 1)/(alpha**2*beta**3 + 2*alpha**2*beta**2 + alpha**2*beta + 2*alpha*beta**4 + 6*alpha*beta**3 + 6*alpha*beta**2 + 2*alpha*beta + beta**5 + 4*beta**4 + 6*beta**3 + 4*beta**2 + beta)*(x**2) + 0*(x) + 0 )
def d0_prime(x, alpha, beta):
    return ( 7*(6*alpha**2 + 4*alpha*beta + 8*alpha + 2*beta + 2)/(alpha**9 + 3*alpha**8*beta + 6*alpha**8 + 3*alpha**7*beta**2 + 15*alpha**7*beta + 15*alpha**7 + alpha**6*beta**3 + 
12*alpha**6*beta**2 + 30*alpha**6*beta + 20*alpha**6 + 3*alpha**5*beta**3 + 18*alpha**5*beta**2 + 30*alpha**5*beta + 15*alpha**5 + 3*alpha**4*beta**3 + 12*alpha**4*beta**2 + 15*alpha**4*beta + 6*alpha**4 + alpha**3*beta**3 + 3*alpha**3*beta**2 + 3*alpha**3*beta + alpha**3)*(x**6) + 6*(7*alpha**3 - 7*alpha**2*beta + 28*alpha**2 - 8*alpha*beta**2 - alpha*beta + 27*alpha - 4*beta**2 + 2*beta + 6)/(alpha**9 + 3*alpha**8*beta + 6*alpha**8 + 3*alpha**7*beta**2 + 15*alpha**7*beta + 15*alpha**7 + alpha**6*beta**3 + 12*alpha**6*beta**2 + 30*alpha**6*beta + 20*alpha**6 + 3*alpha**5*beta**3 + 18*alpha**5*beta**2 + 30*alpha**5*beta + 15*alpha**5 + 3*alpha**4*beta**3 + 12*alpha**4*beta**2 + 15*alpha**4*beta + 6*alpha**4 + alpha**3*beta**3 + 3*alpha**3*beta**2 + 3*alpha**3*beta + alpha**3)*(x**5) + 5*(-14*alpha**3*beta + 14*alpha**3 - 4*alpha**2*beta**2 - 46*alpha**2*beta + 38*alpha**2 + 4*alpha*beta**3 - 22*alpha*beta**2 - 36*alpha*beta + 30*alpha + 2*beta**3 - 10*beta**2 - 6*beta + 6)/(alpha**9 + 3*alpha**8*beta + 6*alpha**8 + 3*alpha**7*beta**2 + 15*alpha**7*beta + 15*alpha**7 + alpha**6*beta**3 + 12*alpha**6*beta**2 + 30*alpha**6*beta + 20*alpha**6 + 3*alpha**5*beta**3 + 18*alpha**5*beta**2 + 30*alpha**5*beta + 15*alpha**5 + 3*alpha**4*beta**3 + 12*alpha**4*beta**2 + 15*alpha**4*beta + 6*alpha**4 + alpha**3*beta**3 + 3*alpha**3*beta**2 + 3*alpha**3*beta + alpha**3)*(x**4) + 4*(7*alpha**3*beta**2 - 28*alpha**3*beta + 7*alpha**3 + 5*alpha**2*beta**3 + 8*alpha**2*beta**2 - 71*alpha**2*beta + 16*alpha**2 + 15*alpha*beta**3 - 9*alpha*beta**2 - 53*alpha*beta + 11*alpha + 6*beta**3 - 6*beta**2 - 10*beta + 2)/(alpha**9 + 3*alpha**8*beta + 6*alpha**8 + 3*alpha**7*beta**2 + 15*alpha**7*beta + 15*alpha**7 + alpha**6*beta**3 + 12*alpha**6*beta**2 + 30*alpha**6*beta + 20*alpha**6 + 3*alpha**5*beta**3 + 18*alpha**5*beta**2 + 30*alpha**5*beta + 15*alpha**5 + 3*alpha**4*beta**3 + 12*alpha**4*beta**2 + 15*alpha**4*beta + 6*alpha**4 + alpha**3*beta**3 + 3*alpha**3*beta**2 + 3*alpha**3*beta + alpha**3)*(x**3) + 3*(14*alpha**3*beta**2 - 14*alpha**3*beta + 10*alpha**2*beta**3 + 28*alpha**2*beta**2 - 32*alpha**2*beta + 18*alpha*beta**3 + 16*alpha*beta**2 - 22*alpha*beta + 6*beta**3 + 2*beta**2 - 4*beta)/(alpha**9 + 3*alpha**8*beta + 6*alpha**8 + 3*alpha**7*beta**2 + 15*alpha**7*beta + 15*alpha**7 + alpha**6*beta**3 + 12*alpha**6*beta**2 + 30*alpha**6*beta 
+ 20*alpha**6 + 3*alpha**5*beta**3 + 18*alpha**5*beta**2 + 30*alpha**5*beta + 15*alpha**5 + 3*alpha**4*beta**3 + 12*alpha**4*beta**2 + 15*alpha**4*beta + 6*alpha**4 + alpha**3*beta**3 + 3*alpha**3*beta**2 + 3*alpha**3*beta + alpha**3)*(x**2) + 2*(7*alpha**2*beta**2 + 5*alpha*beta**3 + 9*alpha*beta**2 + 2*beta**3 + 2*beta**2)/(alpha**8 + 3*alpha**7*beta + 5*alpha**7 + 3*alpha**6*beta**2 + 12*alpha**6*beta + 10*alpha**6 + alpha**5*beta**3 + 9*alpha**5*beta**2 + 18*alpha**5*beta + 10*alpha**5 + 2*alpha**4*beta**3 + 9*alpha**4*beta**2 + 12*alpha**4*beta + 5*alpha**4 + alpha**3*beta**3 + 3*alpha**3*beta**2 + 3*alpha**3*beta + alpha**3)*x + 0 )
def d1_prime(x, alpha, beta):
    return ( 7*1/(alpha**6 + 2*alpha**5*beta + 4*alpha**5 + alpha**4*beta**2 + 6*alpha**4*beta + 6*alpha**4 + 2*alpha**3*beta**2 + 6*alpha**3*beta + 4*alpha**3 + alpha**2*beta**2 + 2*alpha**2*beta + alpha**2)*(x**6) + 6*(alpha - 2*beta + 3)/(alpha**6 + 2*alpha**5*beta + 4*alpha**5 + alpha**4*beta**2 + 6*alpha**4*beta + 6*alpha**4 + 2*alpha**3*beta**2 + 6*alpha**3*beta + 4*alpha**3 + alpha**2*beta**2 + 2*alpha**2*beta + alpha**2)*(x**5) + 5*(-2*alpha*beta + 2*alpha + beta**2 - 6*beta + 3)/(alpha**6 + 2*alpha**5*beta + 4*alpha**5 + alpha**4*beta**2 + 6*alpha**4*beta + 6*alpha**4 + 2*alpha**3*beta**2 + 6*alpha**3*beta + 4*alpha**3 + alpha**2*beta**2 + 2*alpha**2*beta + alpha**2)*(x**4) + 4*(alpha*beta**2 - 4*alpha*beta + 
alpha + 3*beta**2 - 6*beta + 1)/(alpha**6 + 2*alpha**5*beta + 4*alpha**5 + alpha**4*beta**2 + 6*alpha**4*beta + 6*alpha**4 + 2*alpha**3*beta**2 + 6*alpha**3*beta + 4*alpha**3 + alpha**2*beta**2 + 2*alpha**2*beta + alpha**2)*(x**3) + 3*(2*alpha*beta**2 - 2*alpha*beta + 3*beta**2 - 2*beta)/(alpha**6 + 2*alpha**5*beta + 4*alpha**5 + alpha**4*beta**2 + 6*alpha**4*beta + 6*alpha**4 + 2*alpha**3*beta**2 + 6*alpha**3*beta + 4*alpha**3 + alpha**2*beta**2 + 2*alpha**2*beta + alpha**2)*(x**2) + 2*beta**2/(alpha**5 + 2*alpha**4*beta + 3*alpha**4 + alpha**3*beta**2 + 4*alpha**3*beta + 3*alpha**3 + alpha**2*beta**2 + 2*alpha**2*beta + alpha**2)*x + 0 )
def d2_prime(x, alpha, beta):
    return ( 7*(2*alpha*beta + 4*alpha - 2*beta - 2)/(alpha**3*beta**3 + 3*alpha**3*beta**2 + 3*alpha**3*beta + alpha**3)*(x**6) + 6*(4*alpha**2*beta + 8*alpha**2 - 4*alpha*beta**2 - 5*alpha*beta + 9*alpha + 4*beta**2 - 2*beta - 6)/(alpha**3*beta**3 + 3*alpha**3*beta**2 + 3*alpha**3*beta + alpha**3)*(x**5) + 5*(2*alpha**3*beta + 4*alpha**3 - 8*alpha**2*beta**2 - 8*alpha**2*beta + 16*alpha**2 + 2*alpha*beta**3 - 2*alpha*beta**2 - 18*alpha*beta + 6*alpha - 2*beta**3 + 10*beta**2 + 6*beta - 6)/(alpha**3*beta**3 + 3*alpha**3*beta**2 + 3*alpha**3*beta + alpha**3)*(x**4) + 4*(-4*alpha**3*beta**2 - 5*alpha**3*beta + 5*alpha**3 + 4*alpha**2*beta**3 - 8*alpha**2*beta**2 - 28*alpha**2*beta + 8*alpha**2 + 3*alpha*beta**3 + 9*alpha*beta**2 - 13*alpha*beta + alpha - 6*beta**3 + 6*beta**2 + 10*beta - 2)/(alpha**3*beta**3 + 3*alpha**3*beta**2 + 3*alpha**3*beta + alpha**3)*(x**3) + 3*(2*alpha**3*beta**3 - 2*alpha**3*beta**2 - 10*alpha**3*beta + 8*alpha**2*beta**3 + 8*alpha**2*beta**2 - 16*alpha**2*beta + 8*alpha*beta**2 - 2*alpha*beta - 6*beta**3 - 2*beta**2 + 4*beta)/(alpha**3*beta**3 + 3*alpha**3*beta**2 + 3*alpha**3*beta + alpha**3)*(x**2) + 2*(3*alpha**3*beta**3 + 5*alpha**3*beta**2 + 4*alpha**2*beta**3 + 8*alpha**2*beta**2 - alpha*beta**3 + alpha*beta**2 - 2*beta**3 - 2*beta**2)/(alpha**3*beta**3 + 3*alpha**3*beta**2 + 3*alpha**3*beta + alpha**3)*x + 0 )
def d3_prime(x, alpha, beta):
    return ( 7*1/(alpha**2*beta**2 + 2*alpha**2*beta + alpha**2)*(x**6) + 6*(2*alpha - 2*beta + 3)/(alpha**2*beta**2 + 2*alpha**2*beta + alpha**2)*(x**5) + 5*(alpha**2 - 4*alpha*beta + 4*alpha + beta**2 - 6*beta + 3)/(alpha**2*beta**2 + 2*alpha**2*beta + alpha**2)*(x**4) + 4*(-2*alpha**2*beta + alpha**2 + 2*alpha*beta**2 - 8*alpha*beta + 2*alpha + 3*beta**2 - 6*beta + 1)/(alpha**2*beta**2 + 2*alpha**2*beta + alpha**2)*(x**3) + 3*(alpha**2*beta**2 - 2*alpha**2*beta + 4*alpha*beta**2 - 4*alpha*beta + 3*beta**2 - 2*beta)/(alpha**2*beta**2 + 2*alpha**2*beta + alpha**2)*(x**2) + 2*(alpha**2*beta**2 + 2*alpha*beta**2 + beta**2)/(alpha**2*beta**2 + 2*alpha**2*beta + alpha**2)*x + 0 )
def d4_prime(x, alpha, beta):
    return ( 7*(-2*alpha*beta + 2*alpha - 4*beta + 2)/(alpha**3*beta**3 + 3*alpha**2*beta**3 + 3*alpha*beta**3 + beta**3)*(x**6) + 6*(-4*alpha**2*beta + 4*alpha**2 + 4*alpha*beta**2 
- 19*alpha*beta + 12*alpha + 8*beta**2 - 19*beta + 8)/(alpha**3*beta**3 + 3*alpha**2*beta**3 + 3*alpha*beta**3 + beta**3)*(x**5) + 5*(-2*alpha**3*beta + 2*alpha**3 + 8*alpha**2*beta**2 - 22*alpha**2*beta + 14*alpha**2 - 2*alpha*beta**3 + 32*alpha*beta**2 - 54*alpha*beta + 24*alpha - 4*beta**3 + 32*beta**2 - 36*beta + 12)/(alpha**3*beta**3 + 3*alpha**2*beta**3 + 
3*alpha*beta**3 + beta**3)*(x**4) + 4*(4*alpha**3*beta**2 - 7*alpha**3*beta + 4*alpha**3 - 4*alpha**2*beta**3 + 32*alpha**2*beta**2 - 41*alpha**2*beta + 16*alpha**2 - 15*alpha*beta**3 + 72*alpha*beta**2 - 68*alpha*beta + 20*alpha - 15*beta**3 + 48*beta**2 - 34*beta + 8)/(alpha**3*beta**3 + 3*alpha**2*beta**3 + 3*alpha*beta**3 + beta**3)*(x**3) + 3*(-2*alpha**3*beta**3 + 8*alpha**3*beta**2 - 8*alpha**3*beta + 2*alpha**3 - 14*alpha**2*beta**3 + 40*alpha**2*beta**2 - 32*alpha**2*beta + 6*alpha**2 - 30*alpha*beta**3 + 64*alpha*beta**2 - 40*alpha*beta + 6*alpha - 20*beta**3 + 32*beta**2 - 16*beta + 2)/(alpha**3*beta**3 + 3*alpha**2*beta**3 + 3*alpha*beta**3 + beta**3)*(x**2) + 2*(-3*alpha**2*beta**2 + 4*alpha**2*beta - 3*alpha**2 - 10*alpha*beta**2 + 12*alpha*beta - 6*alpha - 10*beta**2 + 8*beta - 3)/(alpha**2*beta**2 + 2*alpha*beta**2 + beta**2)*x + 0 )
def d5_prime(x, alpha, beta):
    return ( 7*1/(alpha**2*beta**2 + 2*alpha*beta**2 + beta**2)*(x**6) + 6*(2*alpha - 2*beta + 4)/(alpha**2*beta**2 + 2*alpha*beta**2 + beta**2)*(x**5) + 5*(alpha**2 - 4*alpha*beta + 6*alpha + beta**2 - 8*beta + 6)/(alpha**2*beta**2 + 2*alpha*beta**2 + beta**2)*(x**4) + 4*(-2*alpha**2*beta + 2*alpha**2 + 2*alpha*beta**2 - 12*alpha*beta + 6*alpha + 4*beta**2 - 12*beta + 4)/(alpha**2*beta**2 + 2*alpha*beta**2 + beta**2)*(x**3) + 3*(alpha**2*beta**2 - 4*alpha**2*beta + alpha**2 + 6*alpha*beta**2 - 12*alpha*beta + 2*alpha + 6*beta**2 - 8*beta + 1)/(alpha**2*beta**2 + 2*alpha*beta**2 + beta**2)*(x**2) + 2*(2*alpha*beta - 2*alpha + 4*beta - 2)/(alpha*beta + beta)*x + 1 )
def d6_prime(x, alpha, beta):
    return ( 7*(-4*alpha*beta - 2*alpha - 6*beta**2 - 8*beta - 2)/(alpha**3*beta**6 + 3*alpha**3*beta**5 + 3*alpha**3*beta**4 + alpha**3*beta**3 + 3*alpha**2*beta**7 + 12*alpha**2*beta**6 + 18*alpha**2*beta**5 + 12*alpha**2*beta**4 + 3*alpha**2*beta**3 + 3*alpha*beta**8 + 15*alpha*beta**7 + 30*alpha*beta**6 + 30*alpha*beta**5 + 15*alpha*beta**4 + 3*alpha*beta**3 + beta**9 + 6*beta**8 + 15*beta**7 + 20*beta**6 + 15*beta**5 + 6*beta**4 + beta**3)*(x**6) + 6*(-8*alpha**2*beta - 4*alpha**2 - 7*alpha*beta**2 - 29*alpha*beta - 12*alpha + 7*beta**3 - 14*beta**2 - 29*beta - 8)/(alpha**3*beta**6 + 3*alpha**3*beta**5 + 3*alpha**3*beta**4 + alpha**3*beta**3 + 3*alpha**2*beta**7 + 12*alpha**2*beta**6 + 18*alpha**2*beta**5 + 12*alpha**2*beta**4 + 3*alpha**2*beta**3 + 3*alpha*beta**8 + 15*alpha*beta**7 + 30*alpha*beta**6 + 30*alpha*beta**5 + 15*alpha*beta**4 + 3*alpha*beta**3 + beta**9 + 6*beta**8 + 15*beta**7 
+ 20*beta**6 + 15*beta**5 + 6*beta**4 + beta**3)*(x**5) + 5*(-4*alpha**3*beta - 2*alpha**3 + 4*alpha**2*beta**2 - 26*alpha**2*beta - 14*alpha**2 + 14*alpha*beta**3 + 4*alpha*beta**2 
- 54*alpha*beta - 24*alpha + 28*beta**3 + 4*beta**2 - 36*beta - 12)/(alpha**3*beta**6 + 3*alpha**3*beta**5 + 3*alpha**3*beta**4 + alpha**3*beta**3 + 3*alpha**2*beta**7 + 12*alpha**2*beta**6 + 18*alpha**2*beta**5 + 12*alpha**2*beta**4 + 3*alpha**2*beta**3 + 3*alpha*beta**8 + 15*alpha*beta**7 + 30*alpha*beta**6 + 30*alpha*beta**5 + 15*alpha*beta**4 + 3*alpha*beta**3 + beta**9 + 6*beta**8 + 15*beta**7 + 20*beta**6 + 15*beta**5 + 6*beta**4 + beta**3)*(x**4) + 4*(5*alpha**3*beta**2 - 5*alpha**3*beta - 4*alpha**3 + 7*alpha**2*beta**3 + 28*alpha**2*beta**2 - 19*alpha**2*beta - 16*alpha**2 + 42*alpha*beta**3 + 54*alpha*beta**2 - 28*alpha*beta - 20*alpha + 42*beta**3 + 36*beta**2 - 14*beta - 8)/(alpha**3*beta**6 + 3*alpha**3*beta**5 + 3*alpha**3*beta**4 + alpha**3*beta**3 + 3*alpha**2*beta**7 + 12*alpha**2*beta**6 + 18*alpha**2*beta**5 + 12*alpha**2*beta**4 + 3*alpha**2*beta**3 + 3*alpha*beta**8 + 15*alpha*beta**7 + 30*alpha*beta**6 + 30*alpha*beta**5 + 15*alpha*beta**4 + 3*alpha*beta**3 + beta**9 + 6*beta**8 + 15*beta**7 + 20*beta**6 + 15*beta**5 + 6*beta**4 + beta**3)*(x**3) + 3*(10*alpha**3*beta**2 + 2*alpha**3*beta - 2*alpha**3 + 14*alpha**2*beta**3 + 44*alpha**2*beta**2 + 8*alpha**2*beta - 6*alpha**2 + 42*alpha*beta**3 + 68*alpha*beta**2 + 10*alpha*beta - 6*alpha + 28*beta**3 + 34*beta**2 + 4*beta - 2)/(alpha**3*beta**6 + 3*alpha**3*beta**5 + 3*alpha**3*beta**4 + alpha**3*beta**3 + 3*alpha**2*beta**7 + 12*alpha**2*beta**6 + 18*alpha**2*beta**5 + 12*alpha**2*beta**4 + 3*alpha**2*beta**3 + 3*alpha*beta**8 + 15*alpha*beta**7 + 30*alpha*beta**6 + 30*alpha*beta**5 + 15*alpha*beta**4 + 3*alpha*beta**3 + beta**9 + 6*beta**8 + 15*beta**7 + 20*beta**6 + 15*beta**5 + 6*beta**4 + beta**3)*(x**2) + 2*(5*alpha**3*beta + 3*alpha**3 + 7*alpha**2*beta**2 + 20*alpha**2*beta + 9*alpha**2 + 14*alpha*beta**2 + 25*alpha*beta + 9*alpha + 7*beta**2 + 10*beta + 3)/(alpha**3*beta**5 + 3*alpha**3*beta**4 + 3*alpha**3*beta**3 + alpha**3*beta**2 + 3*alpha**2*beta**6 + 12*alpha**2*beta**5 + 18*alpha**2*beta**4 + 12*alpha**2*beta**3 + 3*alpha**2*beta**2 + 3*alpha*beta**7 + 15*alpha*beta**6 + 30*alpha*beta**5 + 30*alpha*beta**4 + 15*alpha*beta**3 + 3*alpha*beta**2 + beta**8 + 6*beta**7 + 15*beta**6 + 20*beta**5 + 15*beta**4 + 6*beta**3 + beta**2)*x + 0 )
def d7_prime(x, alpha, beta):
    return ( 7*1/(alpha**2*beta**4 + 2*alpha**2*beta**3 + alpha**2*beta**2 + 2*alpha*beta**5 + 6*alpha*beta**4 + 6*alpha*beta**3 + 2*alpha*beta**2 + beta**6 + 4*beta**5 + 6*beta**4 + 4*beta**3 + beta**2)*(x**6) + 6*(2*alpha - beta + 4)/(alpha**2*beta**4 + 2*alpha**2*beta**3 + alpha**2*beta**2 + 2*alpha*beta**5 + 6*alpha*beta**4 + 6*alpha*beta**3 + 2*alpha*beta**2 + beta**6 + 4*beta**5 + 6*beta**4 + 4*beta**3 + beta**2)*(x**5) + 5*(alpha**2 - 2*alpha*beta + 6*alpha - 4*beta + 6)/(alpha**2*beta**4 + 2*alpha**2*beta**3 + alpha**2*beta**2 + 2*alpha*beta**5 + 6*alpha*beta**4 + 6*alpha*beta**3 + 2*alpha*beta**2 + beta**6 + 4*beta**5 + 6*beta**4 + 4*beta**3 + beta**2)*(x**4) + 4*(-alpha**2*beta + 2*alpha**2 - 6*alpha*beta + 6*alpha - 6*beta + 4)/(alpha**2*beta**4 + 2*alpha**2*beta**3 + alpha**2*beta**2 + 2*alpha*beta**5 + 6*alpha*beta**4 + 6*alpha*beta**3 + 2*alpha*beta**2 + beta**6 + 4*beta**5 + 6*beta**4 + 4*beta**3 + beta**2)*(x**3) + 3*(-2*alpha**2*beta + alpha**2 - 6*alpha*beta + 2*alpha - 4*beta + 1)/(alpha**2*beta**4 + 2*alpha**2*beta**3 + alpha**2*beta**2 + 2*alpha*beta**5 + 6*alpha*beta**4 + 6*alpha*beta**3 + 2*alpha*beta**2 + beta**6 + 4*beta**5 + 6*beta**4 + 4*beta**3 + beta**2)*(x**2) + 2*(-alpha**2 - 2*alpha - 1)/(alpha**2*beta**3 + 2*alpha**2*beta**2 + alpha**2*beta + 2*alpha*beta**4 + 6*alpha*beta**3 + 6*alpha*beta**2 + 2*alpha*beta + beta**5 + 4*beta**4 + 6*beta**3 + 4*beta**2 + beta)*x + 0 )

#########################################################################################################


def d0_horner(x, alpha, beta) -> float:
    return ( x**2*(x*(x*(x*(x*(x*(alpha*(6*alpha + 4*beta + 8) + 2*beta + 2)/(alpha**9 + 3*alpha**8*beta + 6*alpha**8 + 3*alpha**7*beta**2 + 15*alpha**7*beta + 15*alpha**7 + alpha**6*beta**3 + 12*alpha**6*beta**2 + 30*alpha**6*beta + 20*alpha**6 + 3*alpha**5*beta**3 + 18*alpha**5*beta**2 + 30*alpha**5*beta + 15*alpha**5 + 3*alpha**4*beta**3 + 12*alpha**4*beta**2 + 15*alpha**4*beta + 6*alpha**4 + 
alpha**3*beta**3 + 3*alpha**3*beta**2 + 3*alpha**3*beta + alpha**3) + (alpha*(alpha*(7*alpha - 7*beta + 28) + beta*(-8*beta - 1) + 27) + beta*(2 - 4*beta) + 6)/(alpha**9 + 3*alpha**8*beta + 6*alpha**8 + 3*alpha**7*beta**2 + 15*alpha**7*beta + 15*alpha**7 + alpha**6*beta**3 + 12*alpha**6*beta**2 + 30*alpha**6*beta + 20*alpha**6 + 3*alpha**5*beta**3 + 18*alpha**5*beta**2 + 30*alpha**5*beta + 15*alpha**5 + 3*alpha**4*beta**3 + 12*alpha**4*beta**2 + 15*alpha**4*beta + 6*alpha**4 + alpha**3*beta**3 + 3*alpha**3*beta**2 + 3*alpha**3*beta + alpha**3)) + (alpha*(alpha*(alpha*(14 - 14*beta) + beta*(-4*beta - 46) + 38) + beta*(beta*(4*beta - 22) - 36) + 30) + beta*(beta*(2*beta - 10) - 6) + 6)/(alpha**9 + 3*alpha**8*beta + 6*alpha**8 + 3*alpha**7*beta**2 + 15*alpha**7*beta + 15*alpha**7 + alpha**6*beta**3 + 12*alpha**6*beta**2 + 30*alpha**6*beta + 20*alpha**6 + 3*alpha**5*beta**3 + 18*alpha**5*beta**2 + 30*alpha**5*beta + 15*alpha**5 + 3*alpha**4*beta**3 + 12*alpha**4*beta**2 + 15*alpha**4*beta + 6*alpha**4 + alpha**3*beta**3 + 3*alpha**3*beta**2 + 3*alpha**3*beta + alpha**3)) + (alpha*(alpha*(alpha*(beta*(7*beta - 28) + 7) + beta*(beta*(5*beta + 8) - 71) + 16) + beta*(beta*(15*beta - 9) 
- 53) + 11) + beta*(beta*(6*beta - 6) - 10) + 2)/(alpha**9 + 3*alpha**8*beta + 6*alpha**8 + 3*alpha**7*beta**2 + 15*alpha**7*beta + 15*alpha**7 + alpha**6*beta**3 + 12*alpha**6*beta**2 + 30*alpha**6*beta + 20*alpha**6 + 3*alpha**5*beta**3 + 18*alpha**5*beta**2 + 30*alpha**5*beta + 15*alpha**5 + 3*alpha**4*beta**3 + 12*alpha**4*beta**2 + 15*alpha**4*beta + 6*alpha**4 + alpha**3*beta**3 + 3*alpha**3*beta**2 + 3*alpha**3*beta + alpha**3)) + (alpha*(alpha*(alpha*beta*(14*beta - 14) + beta*(beta*(10*beta + 28) - 32)) + beta*(beta*(18*beta + 16) - 22)) + beta*(beta*(6*beta + 2) - 4))/(alpha**9 + 3*alpha**8*beta + 6*alpha**8 + 3*alpha**7*beta**2 + 15*alpha**7*beta + 15*alpha**7 + alpha**6*beta**3 + 12*alpha**6*beta**2 + 30*alpha**6*beta + 20*alpha**6 + 3*alpha**5*beta**3 + 18*alpha**5*beta**2 + 30*alpha**5*beta + 15*alpha**5 + 3*alpha**4*beta**3 + 12*alpha**4*beta**2 + 15*alpha**4*beta + 6*alpha**4 + alpha**3*beta**3 + 3*alpha**3*beta**2 + 3*alpha**3*beta + alpha**3)) + (alpha*(7*alpha*beta**2 + beta**2*(5*beta + 9)) + beta**2*(2*beta + 2))/(alpha**8 + 3*alpha**7*beta + 5*alpha**7 + 3*alpha**6*beta**2 + 12*alpha**6*beta + 10*alpha**6 + alpha**5*beta**3 + 9*alpha**5*beta**2 + 18*alpha**5*beta + 10*alpha**5 + 2*alpha**4*beta**3 + 9*alpha**4*beta**2 + 12*alpha**4*beta + 5*alpha**4 + alpha**3*beta**3 + 3*alpha**3*beta**2 + 3*alpha**3*beta + alpha**3)) )

def d1_horner(x, alpha, beta) -> float:
    return ( x**2*(beta**2/(alpha**5 + 2*alpha**4*beta + 3*alpha**4 + alpha**3*beta**2 + 4*alpha**3*beta + 3*alpha**3 + alpha**2*beta**2 + 2*alpha**2*beta + alpha**2) + x*(x*(x*(x*(x/(alpha**6 + 2*alpha**5*beta + 4*alpha**5 + alpha**4*beta**2 + 6*alpha**4*beta + 6*alpha**4 + 2*alpha**3*beta**2 + 6*alpha**3*beta + 4*alpha**3 + alpha**2*beta**2 + 2*alpha**2*beta + alpha**2) + (alpha - 2*beta + 3)/(alpha**6 + 2*alpha**5*beta + 4*alpha**5 + alpha**4*beta**2 + 6*alpha**4*beta + 6*alpha**4 + 2*alpha**3*beta**2 + 6*alpha**3*beta + 4*alpha**3 + alpha**2*beta**2 + 2*alpha**2*beta + alpha**2)) + (alpha*(2 - 2*beta) + beta*(beta - 6) + 3)/(alpha**6 + 2*alpha**5*beta + 4*alpha**5 + alpha**4*beta**2 + 6*alpha**4*beta + 6*alpha**4 + 2*alpha**3*beta**2 + 6*alpha**3*beta + 4*alpha**3 + alpha**2*beta**2 
+ 2*alpha**2*beta + alpha**2)) + (alpha*(beta*(beta - 4) + 1) + beta*(3*beta - 6) + 1)/(alpha**6 + 2*alpha**5*beta + 4*alpha**5 + alpha**4*beta**2 + 6*alpha**4*beta + 6*alpha**4 + 2*alpha**3*beta**2 + 6*alpha**3*beta + 4*alpha**3 + alpha**2*beta**2 + 2*alpha**2*beta + alpha**2)) + (alpha*beta*(2*beta - 2) + beta*(3*beta - 2))/(alpha**6 + 2*alpha**5*beta + 4*alpha**5 + alpha**4*beta**2 + 6*alpha**4*beta + 6*alpha**4 + 2*alpha**3*beta**2 + 6*alpha**3*beta + 4*alpha**3 + alpha**2*beta**2 + 2*alpha**2*beta + alpha**2))) )

def d2_horner(x, alpha, beta) -> float:
    return ( x**2*(x*(x*(x*(x*(x*(alpha*(2*beta + 4) - 2*beta - 2)/(alpha**3*beta**3 + 3*alpha**3*beta**2 + 3*alpha**3*beta + alpha**3) + (alpha*(alpha*(4*beta + 8) + beta*(-4*beta - 5) + 9) + beta*(4*beta - 2) - 6)/(alpha**3*beta**3 + 3*alpha**3*beta**2 + 3*alpha**3*beta + alpha**3)) + (alpha*(alpha*(alpha*(2*beta + 4) + beta*(-8*beta - 8) + 16) + beta*(beta*(2*beta - 2) - 18) + 6) + beta*(beta*(10 - 2*beta) + 6) - 6)/(alpha**3*beta**3 + 3*alpha**3*beta**2 + 3*alpha**3*beta + alpha**3)) + (alpha*(alpha*(alpha*(beta*(-4*beta - 5) + 5) + beta*(beta*(4*beta - 8) - 28) + 8) + beta*(beta*(3*beta + 9) - 13) + 1) + beta*(beta*(6 - 6*beta) + 10) - 2)/(alpha**3*beta**3 + 3*alpha**3*beta**2 + 3*alpha**3*beta + alpha**3)) + (alpha*(alpha*(alpha*beta*(beta*(2*beta - 2) - 10) + beta*(beta*(8*beta + 
8) - 16)) + beta*(8*beta - 2)) + beta*(beta*(-6*beta - 2) + 4))/(alpha**3*beta**3 + 3*alpha**3*beta**2 + 3*alpha**3*beta + alpha**3)) + (alpha*(alpha*(alpha*beta**2*(3*beta + 5) + beta**2*(4*beta + 8)) + beta**2*(1 - beta)) + beta**2*(-2*beta - 2))/(alpha**3*beta**3 + 3*alpha**3*beta**2 + 3*alpha**3*beta + alpha**3)) )

def d3_horner(x, alpha, beta) -> float:
    return ( x**2*(x*(x*(x*(x*(x/(alpha**2*beta**2 + 2*alpha**2*beta + alpha**2) + (2*alpha - 2*beta + 3)/(alpha**2*beta**2 + 2*alpha**2*beta + alpha**2)) + (alpha*(alpha - 4*beta + 4) + beta*(beta - 6) + 3)/(alpha**2*beta**2 + 2*alpha**2*beta + alpha**2)) + (alpha*(alpha*(1 - 2*beta) + beta*(2*beta - 8) + 2) + beta*(3*beta - 6) + 1)/(alpha**2*beta**2 + 2*alpha**2*beta + alpha**2)) + (alpha*(alpha*beta*(beta - 2) + beta*(4*beta - 4)) + beta*(3*beta - 2))/(alpha**2*beta**2 + 2*alpha**2*beta + alpha**2)) + (alpha*(alpha*beta**2 
+ 2*beta**2) + beta**2)/(alpha**2*beta**2 + 2*alpha**2*beta + alpha**2)) )

def d4_horner(x, alpha, beta) -> float:
    return ( x**2*(x*(x*(x*(x*(x*(alpha*(2 - 2*beta) - 4*beta + 2)/(alpha**3*beta**3 + 3*alpha**2*beta**3 + 3*alpha*beta**3 + beta**3) + (alpha*(alpha*(4 - 4*beta) + beta*(4*beta - 19) + 12) + beta*(8*beta - 19) + 8)/(alpha**3*beta**3 + 3*alpha**2*beta**3 + 3*alpha*beta**3 + beta**3)) + (alpha*(alpha*(alpha*(2 - 2*beta) + beta*(8*beta - 22) + 14) + beta*(beta*(32 - 2*beta) - 54) + 24) + beta*(beta*(32 - 4*beta) - 36) + 12)/(alpha**3*beta**3 + 3*alpha**2*beta**3 + 3*alpha*beta**3 + beta**3)) + (alpha*(alpha*(alpha*(beta*(4*beta - 7) + 4) + beta*(beta*(32 - 4*beta) - 41) + 16) + beta*(beta*(72 - 15*beta) - 68) + 20) + beta*(beta*(48 - 15*beta) - 34) + 8)/(alpha**3*beta**3 + 3*alpha**2*beta**3 + 3*alpha*beta**3 + beta**3)) + (alpha*(alpha*(alpha*(beta*(beta*(8 - 2*beta) - 8) + 2) + beta*(beta*(40 - 14*beta) - 32) + 6) + beta*(beta*(64 - 30*beta) - 40) + 6) + beta*(beta*(32 - 20*beta) - 16) + 2)/(alpha**3*beta**3 + 3*alpha**2*beta**3 + 3*alpha*beta**3 + beta**3)) + (alpha*(alpha*(beta*(4 - 3*beta) - 3) + beta*(12 - 10*beta) - 6) + beta*(8 - 10*beta) - 
3)/(alpha**2*beta**2 + 2*alpha*beta**2 + beta**2)) + 1 )

def d5_horner(x, alpha, beta) -> float:
    return ( x*(x*(x*(x*(x*(x*(x/(alpha**2*beta**2 + 2*alpha*beta**2 + beta**2) + (2*alpha - 2*beta + 4)/(alpha**2*beta**2 + 2*alpha*beta**2 + beta**2)) + (alpha*(alpha - 4*beta + 6) + beta*(beta - 8) + 6)/(alpha**2*beta**2 + 2*alpha*beta**2 + beta**2)) + (alpha*(alpha*(2 - 2*beta) + beta*(2*beta - 12) + 6) + beta*(4*beta - 12) + 4)/(alpha**2*beta**2 + 2*alpha*beta**2 + beta**2)) + (alpha*(alpha*(beta*(beta - 4) + 1) + beta*(6*beta - 12) + 2) + beta*(6*beta - 8) + 1)/(alpha**2*beta**2 + 2*alpha*beta**2 + beta**2)) + (alpha*(2*beta - 2) + 4*beta - 2)/(alpha*beta + beta)) + 1) )

def d6_horner(x, alpha, beta) -> float:
    return ( x**2*(x*(x*(x*(x*(x*(alpha*(-4*beta - 2) + beta*(-6*beta - 8) - 2)/(alpha**3*beta**6 + 3*alpha**3*beta**5 + 3*alpha**3*beta**4 + alpha**3*beta**3 + 3*alpha**2*beta**7 + 12*alpha**2*beta**6 + 18*alpha**2*beta**5 + 12*alpha**2*beta**4 + 3*alpha**2*beta**3 + 3*alpha*beta**8 + 15*alpha*beta**7 + 30*alpha*beta**6 + 30*alpha*beta**5 + 15*alpha*beta**4 + 3*alpha*beta**3 + beta**9 + 6*beta**8 + 15*beta**7 + 20*beta**6 + 15*beta**5 + 6*beta**4 + beta**3) + (alpha*(alpha*(-8*beta - 4) + beta*(-7*beta - 29) - 12) + beta*(beta*(7*beta - 14) - 29) - 8)/(alpha**3*beta**6 + 3*alpha**3*beta**5 + 3*alpha**3*beta**4 + alpha**3*beta**3 + 3*alpha**2*beta**7 + 12*alpha**2*beta**6 + 18*alpha**2*beta**5 + 12*alpha**2*beta**4 + 3*alpha**2*beta**3 + 3*alpha*beta**8 + 15*alpha*beta**7 + 30*alpha*beta**6 + 30*alpha*beta**5 + 15*alpha*beta**4 + 3*alpha*beta**3 + beta**9 + 6*beta**8 + 15*beta**7 + 20*beta**6 + 15*beta**5 + 6*beta**4 
+ beta**3)) + (alpha*(alpha*(alpha*(-4*beta - 2) + beta*(4*beta - 26) - 14) + beta*(beta*(14*beta + 4) - 54) - 24) + beta*(beta*(28*beta + 4) - 36) - 12)/(alpha**3*beta**6 + 3*alpha**3*beta**5 + 3*alpha**3*beta**4 + alpha**3*beta**3 + 3*alpha**2*beta**7 + 12*alpha**2*beta**6 + 18*alpha**2*beta**5 + 12*alpha**2*beta**4 + 3*alpha**2*beta**3 + 3*alpha*beta**8 + 15*alpha*beta**7 + 30*alpha*beta**6 + 
30*alpha*beta**5 + 15*alpha*beta**4 + 3*alpha*beta**3 + beta**9 + 6*beta**8 + 15*beta**7 + 20*beta**6 + 15*beta**5 + 6*beta**4 + beta**3)) + (alpha*(alpha*(alpha*(beta*(5*beta - 5) - 4) + beta*(beta*(7*beta + 28) - 19) - 16) + beta*(beta*(42*beta + 54) - 28) - 20) + beta*(beta*(42*beta + 36) - 14) - 8)/(alpha**3*beta**6 + 3*alpha**3*beta**5 + 3*alpha**3*beta**4 + alpha**3*beta**3 + 3*alpha**2*beta**7 + 12*alpha**2*beta**6 + 18*alpha**2*beta**5 + 12*alpha**2*beta**4 + 3*alpha**2*beta**3 + 3*alpha*beta**8 + 15*alpha*beta**7 + 30*alpha*beta**6 + 30*alpha*beta**5 + 15*alpha*beta**4 + 3*alpha*beta**3 + beta**9 + 6*beta**8 + 15*beta**7 + 20*beta**6 + 15*beta**5 + 6*beta**4 + beta**3)) + (alpha*(alpha*(alpha*(beta*(10*beta + 2) - 2) + beta*(beta*(14*beta + 44) + 8) - 6) + beta*(beta*(42*beta + 
68) + 10) - 6) + beta*(beta*(28*beta + 34) + 4) - 2)/(alpha**3*beta**6 + 3*alpha**3*beta**5 + 3*alpha**3*beta**4 + alpha**3*beta**3 + 3*alpha**2*beta**7 + 12*alpha**2*beta**6 + 18*alpha**2*beta**5 + 12*alpha**2*beta**4 + 3*alpha**2*beta**3 + 3*alpha*beta**8 + 15*alpha*beta**7 + 30*alpha*beta**6 + 30*alpha*beta**5 + 15*alpha*beta**4 + 3*alpha*beta**3 + beta**9 + 6*beta**8 + 15*beta**7 + 20*beta**6 + 15*beta**5 + 6*beta**4 + beta**3)) + (alpha*(alpha*(alpha*(5*beta + 3) + beta*(7*beta + 20) + 9) + beta*(14*beta + 25) + 9) + beta*(7*beta + 10) + 3)/(alpha**3*beta**5 + 3*alpha**3*beta**4 + 3*alpha**3*beta**3 + alpha**3*beta**2 + 3*alpha**2*beta**6 + 12*alpha**2*beta**5 + 18*alpha**2*beta**4 + 12*alpha**2*beta**3 + 3*alpha**2*beta**2 + 3*alpha*beta**7 + 15*alpha*beta**6 + 30*alpha*beta**5 + 30*alpha*beta**4 + 15*alpha*beta**3 + 3*alpha*beta**2 + beta**8 + 6*beta**7 + 15*beta**6 + 20*beta**5 + 15*beta**4 + 6*beta**3 + beta**2)) )

def d7_horner(x, alpha, beta) -> float:
    return ( x**2*(x*(x*(x*(x*(x/(alpha**2*beta**4 + 2*alpha**2*beta**3 + alpha**2*beta**2 + 2*alpha*beta**5 + 6*alpha*beta**4 + 6*alpha*beta**3 + 2*alpha*beta**2 + beta**6 + 4*beta**5 + 6*beta**4 + 4*beta**3 + beta**2) + (2*alpha - beta + 4)/(alpha**2*beta**4 + 2*alpha**2*beta**3 + alpha**2*beta**2 + 2*alpha*beta**5 + 6*alpha*beta**4 + 6*alpha*beta**3 + 2*alpha*beta**2 + beta**6 + 4*beta**5 + 6*beta**4 + 4*beta**3 + beta**2)) + (alpha*(alpha - 2*beta + 6) - 4*beta + 6)/(alpha**2*beta**4 + 2*alpha**2*beta**3 + alpha**2*beta**2 + 2*alpha*beta**5 + 6*alpha*beta**4 + 6*alpha*beta**3 + 2*alpha*beta**2 + beta**6 + 4*beta**5 + 6*beta**4 + 4*beta**3 + beta**2)) + 
(alpha*(alpha*(2 - beta) - 6*beta + 6) - 6*beta + 4)/(alpha**2*beta**4 + 2*alpha**2*beta**3 + alpha**2*beta**2 + 2*alpha*beta**5 + 6*alpha*beta**4 + 6*alpha*beta**3 + 2*alpha*beta**2 + beta**6 + 4*beta**5 + 6*beta**4 + 4*beta**3 + beta**2)) + (alpha*(alpha*(1 - 2*beta) - 6*beta + 2) - 4*beta + 1)/(alpha**2*beta**4 + 2*alpha**2*beta**3 + alpha**2*beta**2 + 2*alpha*beta**5 + 6*alpha*beta**4 + 6*alpha*beta**3 + 2*alpha*beta**2 + beta**6 + 4*beta**5 + 6*beta**4 + 4*beta**3 + beta**2)) + (alpha*(-alpha - 2) - 1)/(alpha**2*beta**3 + 2*alpha**2*beta**2 + alpha**2*beta + 2*alpha*beta**4 + 6*alpha*beta**3 + 6*alpha*beta**2 + 2*alpha*beta + beta**5 + 4*beta**4 + 6*beta**3 + 4*beta**2 + beta)) )

def d0_prime_horner(x, alpha, beta) -> float:
    return ( x*(x*(x*(x*(x*(x*(alpha*(42*alpha + 28*beta + 56) + 14*beta + 14)/(alpha**9 + 3*alpha**8*beta + 6*alpha**8 + 3*alpha**7*beta**2 + 15*alpha**7*beta + 15*alpha**7 + alpha**6*beta**3 + 12*alpha**6*beta**2 + 30*alpha**6*beta + 20*alpha**6 + 3*alpha**5*beta**3 + 18*alpha**5*beta**2 + 30*alpha**5*beta + 15*alpha**5 + 3*alpha**4*beta**3 + 12*alpha**4*beta**2 + 15*alpha**4*beta + 6*alpha**4 
+ alpha**3*beta**3 + 3*alpha**3*beta**2 + 3*alpha**3*beta + alpha**3) + (alpha*(alpha*(42*alpha - 42*beta + 168) + beta*(-48*beta - 6) + 162) + beta*(12 - 24*beta) + 36)/(alpha**9 + 3*alpha**8*beta + 6*alpha**8 + 3*alpha**7*beta**2 + 15*alpha**7*beta + 15*alpha**7 + alpha**6*beta**3 + 12*alpha**6*beta**2 + 30*alpha**6*beta + 20*alpha**6 + 3*alpha**5*beta**3 + 18*alpha**5*beta**2 + 30*alpha**5*beta + 15*alpha**5 + 3*alpha**4*beta**3 + 12*alpha**4*beta**2 + 15*alpha**4*beta + 6*alpha**4 + alpha**3*beta**3 + 3*alpha**3*beta**2 + 
3*alpha**3*beta + alpha**3)) + (alpha*(alpha*(alpha*(70 - 70*beta) + beta*(-20*beta - 230) + 190) + beta*(beta*(20*beta - 110) - 180) + 150) + beta*(beta*(10*beta - 50) - 30) + 30)/(alpha**9 + 3*alpha**8*beta + 6*alpha**8 + 3*alpha**7*beta**2 + 15*alpha**7*beta + 15*alpha**7 + alpha**6*beta**3 + 12*alpha**6*beta**2 + 30*alpha**6*beta + 20*alpha**6 + 3*alpha**5*beta**3 + 18*alpha**5*beta**2 + 30*alpha**5*beta + 15*alpha**5 + 3*alpha**4*beta**3 + 12*alpha**4*beta**2 + 15*alpha**4*beta + 6*alpha**4 + alpha**3*beta**3 + 3*alpha**3*beta**2 + 3*alpha**3*beta + alpha**3)) + (alpha*(alpha*(alpha*(beta*(28*beta - 112) + 28) + beta*(beta*(20*beta + 32) - 284) + 64) + beta*(beta*(60*beta - 36) - 212) + 44) + beta*(beta*(24*beta - 24) - 40) + 8)/(alpha**9 + 3*alpha**8*beta + 6*alpha**8 + 3*alpha**7*beta**2 + 15*alpha**7*beta + 15*alpha**7 + alpha**6*beta**3 + 12*alpha**6*beta**2 + 30*alpha**6*beta + 20*alpha**6 + 3*alpha**5*beta**3 + 18*alpha**5*beta**2 + 30*alpha**5*beta + 15*alpha**5 + 3*alpha**4*beta**3 + 12*alpha**4*beta**2 + 15*alpha**4*beta + 6*alpha**4 
+ alpha**3*beta**3 + 3*alpha**3*beta**2 + 3*alpha**3*beta + alpha**3)) + (alpha*(alpha*(alpha*beta*(42*beta - 42) + beta*(beta*(30*beta + 84) - 96)) + beta*(beta*(54*beta + 48) - 66)) + beta*(beta*(18*beta + 6) - 12))/(alpha**9 + 3*alpha**8*beta + 6*alpha**8 + 3*alpha**7*beta**2 + 15*alpha**7*beta + 15*alpha**7 + alpha**6*beta**3 + 12*alpha**6*beta**2 + 30*alpha**6*beta + 20*alpha**6 + 3*alpha**5*beta**3 + 18*alpha**5*beta**2 + 30*alpha**5*beta + 15*alpha**5 + 3*alpha**4*beta**3 + 12*alpha**4*beta**2 + 15*alpha**4*beta + 6*alpha**4 + alpha**3*beta**3 + 3*alpha**3*beta**2 + 3*alpha**3*beta + alpha**3)) + (alpha*(14*alpha*beta**2 + beta**2*(10*beta + 18)) + beta**2*(4*beta + 4))/(alpha**8 + 3*alpha**7*beta + 5*alpha**7 + 3*alpha**6*beta**2 + 12*alpha**6*beta + 10*alpha**6 + alpha**5*beta**3 + 9*alpha**5*beta**2 + 18*alpha**5*beta + 10*alpha**5 + 2*alpha**4*beta**3 + 9*alpha**4*beta**2 + 12*alpha**4*beta + 5*alpha**4 + alpha**3*beta**3 + 3*alpha**3*beta**2 + 3*alpha**3*beta + alpha**3)) )

def d1_prime_horner(x, alpha, beta) -> float:
    return ( x*(2*beta**2/(alpha**5 + 2*alpha**4*beta + 3*alpha**4 + alpha**3*beta**2 + 4*alpha**3*beta + 3*alpha**3 + alpha**2*beta**2 + 2*alpha**2*beta + alpha**2) + x*(x*(x*(x*(7*x/(alpha**6 + 2*alpha**5*beta + 4*alpha**5 + alpha**4*beta**2 + 6*alpha**4*beta + 6*alpha**4 + 2*alpha**3*beta**2 + 6*alpha**3*beta + 4*alpha**3 + alpha**2*beta**2 + 2*alpha**2*beta + alpha**2) + (6*alpha - 12*beta + 
18)/(alpha**6 + 2*alpha**5*beta + 4*alpha**5 + alpha**4*beta**2 + 6*alpha**4*beta + 6*alpha**4 + 2*alpha**3*beta**2 + 6*alpha**3*beta + 4*alpha**3 + alpha**2*beta**2 + 2*alpha**2*beta + alpha**2)) + (alpha*(10 - 10*beta) + beta*(5*beta - 30) + 15)/(alpha**6 + 2*alpha**5*beta + 4*alpha**5 + alpha**4*beta**2 + 6*alpha**4*beta + 6*alpha**4 + 2*alpha**3*beta**2 + 6*alpha**3*beta + 4*alpha**3 + alpha**2*beta**2 + 2*alpha**2*beta + alpha**2)) + (alpha*(beta*(4*beta - 16) + 4) + beta*(12*beta - 24) + 4)/(alpha**6 + 2*alpha**5*beta + 
4*alpha**5 + alpha**4*beta**2 + 6*alpha**4*beta + 6*alpha**4 + 2*alpha**3*beta**2 + 6*alpha**3*beta + 4*alpha**3 + alpha**2*beta**2 + 2*alpha**2*beta + alpha**2)) + (alpha*beta*(6*beta - 6) + beta*(9*beta - 6))/(alpha**6 + 2*alpha**5*beta + 4*alpha**5 + alpha**4*beta**2 + 6*alpha**4*beta + 6*alpha**4 + 2*alpha**3*beta**2 + 6*alpha**3*beta + 4*alpha**3 + alpha**2*beta**2 + 2*alpha**2*beta + alpha**2))) )

def d2_prime_horner(x, alpha, beta) -> float:
    return ( x*(x*(x*(x*(x*(x*(alpha*(14*beta + 28) - 14*beta - 14)/(alpha**3*beta**3 + 3*alpha**3*beta**2 + 3*alpha**3*beta + alpha**3) + (alpha*(alpha*(24*beta + 48) + beta*(-24*beta - 30) + 54) + beta*(24*beta - 12) - 36)/(alpha**3*beta**3 + 3*alpha**3*beta**2 + 
3*alpha**3*beta + alpha**3)) + (alpha*(alpha*(alpha*(10*beta + 20) + beta*(-40*beta - 40) + 80) + beta*(beta*(10*beta - 10) - 90) + 30) + beta*(beta*(50 - 10*beta) + 30) - 30)/(alpha**3*beta**3 + 3*alpha**3*beta**2 + 3*alpha**3*beta + alpha**3)) + (alpha*(alpha*(alpha*(beta*(-16*beta - 20) + 20) + beta*(beta*(16*beta - 32) - 112) + 32) + beta*(beta*(12*beta + 36) - 52) + 4) + beta*(beta*(24 - 24*beta) + 40) - 8)/(alpha**3*beta**3 + 3*alpha**3*beta**2 + 3*alpha**3*beta + alpha**3)) + (alpha*(alpha*(alpha*beta*(beta*(6*beta - 6) - 30) + beta*(beta*(24*beta + 24) - 48)) + beta*(24*beta - 6)) + beta*(beta*(-18*beta - 6) + 12))/(alpha**3*beta**3 + 3*alpha**3*beta**2 + 3*alpha**3*beta + alpha**3)) + (alpha*(alpha*(alpha*beta**2*(6*beta + 10) + beta**2*(8*beta + 16)) + beta**2*(2 - 2*beta)) + beta**2*(-4*beta - 4))/(alpha**3*beta**3 + 3*alpha**3*beta**2 + 3*alpha**3*beta + alpha**3)) )

def d3_prime_horner(x, alpha, beta) -> float:
    return ( x*(x*(x*(x*(x*(7*x/(alpha**2*beta**2 + 2*alpha**2*beta + alpha**2) + (12*alpha - 12*beta + 18)/(alpha**2*beta**2 + 2*alpha**2*beta + alpha**2)) + (alpha*(5*alpha - 20*beta + 20) + beta*(5*beta - 30) + 15)/(alpha**2*beta**2 + 2*alpha**2*beta + alpha**2)) + (alpha*(alpha*(4 - 8*beta) + beta*(8*beta - 32) + 8) + beta*(12*beta - 24) + 4)/(alpha**2*beta**2 + 2*alpha**2*beta + alpha**2)) + (alpha*(alpha*beta*(3*beta - 6) + beta*(12*beta - 12)) + beta*(9*beta - 6))/(alpha**2*beta**2 + 2*alpha**2*beta + alpha**2)) + (alpha*(2*alpha*beta**2 + 4*beta**2) + 2*beta**2)/(alpha**2*beta**2 + 2*alpha**2*beta + alpha**2)) )

def d4_prime_horner(x, alpha, beta) -> float:
    return ( x*(x*(x*(x*(x*(x*(alpha*(14 - 14*beta) - 28*beta + 14)/(alpha**3*beta**3 + 3*alpha**2*beta**3 + 3*alpha*beta**3 + beta**3) + (alpha*(alpha*(24 - 24*beta) + beta*(24*beta - 114) + 72) + beta*(48*beta - 114) + 48)/(alpha**3*beta**3 + 3*alpha**2*beta**3 + 
3*alpha*beta**3 + beta**3)) + (alpha*(alpha*(alpha*(10 - 10*beta) + beta*(40*beta - 110) + 70) + beta*(beta*(160 - 10*beta) - 270) + 
120) + beta*(beta*(160 - 20*beta) - 180) + 60)/(alpha**3*beta**3 + 3*alpha**2*beta**3 + 3*alpha*beta**3 + beta**3)) + (alpha*(alpha*(alpha*(beta*(16*beta - 28) + 16) + beta*(beta*(128 - 16*beta) - 164) + 64) + beta*(beta*(288 - 60*beta) - 272) + 80) + beta*(beta*(192 - 60*beta) - 136) + 32)/(alpha**3*beta**3 + 3*alpha**2*beta**3 + 3*alpha*beta**3 + beta**3)) + (alpha*(alpha*(alpha*(beta*(beta*(24 - 6*beta) - 24) + 6) + beta*(beta*(120 - 42*beta) - 96) + 18) + beta*(beta*(192 - 90*beta) - 120) + 18) + beta*(beta*(96 - 60*beta) 
- 48) + 6)/(alpha**3*beta**3 + 3*alpha**2*beta**3 + 3*alpha*beta**3 + beta**3)) + (alpha*(alpha*(beta*(8 - 6*beta) - 6) + beta*(24 - 
20*beta) - 12) + beta*(16 - 20*beta) - 6)/(alpha**2*beta**2 + 2*alpha*beta**2 + beta**2)) )

def d5_prime_horner(x, alpha, beta) -> float:
    return ( x*(x*(x*(x*(x*(7*x/(alpha**2*beta**2 + 2*alpha*beta**2 + beta**2) + (12*alpha - 12*beta + 24)/(alpha**2*beta**2 + 2*alpha*beta**2 + beta**2)) + (alpha*(5*alpha - 20*beta + 30) + beta*(5*beta - 40) + 30)/(alpha**2*beta**2 + 2*alpha*beta**2 + beta**2)) + 
(alpha*(alpha*(8 - 8*beta) + beta*(8*beta - 48) + 24) + beta*(16*beta - 48) + 16)/(alpha**2*beta**2 + 2*alpha*beta**2 + beta**2)) + (alpha*(alpha*(beta*(3*beta - 12) + 3) + beta*(18*beta - 36) + 6) + beta*(18*beta - 24) + 3)/(alpha**2*beta**2 + 2*alpha*beta**2 + beta**2)) + (alpha*(4*beta - 4) + 8*beta - 4)/(alpha*beta + beta)) + 1 )

def d6_prime_horner(x, alpha, beta) -> float:
    return ( x*(x*(x*(x*(x*(x*(alpha*(-28*beta - 14) + beta*(-42*beta - 56) - 14)/(alpha**3*beta**6 + 3*alpha**3*beta**5 + 3*alpha**3*beta**4 + alpha**3*beta**3 + 3*alpha**2*beta**7 + 12*alpha**2*beta**6 + 18*alpha**2*beta**5 + 12*alpha**2*beta**4 + 3*alpha**2*beta**3 + 3*alpha*beta**8 + 15*alpha*beta**7 + 30*alpha*beta**6 + 30*alpha*beta**5 + 15*alpha*beta**4 + 3*alpha*beta**3 + beta**9 + 6*beta**8 + 15*beta**7 + 20*beta**6 + 15*beta**5 + 6*beta**4 + beta**3) + (alpha*(alpha*(-48*beta - 24) + beta*(-42*beta - 174) - 72) + beta*(beta*(42*beta - 84) - 174) - 48)/(alpha**3*beta**6 + 3*alpha**3*beta**5 + 3*alpha**3*beta**4 + alpha**3*beta**3 + 3*alpha**2*beta**7 + 12*alpha**2*beta**6 + 18*alpha**2*beta**5 + 12*alpha**2*beta**4 + 3*alpha**2*beta**3 + 3*alpha*beta**8 + 15*alpha*beta**7 + 30*alpha*beta**6 + 30*alpha*beta**5 + 15*alpha*beta**4 + 3*alpha*beta**3 + beta**9 + 6*beta**8 + 15*beta**7 + 20*beta**6 + 15*beta**5 + 6*beta**4 + beta**3)) + (alpha*(alpha*(alpha*(-20*beta - 10) + beta*(20*beta - 130) - 70) + beta*(beta*(70*beta + 20) - 270) - 120) + 
beta*(beta*(140*beta + 20) - 180) - 60)/(alpha**3*beta**6 + 3*alpha**3*beta**5 + 3*alpha**3*beta**4 + alpha**3*beta**3 + 3*alpha**2*beta**7 + 12*alpha**2*beta**6 + 18*alpha**2*beta**5 + 12*alpha**2*beta**4 + 3*alpha**2*beta**3 + 3*alpha*beta**8 + 15*alpha*beta**7 + 
30*alpha*beta**6 + 30*alpha*beta**5 + 15*alpha*beta**4 + 3*alpha*beta**3 + beta**9 + 6*beta**8 + 15*beta**7 + 20*beta**6 + 15*beta**5 + 6*beta**4 + beta**3)) + (alpha*(alpha*(alpha*(beta*(20*beta - 20) - 16) + beta*(beta*(28*beta + 112) - 76) - 64) + beta*(beta*(168*beta + 216) - 112) - 80) + beta*(beta*(168*beta + 144) - 56) - 32)/(alpha**3*beta**6 + 3*alpha**3*beta**5 + 3*alpha**3*beta**4 + alpha**3*beta**3 + 3*alpha**2*beta**7 + 12*alpha**2*beta**6 + 18*alpha**2*beta**5 + 12*alpha**2*beta**4 + 3*alpha**2*beta**3 + 3*alpha*beta**8 + 15*alpha*beta**7 + 30*alpha*beta**6 + 30*alpha*beta**5 + 15*alpha*beta**4 + 3*alpha*beta**3 + beta**9 + 6*beta**8 + 15*beta**7 + 20*beta**6 + 15*beta**5 + 6*beta**4 + beta**3)) + (alpha*(alpha*(alpha*(beta*(30*beta + 6) - 6) + beta*(beta*(42*beta + 132) + 24) - 18) + beta*(beta*(126*beta + 204) + 30) - 18) + beta*(beta*(84*beta + 102) + 12) - 6)/(alpha**3*beta**6 + 3*alpha**3*beta**5 + 3*alpha**3*beta**4 + alpha**3*beta**3 + 3*alpha**2*beta**7 + 12*alpha**2*beta**6 + 18*alpha**2*beta**5 + 12*alpha**2*beta**4 + 3*alpha**2*beta**3 + 3*alpha*beta**8 + 15*alpha*beta**7 + 30*alpha*beta**6 + 30*alpha*beta**5 + 15*alpha*beta**4 + 3*alpha*beta**3 + beta**9 + 6*beta**8 + 15*beta**7 + 20*beta**6 + 15*beta**5 + 6*beta**4 + beta**3)) + (alpha*(alpha*(alpha*(10*beta + 6) + beta*(14*beta + 40) + 18) + beta*(28*beta + 50) + 18) + beta*(14*beta + 20) + 6)/(alpha**3*beta**5 + 3*alpha**3*beta**4 + 3*alpha**3*beta**3 + alpha**3*beta**2 + 3*alpha**2*beta**6 + 12*alpha**2*beta**5 + 18*alpha**2*beta**4 + 12*alpha**2*beta**3 + 3*alpha**2*beta**2 + 3*alpha*beta**7 + 15*alpha*beta**6 + 30*alpha*beta**5 + 30*alpha*beta**4 + 15*alpha*beta**3 + 3*alpha*beta**2 + beta**8 + 6*beta**7 + 15*beta**6 + 
20*beta**5 + 15*beta**4 + 6*beta**3 + beta**2)) )

def d7_prime_horner(x, alpha, beta) -> float:
    return ( x*(x*(x*(x*(x*(7*x/(alpha**2*beta**4 + 2*alpha**2*beta**3 + alpha**2*beta**2 + 2*alpha*beta**5 + 6*alpha*beta**4 + 6*alpha*beta**3 + 2*alpha*beta**2 + beta**6 + 4*beta**5 + 6*beta**4 + 4*beta**3 + beta**2) + (12*alpha - 6*beta + 24)/(alpha**2*beta**4 + 
2*alpha**2*beta**3 + alpha**2*beta**2 + 2*alpha*beta**5 + 6*alpha*beta**4 + 6*alpha*beta**3 + 2*alpha*beta**2 + beta**6 + 4*beta**5 + 6*beta**4 + 4*beta**3 + beta**2)) + (alpha*(5*alpha - 10*beta + 30) - 20*beta + 30)/(alpha**2*beta**4 + 2*alpha**2*beta**3 + alpha**2*beta**2 + 2*alpha*beta**5 + 6*alpha*beta**4 + 6*alpha*beta**3 + 2*alpha*beta**2 + beta**6 + 4*beta**5 + 6*beta**4 + 4*beta**3 + beta**2)) + (alpha*(alpha*(8 - 4*beta) - 24*beta + 24) - 24*beta + 16)/(alpha**2*beta**4 + 2*alpha**2*beta**3 + alpha**2*beta**2 + 2*alpha*beta**5 + 6*alpha*beta**4 + 6*alpha*beta**3 + 2*alpha*beta**2 + beta**6 + 4*beta**5 + 6*beta**4 + 4*beta**3 + beta**2)) + (alpha*(alpha*(3 - 6*beta) - 18*beta + 6) - 12*beta + 3)/(alpha**2*beta**4 + 2*alpha**2*beta**3 + alpha**2*beta**2 + 2*alpha*beta**5 + 6*alpha*beta**4 + 6*alpha*beta**3 + 2*alpha*beta**2 + beta**6 + 4*beta**5 + 6*beta**4 + 4*beta**3 + beta**2)) + (alpha*(-2*alpha - 4) - 2)/(alpha**2*beta**3 + 2*alpha**2*beta**2 + alpha**2*beta + 2*alpha*beta**4 + 6*alpha*beta**3 + 6*alpha*beta**2 + 2*alpha*beta + beta**5 + 4*beta**4 + 6*beta**3 + 4*beta**2 + beta)) )


#########################################################################################################

# I use a class to represent the Hermite Birkhoff interpolant
# we will have an instance of this class on each step
class HB8:
    def __init__(   
        self, 
        x_i_minus_2, x_i_minus_1, x_i, x_i_plus_1,
        y_i_minus_2, f_i_minus_2,
        y_i_minus_1, f_i_minus_1,
        y_i, f_i,
        y_i_plus_1, f_i_plus_1,
        monitor
    ):
        h_i = x_i - x_i_minus_1
        h_i_minus_1 = x_i_minus_1 - x_i_minus_2
        h_i_plus_1 = x_i_plus_1 - x_i
        
        self.alpha = h_i_minus_1 / h_i
        self.beta = h_i_plus_1 / h_i

        monitor.different_values_alpha.add(self.alpha)
        monitor.different_values_beta.add(self.beta)

        self.h_i = h_i
        self.x_i = x_i

        # we also store x_i_minus_1 and x_i_plus_1 so that we can build the final interpolant
        self.x_i_plus_1 = x_i_plus_1
        self.x_i_minus_1 = x_i_minus_1
        self.x_i_minus_2 = x_i_minus_2

        self.y_i_minus_2 = y_i_minus_2 
        self.f_i_minus_2 = f_i_minus_2

        self.y_i_minus_1 = y_i_minus_1 
        self.f_i_minus_1 = f_i_minus_1
        
        self.y_i = y_i
        self.f_i = f_i

        self.y_i_plus_1 = y_i_plus_1
        self.f_i_plus_1 = f_i_plus_1

        xs = get_Chebyshev_nodes(x_i_minus_2, x_i_plus_1, 8)
        ys = [self.eval(x) for x in xs]
        self.eval_bary_interp = BarycentricInterpolator(xs, ys)

        y_primes = [self.prime(x) for x in xs]
        self.prime_bary_interp = BarycentricInterpolator(xs, y_primes) 

    def eval(self, x):
        pheta = (x - self.x_i) / self.h_i  # x = t_i + pheta*h_i so pheta = (x - t_i) / h_i
        return (  
                             d0(pheta, self.alpha, self.beta) * self.y_i_minus_2 
                + self.h_i * d1(pheta, self.alpha, self.beta) * self.f_i_minus_2

                           + d2(pheta, self.alpha, self.beta) * self.y_i_minus_1 
                + self.h_i * d3(pheta, self.alpha, self.beta) * self.f_i_minus_1
                
                           + d4(pheta, self.alpha, self.beta) * self.y_i         
                + self.h_i * d5(pheta, self.alpha, self.beta) * self.f_i 
                
                           + d6(pheta, self.alpha, self.beta) * self.y_i_plus_1  
                + self.h_i * d7(pheta, self.alpha, self.beta) * self.f_i_plus_1
        )

    def prime(self, x):
        pheta = (x - self.x_i) / self.h_i  # x = t_i + pheta*h_i so pheta = (x - t_i) / h_i
        return (  
              d0_prime(pheta, self.alpha, self.beta) * self.y_i_minus_2 / self.h_i 
            + d1_prime(pheta, self.alpha, self.beta) * self.f_i_minus_2

            + d2_prime(pheta, self.alpha, self.beta) * self.y_i_minus_1 / self.h_i 
            + d3_prime(pheta, self.alpha, self.beta) * self.f_i_minus_1

            + d4_prime(pheta, self.alpha, self.beta) * self.y_i         / self.h_i 
            + d5_prime(pheta, self.alpha, self.beta) * self.f_i 

            + d6_prime(pheta, self.alpha, self.beta) * self.y_i_plus_1  / self.h_i 
            + d7_prime(pheta, self.alpha, self.beta) * self.f_i_plus_1
        )

    def eval_horner(self, x):
        pheta = (x - self.x_i) / self.h_i  # x = t_i + pheta*h_i so pheta = (x - t_i) / h_i
        return (  
                             d0_horner(pheta, self.alpha, self.beta) * self.y_i_minus_2 
                + self.h_i * d1_horner(pheta, self.alpha, self.beta) * self.f_i_minus_2

                           + d2_horner(pheta, self.alpha, self.beta) * self.y_i_minus_1 
                + self.h_i * d3_horner(pheta, self.alpha, self.beta) * self.f_i_minus_1
                
                           + d4_horner(pheta, self.alpha, self.beta) * self.y_i         
                + self.h_i * d5_horner(pheta, self.alpha, self.beta) * self.f_i 
                
                           + d6_horner(pheta, self.alpha, self.beta) * self.y_i_plus_1  
                + self.h_i * d7_horner(pheta, self.alpha, self.beta) * self.f_i_plus_1
        )

    def prime_horner(self, x):
        pheta = (x - self.x_i) / self.h_i  # x = t_i + pheta*h_i so pheta = (x - t_i) / h_i
        return (  
              d0_prime_horner(pheta, self.alpha, self.beta) * self.y_i_minus_2 / self.h_i 
            + d1_prime_horner(pheta, self.alpha, self.beta) * self.f_i_minus_2

            + d2_prime_horner(pheta, self.alpha, self.beta) * self.y_i_minus_1 / self.h_i 
            + d3_prime_horner(pheta, self.alpha, self.beta) * self.f_i_minus_1

            + d4_prime_horner(pheta, self.alpha, self.beta) * self.y_i         / self.h_i 
            + d5_prime_horner(pheta, self.alpha, self.beta) * self.f_i 

            + d6_prime_horner(pheta, self.alpha, self.beta) * self.y_i_plus_1  / self.h_i 
            + d7_prime_horner(pheta, self.alpha, self.beta) * self.f_i_plus_1
        )
    
    def eval_bary(self, x) -> float:
        return self.eval_bary_interp(x)

    def prime_bary(self, x) -> float:
        return self.prime_bary_interp(x)