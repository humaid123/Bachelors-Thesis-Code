
from math import cos, pi
from scipy.interpolate import BarycentricInterpolator

def get_Chebyshev_nodes(a, b, n):
    res = []
    for k in range(1, n+1):
        res.append(
            (a+b)/2 + (b-a)/2 * cos( (2*k - 1) / (2*n) * pi)
        )
    return res

class ContinuousSolution:
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
            HB (
                x_i_minus_2, x_i_minus_1, x_i, x_i_plus_1,
                y_i_minus_2, f_i_minus_2,
                y_i_minus_1, f_i_minus_1,
                y_i, f_i,
                y_i_plus_1, f_i_plus_1,
                monitor 
            )
        )
    continuous_sol = ContinuousSolution()
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
        
        interp = HB (
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
    return ( (2*beta + 1)/(10*alpha**4*beta + 5*alpha**4 + 10*alpha**3*beta**2 + 30*alpha**3*beta + 12*alpha**3 + 15*alpha**2*beta**2 + 27*alpha**2*beta + 9*alpha**2 + 5*alpha*beta**2 + 7*alpha*beta + 2*alpha)*(x**5) + (-5*beta**2 + 5*beta + 4)/(20*alpha**4*beta + 10*alpha**4 + 20*alpha**3*beta**2 + 60*alpha**3*beta + 24*alpha**3 + 30*alpha**2*beta**2 + 54*alpha**2*beta + 18*alpha**2 + 10*alpha*beta**2 + 14*alpha*beta + 4*alpha)*(x**4) + (-5*beta**2 - 
beta + 1)/(10*alpha**4*beta + 5*alpha**4 + 10*alpha**3*beta**2 + 30*alpha**3*beta + 12*alpha**3 + 15*alpha**2*beta**2 + 27*alpha**2*beta + 9*alpha**2 + 
5*alpha*beta**2 + 7*alpha*beta + 2*alpha)*(x**3) + (-5*beta**2 - 3*beta)/(20*alpha**4*beta + 10*alpha**4 + 20*alpha**3*beta**2 + 60*alpha**3*beta + 24*alpha**3 + 30*alpha**2*beta**2 + 54*alpha**2*beta + 18*alpha**2 + 10*alpha*beta**2 + 14*alpha*beta + 4*alpha)*(x**2) + 0*x + 0 )
def d1(x, alpha, beta):
    return ( -12/(10*alpha*beta + 5*alpha + 5*beta + 2)*(x**5) + (-15*alpha + 15*beta - 30)/(10*alpha*beta + 5*alpha + 5*beta + 2)*(x**4) + (20*alpha*beta - 20*alpha + 40*beta - 20)/(10*alpha*beta + 5*alpha + 5*beta + 2)*(x**3) + (30*alpha*beta + 30*beta)/(10*alpha*beta + 5*alpha + 5*beta + 2)*(x**2) + 
0*x + 0 )
def d2(x, alpha, beta):
    return ( (-6*alpha*beta - 4*alpha - 2*beta - 1)/(10*alpha**2*beta**2 + 15*alpha**2*beta + 5*alpha**2 + 5*alpha*beta**2 + 7*alpha*beta + 2*alpha)*(x**5) + (-15*alpha**2*beta - 10*alpha**2 + 15*alpha*beta**2 - 20*alpha*beta - 20*alpha + 5*beta**2 - 5*beta - 4)/(20*alpha**2*beta**2 + 30*alpha**2*beta + 10*alpha**2 + 10*alpha*beta**2 + 14*alpha*beta + 4*alpha)*(x**4) + (10*alpha**2*beta**2 - 5*alpha**2 + 20*alpha*beta**2 + 5*alpha*beta - 6*alpha + 5*beta**2 + beta - 1)/(10*alpha**2*beta**2 + 15*alpha**2*beta + 5*alpha**2 + 5*alpha*beta**2 + 7*alpha*beta + 2*alpha)*(x**3) + (20*alpha**2*beta**2 + 15*alpha**2*beta + 25*alpha*beta**2 + 18*alpha*beta + 5*beta**2 + 3*beta)/(20*alpha**2*beta**2 + 30*alpha**2*beta + 10*alpha**2 + 10*alpha*beta**2 + 14*alpha*beta + 4*alpha)*(x**2) + 0*x + 0 )
def d3(x, alpha, beta):
    return ( 12/(10*alpha*beta + 5*alpha + 5*beta + 2)*(x**5) + (15*alpha - 15*beta + 30)/(10*alpha*beta + 5*alpha + 5*beta + 2)*(x**4) + (-20*alpha*beta + 20*alpha - 40*beta + 20)/(10*alpha*beta + 5*alpha + 5*beta + 2)*(x**3) + (-30*alpha*beta - 30*beta)/(10*alpha*beta + 5*alpha + 5*beta + 2)*(x**2) + 
0*x + 1 )
def d4(x, alpha, beta):
    return ( (-6*alpha*beta - 2*alpha - 4*beta - 1)/(10*alpha**2*beta**2 + 5*alpha**2*beta + 15*alpha*beta**2 + 7*alpha*beta + 5*beta**2 + 2*beta)*(x**5) + (-15*alpha**2*beta - 5*alpha**2 + 15*alpha*beta**2 - 40*alpha*beta - 15*alpha + 10*beta**2 - 20*beta - 6)/(20*alpha**2*beta**2 + 10*alpha**2*beta + 
30*alpha*beta**2 + 14*alpha*beta + 10*beta**2 + 4*beta)*(x**4) + (10*alpha**2*beta**2 - 10*alpha**2*beta - 5*alpha**2 + 30*alpha*beta**2 - 15*alpha*beta - 9*alpha + 15*beta**2 - 6*beta - 3)/(10*alpha**2*beta**2 + 5*alpha**2*beta + 15*alpha*beta**2 + 7*alpha*beta + 5*beta**2 + 2*beta)*(x**3) + (40*alpha**2*beta**2 + 5*alpha**2*beta - 5*alpha**2 + 75*alpha*beta**2 + 12*alpha*beta - 7*alpha + 30*beta**2 + 4*beta - 2)/(20*alpha**2*beta**2 + 10*alpha**2*beta + 30*alpha*beta**2 + 14*alpha*beta + 10*beta**2 + 4*beta)*(x**2) + 1*x + 0 )
def d5(x, alpha, beta):
    return ( (2*alpha + 1)/(10*alpha**2*beta**3 + 15*alpha**2*beta**2 + 5*alpha**2*beta + 10*alpha*beta**4 + 30*alpha*beta**3 + 27*alpha*beta**2 + 7*alpha*beta + 5*beta**4 + 12*beta**3 + 9*beta**2 + 2*beta)*(x**5) + (5*alpha**2 + 15*alpha + 6)/(20*alpha**2*beta**3 + 30*alpha**2*beta**2 + 10*alpha**2*beta + 20*alpha*beta**4 + 60*alpha*beta**3 + 54*alpha*beta**2 + 14*alpha*beta + 10*beta**4 + 24*beta**3 + 18*beta**2 + 4*beta)*(x**4) + (5*alpha**2 + 9*alpha + 3)/(10*alpha**2*beta**3 + 15*alpha**2*beta**2 + 5*alpha**2*beta + 10*alpha*beta**4 + 30*alpha*beta**3 + 27*alpha*beta**2 + 7*alpha*beta + 5*beta**4 + 12*beta**3 + 9*beta**2 + 2*beta)*(x**3) + (5*alpha**2 + 7*alpha + 2)/(20*alpha**2*beta**3 + 30*alpha**2*beta**2 + 10*alpha**2*beta + 20*alpha*beta**4 + 60*alpha*beta**3 + 54*alpha*beta**2 + 14*alpha*beta + 10*beta**4 + 24*beta**3 + 18*beta**2 + 4*beta)*(x**2) + 0*x + 0 )
def d0_prime(x, alpha, beta):
    return ( 5*(2*beta + 1)/(10*alpha**4*beta + 5*alpha**4 + 10*alpha**3*beta**2 + 30*alpha**3*beta + 12*alpha**3 + 15*alpha**2*beta**2 + 27*alpha**2*beta + 9*alpha**2 + 5*alpha*beta**2 + 7*alpha*beta + 2*alpha)*(x**4) + 4*(-5*beta**2 + 5*beta + 4)/(20*alpha**4*beta + 10*alpha**4 + 20*alpha**3*beta**2 + 60*alpha**3*beta + 24*alpha**3 + 30*alpha**2*beta**2 + 54*alpha**2*beta + 18*alpha**2 + 10*alpha*beta**2 + 14*alpha*beta + 4*alpha)*(x**3) + 3*(-5*beta**2 - beta + 1)/(10*alpha**4*beta + 5*alpha**4 + 10*alpha**3*beta**2 + 30*alpha**3*beta + 12*alpha**3 + 15*alpha**2*beta**2 + 27*alpha**2*beta + 9*alpha**2 + 5*alpha*beta**2 + 7*alpha*beta + 2*alpha)*(x**2) + 2*(-5*beta**2 - 3*beta)/(20*alpha**4*beta + 10*alpha**4 + 20*alpha**3*beta**2 + 60*alpha**3*beta + 24*alpha**3 + 30*alpha**2*beta**2 + 54*alpha**2*beta + 18*alpha**2 + 10*alpha*beta**2 + 14*alpha*beta + 4*alpha)*x + 0 )
def d1_prime(x, alpha, beta):
    return ( 5*-12/(10*alpha*beta + 5*alpha + 5*beta + 2)*(x**4) + 4*(-15*alpha + 15*beta - 30)/(10*alpha*beta + 5*alpha + 5*beta + 2)*(x**3) + 3*(20*alpha*beta - 20*alpha + 40*beta - 20)/(10*alpha*beta + 5*alpha + 5*beta + 2)*(x**2) + 2*(30*alpha*beta + 30*beta)/(10*alpha*beta + 5*alpha + 5*beta + 2)*x + 0 )
def d2_prime(x, alpha, beta):
    return ( 5*(-6*alpha*beta - 4*alpha - 2*beta - 1)/(10*alpha**2*beta**2 + 15*alpha**2*beta + 5*alpha**2 + 5*alpha*beta**2 + 7*alpha*beta + 2*alpha)*(x**4) + 4*(-15*alpha**2*beta - 10*alpha**2 + 15*alpha*beta**2 - 20*alpha*beta - 20*alpha + 5*beta**2 - 5*beta - 4)/(20*alpha**2*beta**2 + 30*alpha**2*beta + 10*alpha**2 + 10*alpha*beta**2 + 14*alpha*beta + 4*alpha)*(x**3) + 3*(10*alpha**2*beta**2 - 5*alpha**2 + 20*alpha*beta**2 + 5*alpha*beta - 6*alpha 
+ 5*beta**2 + beta - 1)/(10*alpha**2*beta**2 + 15*alpha**2*beta + 5*alpha**2 + 5*alpha*beta**2 + 7*alpha*beta + 2*alpha)*(x**2) + 2*(20*alpha**2*beta**2 + 15*alpha**2*beta + 25*alpha*beta**2 + 18*alpha*beta + 5*beta**2 + 3*beta)/(20*alpha**2*beta**2 + 30*alpha**2*beta + 10*alpha**2 + 10*alpha*beta**2 + 
14*alpha*beta + 4*alpha)*x + 0 )
def d3_prime(x, alpha, beta):
    return ( 5*12/(10*alpha*beta + 5*alpha + 5*beta + 2)*(x**4) + 4*(15*alpha - 15*beta + 30)/(10*alpha*beta + 5*alpha + 5*beta + 2)*(x**3) + 3*(-20*alpha*beta + 20*alpha - 40*beta + 20)/(10*alpha*beta + 5*alpha + 5*beta + 2)*(x**2) + 2*(-30*alpha*beta - 30*beta)/(10*alpha*beta + 5*alpha + 5*beta + 2)*x + 0 )
def d4_prime(x, alpha, beta):
    return ( 5*(-6*alpha*beta - 2*alpha - 4*beta - 1)/(10*alpha**2*beta**2 + 5*alpha**2*beta + 15*alpha*beta**2 + 7*alpha*beta + 5*beta**2 + 2*beta)*(x**4) + 4*(-15*alpha**2*beta - 5*alpha**2 + 15*alpha*beta**2 - 40*alpha*beta - 15*alpha + 10*beta**2 - 20*beta - 6)/(20*alpha**2*beta**2 + 10*alpha**2*beta + 30*alpha*beta**2 + 14*alpha*beta + 10*beta**2 + 4*beta)*(x**3) + 3*(10*alpha**2*beta**2 - 10*alpha**2*beta - 5*alpha**2 + 30*alpha*beta**2 - 15*alpha*beta - 9*alpha + 15*beta**2 - 6*beta - 3)/(10*alpha**2*beta**2 + 5*alpha**2*beta + 15*alpha*beta**2 + 7*alpha*beta + 5*beta**2 + 2*beta)*(x**2) + 2*(40*alpha**2*beta**2 + 5*alpha**2*beta - 5*alpha**2 + 75*alpha*beta**2 + 12*alpha*beta - 7*alpha + 30*beta**2 + 4*beta - 2)/(20*alpha**2*beta**2 + 10*alpha**2*beta + 30*alpha*beta**2 + 14*alpha*beta + 10*beta**2 + 4*beta)*x + 1 )
def d5_prime(x, alpha, beta):
    return ( 5*(2*alpha + 1)/(10*alpha**2*beta**3 + 15*alpha**2*beta**2 + 5*alpha**2*beta + 10*alpha*beta**4 + 30*alpha*beta**3 + 27*alpha*beta**2 + 7*alpha*beta + 5*beta**4 + 12*beta**3 + 9*beta**2 + 2*beta)*(x**4) + 4*(5*alpha**2 + 15*alpha + 6)/(20*alpha**2*beta**3 + 30*alpha**2*beta**2 + 10*alpha**2*beta + 20*alpha*beta**4 + 60*alpha*beta**3 + 54*alpha*beta**2 + 14*alpha*beta + 10*beta**4 + 24*beta**3 + 18*beta**2 + 4*beta)*(x**3) + 3*(5*alpha**2 + 9*alpha + 3)/(10*alpha**2*beta**3 + 15*alpha**2*beta**2 + 5*alpha**2*beta + 10*alpha*beta**4 + 30*alpha*beta**3 + 27*alpha*beta**2 + 7*alpha*beta + 5*beta**4 + 12*beta**3 + 9*beta**2 + 2*beta)*(x**2) + 2*(5*alpha**2 + 7*alpha + 2)/(20*alpha**2*beta**3 + 30*alpha**2*beta**2 + 10*alpha**2*beta + 20*alpha*beta**4 + 60*alpha*beta**3 + 54*alpha*beta**2 + 14*alpha*beta + 10*beta**4 + 24*beta**3 + 18*beta**2 + 4*beta)*x + 0 )
#########################################################################################################
def d0_horner(x, alpha, beta) -> float:
    return ( x**2*(beta*(-5*beta - 3)/(20*alpha**4*beta + 10*alpha**4 + 20*alpha**3*beta**2 + 60*alpha**3*beta + 24*alpha**3 + 30*alpha**2*beta**2 + 54*alpha**2*beta + 18*alpha**2 + 10*alpha*beta**2 + 14*alpha*beta + 4*alpha) + x*(x*(x*(2*beta + 1)/(10*alpha**4*beta + 5*alpha**4 + 10*alpha**3*beta**2 + 
30*alpha**3*beta + 12*alpha**3 + 15*alpha**2*beta**2 + 27*alpha**2*beta + 9*alpha**2 + 5*alpha*beta**2 + 7*alpha*beta + 2*alpha) + (beta*(5 - 5*beta) + 
4)/(20*alpha**4*beta + 10*alpha**4 + 20*alpha**3*beta**2 + 60*alpha**3*beta + 24*alpha**3 + 30*alpha**2*beta**2 + 54*alpha**2*beta + 18*alpha**2 + 10*alpha*beta**2 + 14*alpha*beta + 4*alpha)) + (beta*(-5*beta - 1) + 1)/(10*alpha**4*beta + 5*alpha**4 + 10*alpha**3*beta**2 + 30*alpha**3*beta + 12*alpha**3 + 15*alpha**2*beta**2 + 27*alpha**2*beta + 9*alpha**2 + 5*alpha*beta**2 + 7*alpha*beta + 2*alpha))) )

def d1_horner(x, alpha, beta) -> float:
    return ( x**2*(x*(x*(-12*x/(10*alpha*beta + 5*alpha + 5*beta + 2) + (-15*alpha + 15*beta - 30)/(10*alpha*beta + 5*alpha + 5*beta + 2)) + (alpha*(20*beta - 20) + 40*beta - 20)/(10*alpha*beta + 5*alpha + 5*beta + 2)) + (30*alpha*beta + 30*beta)/(10*alpha*beta + 5*alpha + 5*beta + 2)) )

def d2_horner(x, alpha, beta) -> float:
    return ( x**2*(x*(x*(x*(alpha*(-6*beta - 4) - 2*beta - 1)/(10*alpha**2*beta**2 + 15*alpha**2*beta + 5*alpha**2 + 5*alpha*beta**2 + 7*alpha*beta + 2*alpha) + (alpha*(alpha*(-15*beta - 10) + beta*(15*beta - 20) - 20) + beta*(5*beta - 5) - 4)/(20*alpha**2*beta**2 + 30*alpha**2*beta + 10*alpha**2 + 10*alpha*beta**2 + 14*alpha*beta + 4*alpha)) + (alpha*(alpha*(10*beta**2 - 5) + beta*(20*beta + 5) - 6) + beta*(5*beta + 1) - 1)/(10*alpha**2*beta**2 + 15*alpha**2*beta + 5*alpha**2 + 5*alpha*beta**2 + 7*alpha*beta + 2*alpha)) + (alpha*(alpha*beta*(20*beta + 15) + beta*(25*beta + 18)) + beta*(5*beta + 3))/(20*alpha**2*beta**2 + 30*alpha**2*beta + 10*alpha**2 + 10*alpha*beta**2 + 14*alpha*beta + 4*alpha)) )

def d3_horner(x, alpha, beta) -> float:
    return ( x**2*(x*(x*(12*x/(10*alpha*beta + 5*alpha + 5*beta + 2) + (15*alpha - 15*beta + 30)/(10*alpha*beta + 5*alpha + 5*beta + 2)) + (alpha*(20 - 
20*beta) - 40*beta + 20)/(10*alpha*beta + 5*alpha + 5*beta + 2)) + (-30*alpha*beta - 30*beta)/(10*alpha*beta + 5*alpha + 5*beta + 2)) + 1 )

def d4_horner(x, alpha, beta) -> float:
    return ( x*(x*(x*(x*(x*(alpha*(-6*beta - 2) - 4*beta - 1)/(10*alpha**2*beta**2 + 5*alpha**2*beta + 15*alpha*beta**2 + 7*alpha*beta + 5*beta**2 + 2*beta) + (alpha*(alpha*(-15*beta - 5) + beta*(15*beta - 40) - 15) + beta*(10*beta - 20) - 6)/(20*alpha**2*beta**2 + 10*alpha**2*beta + 30*alpha*beta**2 + 
14*alpha*beta + 10*beta**2 + 4*beta)) + (alpha*(alpha*(beta*(10*beta - 10) - 5) + beta*(30*beta - 15) - 9) + beta*(15*beta - 6) - 3)/(10*alpha**2*beta**2 + 5*alpha**2*beta + 15*alpha*beta**2 + 7*alpha*beta + 5*beta**2 + 2*beta)) + (alpha*(alpha*(beta*(40*beta + 5) - 5) + beta*(75*beta + 12) - 7) + beta*(30*beta + 4) - 2)/(20*alpha**2*beta**2 + 10*alpha**2*beta + 30*alpha*beta**2 + 14*alpha*beta + 10*beta**2 + 4*beta)) + 1) )

def d5_horner(x, alpha, beta) -> float:
    return ( x**2*(x*(x*(x*(2*alpha + 1)/(10*alpha**2*beta**3 + 15*alpha**2*beta**2 + 5*alpha**2*beta + 10*alpha*beta**4 + 30*alpha*beta**3 + 27*alpha*beta**2 + 7*alpha*beta + 5*beta**4 + 12*beta**3 + 9*beta**2 + 2*beta) + (alpha*(5*alpha + 15) + 6)/(20*alpha**2*beta**3 + 30*alpha**2*beta**2 + 10*alpha**2*beta + 20*alpha*beta**4 + 60*alpha*beta**3 + 54*alpha*beta**2 + 14*alpha*beta + 10*beta**4 + 24*beta**3 + 18*beta**2 + 4*beta)) + (alpha*(5*alpha + 9) + 3)/(10*alpha**2*beta**3 + 15*alpha**2*beta**2 + 5*alpha**2*beta + 10*alpha*beta**4 + 30*alpha*beta**3 + 27*alpha*beta**2 + 7*alpha*beta + 5*beta**4 
+ 12*beta**3 + 9*beta**2 + 2*beta)) + (alpha*(5*alpha + 7) + 2)/(20*alpha**2*beta**3 + 30*alpha**2*beta**2 + 10*alpha**2*beta + 20*alpha*beta**4 + 60*alpha*beta**3 + 54*alpha*beta**2 + 14*alpha*beta + 10*beta**4 + 24*beta**3 + 18*beta**2 + 4*beta)) )

def d0_prime_horner(x, alpha, beta) -> float:
    return ( x*(beta*(-5*beta - 3)/(10*alpha**4*beta + 5*alpha**4 + 10*alpha**3*beta**2 + 30*alpha**3*beta + 12*alpha**3 + 15*alpha**2*beta**2 + 27*alpha**2*beta + 9*alpha**2 + 5*alpha*beta**2 + 7*alpha*beta + 2*alpha) + x*(x*(x*(10*beta + 5)/(10*alpha**4*beta + 5*alpha**4 + 10*alpha**3*beta**2 + 30*alpha**3*beta + 12*alpha**3 + 15*alpha**2*beta**2 + 27*alpha**2*beta + 9*alpha**2 + 5*alpha*beta**2 + 7*alpha*beta + 2*alpha) + (beta*(10 - 10*beta) + 8)/(10*alpha**4*beta + 5*alpha**4 + 10*alpha**3*beta**2 + 30*alpha**3*beta + 12*alpha**3 + 15*alpha**2*beta**2 + 27*alpha**2*beta + 9*alpha**2 + 5*alpha*beta**2 + 7*alpha*beta + 2*alpha)) + (beta*(-15*beta - 3) + 3)/(10*alpha**4*beta + 5*alpha**4 + 10*alpha**3*beta**2 + 30*alpha**3*beta + 12*alpha**3 + 15*alpha**2*beta**2 + 27*alpha**2*beta + 9*alpha**2 + 5*alpha*beta**2 + 7*alpha*beta + 2*alpha))) )

def d1_prime_horner(x, alpha, beta) -> float:
    return ( x*(x*(x*(-60*x/(10*alpha*beta + 5*alpha + 5*beta + 2) + (-60*alpha + 60*beta - 120)/(10*alpha*beta + 5*alpha + 5*beta + 2)) + (alpha*(60*beta - 60) + 120*beta - 60)/(10*alpha*beta + 5*alpha + 5*beta + 2)) + (60*alpha*beta + 60*beta)/(10*alpha*beta + 5*alpha + 5*beta + 2)) )

def d2_prime_horner(x, alpha, beta) -> float:
    return ( x*(x*(x*(x*(alpha*(-30*beta - 20) - 10*beta - 5)/(10*alpha**2*beta**2 + 15*alpha**2*beta + 5*alpha**2 + 5*alpha*beta**2 + 7*alpha*beta + 2*alpha) + (alpha*(alpha*(-30*beta - 20) + beta*(30*beta - 40) - 40) + beta*(10*beta - 10) - 8)/(10*alpha**2*beta**2 + 15*alpha**2*beta + 5*alpha**2 + 5*alpha*beta**2 + 7*alpha*beta + 2*alpha)) + (alpha*(alpha*(30*beta**2 - 15) + beta*(60*beta + 15) - 18) + beta*(15*beta + 3) - 3)/(10*alpha**2*beta**2 + 15*alpha**2*beta + 5*alpha**2 + 5*alpha*beta**2 + 7*alpha*beta + 2*alpha)) + (alpha*(alpha*beta*(20*beta + 15) + beta*(25*beta + 18)) + beta*(5*beta + 3))/(10*alpha**2*beta**2 + 15*alpha**2*beta + 5*alpha**2 + 5*alpha*beta**2 + 7*alpha*beta + 2*alpha)) )

def d3_prime_horner(x, alpha, beta) -> float:
    return ( x*(x*(x*(60*x/(10*alpha*beta + 5*alpha + 5*beta + 2) + (60*alpha - 60*beta + 120)/(10*alpha*beta + 5*alpha + 5*beta + 2)) + (alpha*(60 - 60*beta) - 120*beta + 60)/(10*alpha*beta + 5*alpha + 5*beta + 2)) + (-60*alpha*beta - 60*beta)/(10*alpha*beta + 5*alpha + 5*beta + 2)) )

def d4_prime_horner(x, alpha, beta) -> float:
    return ( x*(x*(x*(x*(alpha*(-30*beta - 10) - 20*beta - 5)/(10*alpha**2*beta**2 + 5*alpha**2*beta + 15*alpha*beta**2 + 7*alpha*beta + 5*beta**2 + 2*beta) + (alpha*(alpha*(-30*beta - 10) + beta*(30*beta - 80) - 30) + beta*(20*beta - 40) - 12)/(10*alpha**2*beta**2 + 5*alpha**2*beta + 15*alpha*beta**2 + 7*alpha*beta + 5*beta**2 + 2*beta)) + (alpha*(alpha*(beta*(30*beta - 30) - 15) + beta*(90*beta - 45) - 27) + beta*(45*beta - 18) - 9)/(10*alpha**2*beta**2 + 5*alpha**2*beta + 15*alpha*beta**2 + 7*alpha*beta + 5*beta**2 + 2*beta)) + (alpha*(alpha*(beta*(40*beta + 5) - 5) + beta*(75*beta + 12) - 7) + beta*(30*beta + 4) - 2)/(10*alpha**2*beta**2 + 5*alpha**2*beta + 15*alpha*beta**2 + 7*alpha*beta + 5*beta**2 + 2*beta)) + 1 )

def d5_prime_horner(x, alpha, beta) -> float:
    return ( x*(x*(x*(x*(10*alpha + 5)/(10*alpha**2*beta**3 + 15*alpha**2*beta**2 + 5*alpha**2*beta + 10*alpha*beta**4 + 30*alpha*beta**3 + 27*alpha*beta**2 + 7*alpha*beta + 5*beta**4 + 12*beta**3 + 9*beta**2 + 2*beta) + (alpha*(10*alpha + 30) + 12)/(10*alpha**2*beta**3 + 15*alpha**2*beta**2 + 5*alpha**2*beta + 10*alpha*beta**4 + 30*alpha*beta**3 + 27*alpha*beta**2 + 7*alpha*beta + 5*beta**4 + 12*beta**3 + 9*beta**2 + 2*beta)) + (alpha*(15*alpha + 27) 
+ 9)/(10*alpha**2*beta**3 + 15*alpha**2*beta**2 + 5*alpha**2*beta + 10*alpha*beta**4 + 30*alpha*beta**3 + 27*alpha*beta**2 + 7*alpha*beta + 5*beta**4 + 
12*beta**3 + 9*beta**2 + 2*beta)) + (alpha*(5*alpha + 7) + 2)/(10*alpha**2*beta**3 + 15*alpha**2*beta**2 + 5*alpha**2*beta + 10*alpha*beta**4 + 30*alpha*beta**3 + 27*alpha*beta**2 + 7*alpha*beta + 5*beta**4 + 12*beta**3 + 9*beta**2 + 2*beta)) )


#########################################################################################################

# I use a class to represent the Hermite Birkhoff interpolant
# we will have an instance of this class on each step
class HB:
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
                + self.h_i * d0(pheta, self.alpha, self.beta) * self.f_i_minus_2

                           + d1(pheta, self.alpha, self.beta) * self.y_i_minus_1 
                + self.h_i * d2(pheta, self.alpha, self.beta) * self.f_i_minus_1
                
                           + d3(pheta, self.alpha, self.beta) * self.y_i         
                + self.h_i * d4(pheta, self.alpha, self.beta) * self.f_i 
                 
                + self.h_i * d5(pheta, self.alpha, self.beta) * self.f_i_plus_1
        )

    def prime(self, x):
        pheta = (x - self.x_i) / self.h_i  # x = t_i + pheta*h_i so pheta = (x - t_i) / h_i
        return (  
              d0_prime(pheta, self.alpha, self.beta) * self.f_i_minus_2

            + d1_prime(pheta, self.alpha, self.beta) * self.y_i_minus_1 / self.h_i 
            + d2_prime(pheta, self.alpha, self.beta) * self.f_i_minus_1

            + d3_prime(pheta, self.alpha, self.beta) * self.y_i         / self.h_i 
            + d4_prime(pheta, self.alpha, self.beta) * self.f_i 

            + d5_prime(pheta, self.alpha, self.beta) * self.f_i_plus_1
        )

    def eval_horner(self, x):
        pheta = (x - self.x_i) / self.h_i  # x = t_i + pheta*h_i so pheta = (x - t_i) / h_i
        return (  
                + self.h_i * d0_horner(pheta, self.alpha, self.beta) * self.f_i_minus_2

                           + d1_horner(pheta, self.alpha, self.beta) * self.y_i_minus_1 
                + self.h_i * d2_horner(pheta, self.alpha, self.beta) * self.f_i_minus_1
                
                           + d3_horner(pheta, self.alpha, self.beta) * self.y_i         
                + self.h_i * d4_horner(pheta, self.alpha, self.beta) * self.f_i 
                
                + self.h_i * d5_horner(pheta, self.alpha, self.beta) * self.f_i_plus_1
        )

    def prime_horner(self, x):
        pheta = (x - self.x_i) / self.h_i  # x = t_i + pheta*h_i so pheta = (x - t_i) / h_i
        return (  
              d0_prime_horner(pheta, self.alpha, self.beta) * self.f_i_minus_2

            + d1_prime_horner(pheta, self.alpha, self.beta) * self.y_i_minus_1 / self.h_i 
            + d2_prime_horner(pheta, self.alpha, self.beta) * self.f_i_minus_1

            + d3_prime_horner(pheta, self.alpha, self.beta) * self.y_i         / self.h_i 
            + d4_prime_horner(pheta, self.alpha, self.beta) * self.f_i 

            + d5_prime_horner(pheta, self.alpha, self.beta) * self.f_i_plus_1
        )
    
    def eval_bary(self, x) -> float:
        return self.eval_bary_interp(x)

    def prime_bary(self, x) -> float:
        return self.prime_bary_interp(x)