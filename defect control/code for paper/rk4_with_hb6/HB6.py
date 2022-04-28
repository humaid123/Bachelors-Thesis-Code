# this is the final file that is copy and pasted everytime we need HB6

# these formulae were derived using sympy: see get_hb_coefficients.py and test_hb_coefficients.py for the derivations and the tests
# we create a class for each of these so that the class abstracts a Hermite Birkhoff interpolant
# we then need to instantiate one object of this class for each interval and we get to do Hermite Birkhoff interpolation

from math import cos, pi
from scipy.interpolate import BarycentricInterpolator

def get_Chebyshev_nodes(a, b, n):
    res = []
    for k in range(1, n+1):
        res.append(
            (a + b)/2 + (b-a)/2 * cos( (2*k - 1) / (2*n) * pi)
        )
    return res

def d0(x, alpha):
    return (4*alpha + 2)/(alpha**6 + 3*alpha**5 + 3*alpha**4 + alpha**3)*(x**5) + (5*alpha**2 - 5*alpha - 4)/(alpha**6 + 3*alpha**5 + 3*alpha**4 + alpha**3)*(x**4) + (-10*alpha**2 - 2*alpha + 2)/(alpha**6 + 3*alpha**5 + 3*alpha**4 + alpha**3)*(x**3) + (5*alpha + 3)/(alpha**5 + 3*alpha**4 + 3*alpha**3 + alpha**2)*(x**2) + 0*x + 0
def d1(x, alpha):
    return 1/(alpha**4 + 2*alpha**3 + alpha**2)*(x**5) + (alpha - 2)/(alpha**4 + 2*alpha**3 + alpha**2)*(x**4) + (1 - 2*alpha)/(alpha**4 + 2*alpha**3 + alpha**2)*(x**3) + 1/(alpha**3 + 2*alpha**2 + alpha)*(x**2) + 0*x + 0
def d2(x, alpha):
    return (2*alpha - 2)/alpha**3*(x**5) + (4*alpha**2 - 7*alpha + 4)/alpha**3*(x**4) + (2*alpha**3 - 8*alpha**2 + 8*alpha - 2)/alpha**3*(x**3) + (-3*alpha**2 + 4*alpha - 3)/alpha**2*(x**2) + 0*x + 1
def d3(x, alpha):
    return alpha**(-2)*(x**5) + (2*alpha - 2)/alpha**2*(x**4) + (alpha**2 - 4*alpha + 1)/alpha**2*(x**3) + (2 - 2*alpha)/alpha*(x**2) + 1*x + 0
def d4(x, alpha):
    return (-2*alpha - 4)/(alpha**3 + 3*alpha**2 + 3*alpha + 1)*(x**5) + (-4*alpha**2 - 5*alpha + 5)/(alpha**3 + 3*alpha**2 + 3*alpha + 1)*(x**4) + (-2*alpha**3 + 2*alpha**2 + 10*alpha)/(alpha**3 + 3*alpha**2 + 3*alpha + 1)*(x**3) + (3*alpha**3 + 5*alpha**2)/(alpha**3 + 3*alpha**2 + 3*alpha + 1)*(x**2) + 0*x + 0
def d5(x, alpha):
    return 1/(alpha**2 + 2*alpha + 1)*(x**5) + (2*alpha - 1)/(alpha**2 + 2*alpha + 1)*(x**4) + (alpha**2 - 2*alpha)/(alpha**2 + 2*alpha + 1)*(x**3) + -alpha**2/(alpha**2 + 2*alpha + 1)*(x**2) + 0*x + 0

def d0_prime(x, alpha):
    return 5*(4*alpha + 2)/(alpha**6 + 3*alpha**5 + 3*alpha**4 + alpha**3)*(x**4) + 4*(5*alpha**2 - 5*alpha - 4)/(alpha**6 + 3*alpha**5 + 3*alpha**4 + alpha**3)*(x**3) + 3*(-10*alpha**2 - 2*alpha + 2)/(alpha**6 + 3*alpha**5 + 3*alpha**4 + alpha**3)*(x**2) + 2*(5*alpha + 3)/(alpha**5 + 3*alpha**4 + 3*alpha**3 + alpha**2)*x + 0
def d1_prime(x, alpha):
    return 5*1/(alpha**4 + 2*alpha**3 + alpha**2)*(x**4) + 4*(alpha - 2)/(alpha**4 + 2*alpha**3 + alpha**2)*(x**3) + 3*(1 - 2*alpha)/(alpha**4 + 2*alpha**3 + alpha**2)*(x**2) + 2*1/(alpha**3 + 2*alpha**2 + alpha)*x + 0
def d2_prime(x, alpha):
    return 5*(2*alpha - 2)/alpha**3*(x**4) + 4*(4*alpha**2 - 7*alpha + 4)/alpha**3*(x**3) + 3*(2*alpha**3 - 8*alpha**2 + 8*alpha - 2)/alpha**3*(x**2) + 2*(-3*alpha**2 + 4*alpha - 3)/alpha**2*x + 0
def d3_prime(x, alpha):
    return 5*alpha**(-2)*(x**4) + 4*(2*alpha - 2)/alpha**2*(x**3) + 3*(alpha**2 - 4*alpha + 1)/alpha**2*(x**2) + 2*(2 - 2*alpha)/alpha*x + 1
def d4_prime(x, alpha):
    return 5*(-2*alpha - 4)/(alpha**3 + 3*alpha**2 + 3*alpha + 1)*(x**4) + 4*(-4*alpha**2 - 5*alpha + 5)/(alpha**3 + 3*alpha**2 + 3*alpha + 1)*(x**3) + 3*(-2*alpha**3 + 2*alpha**2 + 10*alpha)/(alpha**3 + 3*alpha**2 + 3*alpha + 1)*(x**2) + 2*(3*alpha**3 + 5*alpha**2)/(alpha**3 + 3*alpha**2 + 3*alpha + 1)*x + 0
def d5_prime(x, alpha):
    return 5*1/(alpha**2 + 2*alpha + 1)*(x**4) + 4*(2*alpha - 1)/(alpha**2 + 2*alpha + 1)*(x**3) + 3*(alpha**2 - 2*alpha)/(alpha**2 + 2*alpha + 1)*(x**2) + 2*-alpha**2/(alpha**2 + 2*alpha + 1)*x + 0


# ==========================================================================================
def d0_horner(x, alpha) -> float:
    return ( x**2*(x*(x*(x*(4*alpha + 2)/(alpha**6 + 3*alpha**5 + 3*alpha**4 + alpha**3) + (alpha*(5*alpha - 5) - 4)/(alpha**6 + 3*alpha**5 + 3*alpha**4 + alpha**3)) + (alpha*(-10*alpha - 2) + 2)/(alpha**6 + 3*alpha**5 + 3*alpha**4 + alpha**3)) + (5*alpha + 3)/(alpha**5 + 3*alpha**4 + 3*alpha**3 + alpha**2)) )

def d1_horner(x, alpha) -> float:
    return ( x**2*(x*(x*(x/(alpha**4 + 2*alpha**3 + alpha**2) + (alpha - 2)/(alpha**4 + 2*alpha**3 + alpha**2)) + (1 - 2*alpha)/(alpha**4 + 2*alpha**3 + alpha**2)) + 1/(alpha**3 + 2*alpha**2 + alpha)) )

def d2_horner(x, alpha) -> float:
    return ( x**2*(x*(x*((4 + (-7 + 4/alpha)/alpha)/alpha + x*(2 - 2/alpha)/alpha**2) + 2 + (-8 + (8 - 2/alpha)/alpha)/alpha) - 3 + (4 - 3/alpha)/alpha) + 1 )

def d3_horner(x, alpha) -> float:
    return ( x*(x*(x*(x*((2 - 2/alpha)/alpha + x/alpha**2) + 1 + (-4 + 1/alpha)/alpha) - 2 + 2/alpha) + 1) )

def d4_horner(x, alpha) -> float:
    return ( x**2*(alpha**2*(3*alpha + 5)/(alpha**3 + 3*alpha**2 + 3*alpha + 1) + x*(alpha*(alpha*(2 - 2*alpha) + 10)/(alpha**3 + 3*alpha**2 + 3*alpha + 1) + x*(x*(-2*alpha - 4)/(alpha**3 + 3*alpha**2 + 3*alpha + 1) + (alpha*(-4*alpha - 5) + 5)/(alpha**3 + 3*alpha**2 + 3*alpha + 1)))) )

def d5_horner(x, alpha) -> float:
    return ( x**2*(-alpha**2/(alpha**2 + 2*alpha + 1) + x*(alpha*(alpha - 2)/(alpha**2 + 2*alpha + 1) + x*(x/(alpha**2 + 2*alpha + 1) + (2*alpha - 1)/(alpha**2 + 2*alpha + 1)))) )

def d0_prime_horner(x, alpha) -> float:
    return ( x*(x*(x*(x*(20*alpha + 10)/(alpha**6 + 3*alpha**5 + 3*alpha**4 + alpha**3) + (alpha*(20*alpha - 20) - 16)/(alpha**6 + 3*alpha**5 + 3*alpha**4 + alpha**3)) + (alpha*(-30*alpha - 6) + 6)/(alpha**6 + 3*alpha**5 + 3*alpha**4 + alpha**3)) + (10*alpha + 6)/(alpha**5 + 3*alpha**4 + 3*alpha**3 + alpha**2)) )

def d1_prime_horner(x, alpha) -> float:
    return ( x*(x*(x*(5*x/(alpha**4 + 2*alpha**3 + alpha**2) + (4*alpha - 8)/(alpha**4 + 2*alpha**3 + alpha**2)) + (3 - 6*alpha)/(alpha**4 + 2*alpha**3 + alpha**2)) + 2/(alpha**3 + 2*alpha**2 + alpha)) )

def d2_prime_horner(x, alpha) -> float:
    return ( x*(x*(x*((16 + (-28 + 16/alpha)/alpha)/alpha + x*(10 - 10/alpha)/alpha**2) + 6 + (-24 + (24 - 6/alpha)/alpha)/alpha) - 6 + (8 - 6/alpha)/alpha) )

def d3_prime_horner(x, alpha) -> float:
    return ( x*(x*(x*((8 - 8/alpha)/alpha + 5*x/alpha**2) + 3 + (-12 + 3/alpha)/alpha) - 4 + 4/alpha) + 1 )

def d4_prime_horner(x, alpha) -> float:
    return ( x*(alpha**2*(6*alpha + 10)/(alpha**3 + 3*alpha**2 + 3*alpha + 1) + x*(alpha*(alpha*(6 - 6*alpha) + 30)/(alpha**3 + 3*alpha**2 + 3*alpha + 1) + x*(x*(-10*alpha - 20)/(alpha**3 + 3*alpha**2 + 3*alpha + 1) + (alpha*(-16*alpha - 20) + 20)/(alpha**3 + 3*alpha**2 + 3*alpha + 1)))) )

def d5_prime_horner(x, alpha) -> float:
    return ( x*(-2*alpha**2/(alpha**2 + 2*alpha + 1) + x*(alpha*(3*alpha - 6)/(alpha**2 + 2*alpha + 1) + x*(5*x/(alpha**2 + 2*alpha + 1) + (8*alpha - 4)/(alpha**2 + 2*alpha + 1)))) )

# ==========================================================================================

# I use a class to represent the Hermite Birkhoff interpolant
# we will have an instance of this class on each step
class HB:
    def __init__(   
        self, 
        x_i_minus_1, x_i, x_i_plus_1,
        y_i_minus_1, f_i_minus_1,
        y_i, f_i,
        y_i_plus_1, f_i_plus_1,
        monitor 
    ):
        h_i = x_i_plus_1 - x_i
        h_i_minus_1 = x_i - x_i_minus_1
        
        self.alpha = h_i_minus_1 / h_i

        monitor.different_values_alpha.add(self.alpha)

        self.h_i = h_i
        self.x_i = x_i

        # we also store x_i_minus_1 and x_i_plus_1 so that we can build the final interpolant
        self.x_i_plus_1 = x_i_plus_1
        self.x_i_minus_1 = x_i_minus_1

        self.y_i_minus_1 = y_i_minus_1 
        self.f_i_minus_1 = f_i_minus_1
        
        self.y_i = y_i
        self.f_i = f_i

        self.y_i_plus_1 = y_i_plus_1
        self.f_i_plus_1 = f_i_plus_1

        # getting barycentric interpolator =================
        xs = get_Chebyshev_nodes(x_i_minus_1, x_i_plus_1, 6)
        ys = [self.eval(x) for x in xs]
        self.eval_bary_interp = BarycentricInterpolator(xs, ys)

        y_primes = [self.prime(x) for x in xs]
        self.prime_bary_interp = BarycentricInterpolator(xs, y_primes) 


    def eval(self, x):
        pheta = (x - self.x_i) / self.h_i  # x = t_i + pheta*h_i so pheta = (x - t_i) / h_i
        return (  
                             d0(pheta, self.alpha) * self.y_i_minus_1 
                + self.h_i * d1(pheta, self.alpha) * self.f_i_minus_1
                
                           + d2(pheta, self.alpha) * self.y_i         
                + self.h_i * d3(pheta, self.alpha) * self.f_i 
                
                           + d4(pheta, self.alpha) * self.y_i_plus_1  
                + self.h_i * d5(pheta, self.alpha) * self.f_i_plus_1
        )

    def prime(self, x):
        pheta = (x - self.x_i) / self.h_i  # x = t_i + pheta*h_i so pheta = (x - t_i) / h_i
        return (  
              d0_prime(pheta, self.alpha) * self.y_i_minus_1 / self.h_i 
            + d1_prime(pheta, self.alpha) * self.f_i_minus_1

            + d2_prime(pheta, self.alpha) * self.y_i         / self.h_i 
            + d3_prime(pheta, self.alpha) * self.f_i 

            + d4_prime(pheta, self.alpha) * self.y_i_plus_1  / self.h_i 
            + d5_prime(pheta, self.alpha) * self.f_i_plus_1
        )

    def eval_horner(self, x):
        pheta = (x - self.x_i) / self.h_i  # x = t_i + pheta*h_i so pheta = (x - t_i) / h_i
        return (  
                             d0_horner(pheta, self.alpha) * self.y_i_minus_1 
                + self.h_i * d1_horner(pheta, self.alpha) * self.f_i_minus_1
                
                           + d2_horner(pheta, self.alpha) * self.y_i         
                + self.h_i * d3_horner(pheta, self.alpha) * self.f_i 
                
                           + d4_horner(pheta, self.alpha) * self.y_i_plus_1  
                + self.h_i * d5_horner(pheta, self.alpha) * self.f_i_plus_1
        )

    def prime_horner(self, x):
        pheta = (x - self.x_i) / self.h_i  # x = t_i + pheta*h_i so pheta = (x - t_i) / h_i
        return (  
              d0_prime_horner(pheta, self.alpha) * self.y_i_minus_1 / self.h_i 
            + d1_prime_horner(pheta, self.alpha) * self.f_i_minus_1

            + d2_prime_horner(pheta, self.alpha) * self.y_i         / self.h_i 
            + d3_prime_horner(pheta, self.alpha) * self.f_i 

            + d4_prime_horner(pheta, self.alpha) * self.y_i_plus_1  / self.h_i 
            + d5_prime_horner(pheta, self.alpha) * self.f_i_plus_1
        )
        
    def eval_bary(self, x) -> float:
        return self.eval_bary_interp(x)

    def prime_bary(self, x) -> float:
        return self.prime_bary_interp(x)

def create_defect_samplings(res, fn_s, monitor):
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
                y_i_plus_1, f_i_plus_1,
                monitor
        )
        result.append( (x_i_minus_1, x_i, x_i_plus_1, interp) )
    return result

def create_continuous_sol_from_results(res, fn_s, monitor):
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
                y_i_plus_1, f_i_plus_1,
                monitor
            )
        )
    continuous_sol = ContinuousSolution()
    continuous_sol.extend(interps)
    return continuous_sol

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
        print(f"ERROR in eval: {x} is outside of the solution range: {first_hb.x_i_minus_1} <= x <= {self.interps[-1].x_i_plus_1}")
        return -1

    def prime(self, x) -> float:
        for hb in self.interps:
            if (hb.x_i <= x <= hb.x_i_plus_1):
                return hb.prime(x)

        first_hb = self.interps[0]
        if (first_hb.x_i_minus_1 <= x <= first_hb.x_i):
            return first_hb.prime(x)
        
        print(f"ERROR in prime: {x} is outside of the solution range: {first_hb.x_i_minus_1} <= x <= {self.interps[-1].x_i_plus_1}")
        return -1
    
    def append(self, interp) -> None:
        self.interps.append(interp)
    
    def extend(self, newInterps) -> None:
        self.interps.extend(newInterps)