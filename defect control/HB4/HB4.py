
# classical Hermite Cubic.
# it has interpolation error of order 4
# taken from Wikipedia => https://en.wikipedia.org/wiki/Cubic_Hermite_spline

from typing import Tuple

def h00(t) -> float:
    return 2*(t**3) - 3*(t**2) + 1

def h10(t) -> float:
    return t**3 - 2*(t**2) + t

def h01(t) -> float:
    return -2*(t**3) + 3*(t**2)

def h11(t) -> float:
    return t**3 - t**2

def h00_prime(t) -> float:
    return 6*(t**2) - 6*t

def h10_prime(t) -> float:
    return 3*(t**2) - 4*t + 1

def h01_prime(t) -> float:
    return -6*(t**2) + 6*t

def h11_prime(t) -> float:
    return 3*(t**2) - 2*t

def h00_horner(x) -> float:
    return ( x**2*(2*x - 3) + 1 )

def h10_horner(x) -> float:
    return ( x*(x*(x - 2) + 1) )

def h01_horner(x) -> float:
    return ( x**2*(3 - 2*x) )

def h11_horner(x) -> float:
    return ( x**2*(x - 1) )

def h00_prime_horner(x) -> float:
    return ( x*(6*x - 6) )

def h10_prime_horner(x) -> float:
    return ( x*(3*x - 4) + 1 )

def h01_prime_horner(x) -> float:
    return ( x*(6 - 6*x) )

def h11_prime_horner(x) -> float:
    return ( x*(3*x - 2) )


class HB:
    def __init__(
        self, 
        x_i, x_i_plus_1, 
        y_i, f_i, 
        y_i_plus_1, f_i_plus_1
    ):
        self.h_i = x_i_plus_1 - x_i
        self.x_i = x_i
        self.x_i_plus_1 = x_i_plus_1
        
        self.y_i = y_i
        self.f_i = f_i

        self.y_i_plus_1 = y_i_plus_1
        self.f_i_plus_1 = f_i_plus_1

    def eval(self, x) -> float:
        t = (x - self.x_i) / self.h_i
        return (
              h00(t) * self.y_i        + self.h_i * h10(t) * self.f_i
            + h01(t) * self.y_i_plus_1 + self.h_i * h11(t) * self.f_i_plus_1
        )
    
    def prime(self, x) -> float:
        t = (x - self.x_i) / self.h_i
        return (
              h00_prime(t) * self.y_i        / self.h_i + h10_prime(t) * self.f_i
            + h01_prime(t) * self.y_i_plus_1 / self.h_i + h11_prime(t) * self.f_i_plus_1
        )
    
    def eval_horner(self, x) -> float:
        t = (x - self.x_i) / self.h_i
        return (
              h00_horner(t) * self.y_i        + self.h_i * h10_horner(t) * self.f_i
            + h01_horner(t) * self.y_i_plus_1 + self.h_i * h11_horner(t) * self.f_i_plus_1
        )
    
    def prime_horner(self, x) -> float:
        t = (x - self.x_i) / self.h_i
        return (
              h00_prime_horner(t) * self.y_i        / self.h_i + h10_prime_horner(t) * self.f_i
            + h01_prime_horner(t) * self.y_i_plus_1 / self.h_i + h11_prime_horner(t) * self.f_i_plus_1
        )

def create_defect_samplings(res, fn_s) -> Tuple[float, float, HB]:
    result = []
    for i in range(len(res) - 1):  
        x_i, y_i                 = res[i]    
        x_i_plus_1, y_i_plus_1   = res[i + 1]
  
        f_i         = fn_s[i]    
        f_i_plus_1  = fn_s[i + 1]
        
        interp = HB (
                x_i, x_i_plus_1,
                y_i, f_i,
                y_i_plus_1, f_i_plus_1 
        )
        result.append( (x_i, x_i_plus_1, interp) )
    return result

class ContinuousSolution:
    def __init__(self) -> None:
        self.interps = []
    
    def eval(self, x) -> float:
        for hb in self.interps:
            if (hb.x_i <= x <= hb.x_i_plus_1):
                return hb.eval(x)

        print(f"ERROR: {x} is outside of the solution range: {self.interps[0].x_i} <= x <= {self.interps[-1].x_i_plus_1}")
        return -1

    def prime(self, x) -> float:
        for hb in self.interps:
            if (hb.x_i <= x <= hb.x_i_plus_1):
                return hb.prime(x)

        print(f"ERROR: {x} is outside of the solution range: {self.interps[0].x_i} <= x <= {self.interps[-1].x_i_plus_1}")
        return -1
    
    def append(self, interp) -> None:
        self.interps.append(interp)
    
    def extend(self, newInterps) -> None:
        self.interps.extend(newInterps)

def create_continuous_sol_from_results(res, fn_s) -> ContinuousSolution:
    interps = []
    
    for i in range(len(res) - 1): 
        x_i, y_i                 = res[i]    
        x_i_plus_1, y_i_plus_1   = res[i + 1]
  
        f_i         = fn_s[i]    
        f_i_plus_1  = fn_s[i + 1]
        
        interps.append(
            HB (
                x_i, x_i_plus_1,
                y_i, f_i,
                y_i_plus_1, f_i_plus_1 
            )
        )
    continuous_sol = ContinuousSolution()
    continuous_sol.extend(interps)
    return continuous_sol

