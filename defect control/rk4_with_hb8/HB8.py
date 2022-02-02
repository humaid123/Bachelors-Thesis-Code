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