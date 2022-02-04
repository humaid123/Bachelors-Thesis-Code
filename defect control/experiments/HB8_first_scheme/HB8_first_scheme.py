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


def d0(x, alpha, beta) -> float:
    return ( (2*alpha**2 + 8*alpha*beta + 2*alpha + 6*beta**2 + 4*beta)/(alpha**6*beta**3 + 6*alpha**5*beta**4 + 3*alpha**5*beta**3 + 15*alpha**4*beta**5 + 15*alpha**4*beta**4 + 3*alpha**4*beta**3 + 20*alpha**3*beta**6 + 30*alpha**3*beta**5 + 12*alpha**3*beta**4 + alpha**3*beta**3 + 15*alpha**2*beta**7 + 30*alpha**2*beta**6 + 18*alpha**2*beta**5 + 3*alpha**2*beta**4 + 6*alpha*beta**8 + 15*alpha*beta**7 + 12*alpha*beta**6 + 3*alpha*beta**5 + beta**9 + 3*beta**8 + 3*beta**7 + beta**6)*(x**7) + (6*alpha**3 + 27*alpha**2*beta + 2*alpha**2 + 28*alpha*beta**2 - alpha*beta - 4*alpha + 7*beta**3 - 7*beta**2 - 8*beta)/(alpha**6*beta**3 + 6*alpha**5*beta**4 + 3*alpha**5*beta**3 + 15*alpha**4*beta**5 + 15*alpha**4*beta**4 + 3*alpha**4*beta**3 + 20*alpha**3*beta**6 + 30*alpha**3*beta**5 + 12*alpha**3*beta**4 + alpha**3*beta**3 + 15*alpha**2*beta**7 + 30*alpha**2*beta**6 + 18*alpha**2*beta**5 + 3*alpha**2*beta**4 + 6*alpha*beta**8 + 15*alpha*beta**7 + 12*alpha*beta**6 + 3*alpha*beta**5 + beta**9 + 3*beta**8 + 3*beta**7 + beta**6)*(x**6) + (6*alpha**4 + 30*alpha**3*beta - 6*alpha**3 + 38*alpha**2*beta**2 - 36*alpha**2*beta - 10*alpha**2 + 14*alpha*beta**3 - 46*alpha*beta**2 - 22*alpha*beta + 2*alpha - 14*beta**3 - 4*beta**2 + 4*beta)/(alpha**6*beta**3 + 6*alpha**5*beta**4 + 3*alpha**5*beta**3 + 15*alpha**4*beta**5 + 15*alpha**4*beta**4 + 3*alpha**4*beta**3 + 20*alpha**3*beta**6 + 30*alpha**3*beta**5 + 12*alpha**3*beta**4 + alpha**3*beta**3 + 15*alpha**2*beta**7 + 30*alpha**2*beta**6 + 18*alpha**2*beta**5 + 3*alpha**2*beta**4 + 6*alpha*beta**8 + 15*alpha*beta**7 + 12*alpha*beta**6 + 3*alpha*beta**5 + beta**9 + 3*beta**8 + 3*beta**7 + beta**6)*(x**5) + (2*alpha**5 + 11*alpha**4*beta - 10*alpha**4 + 16*alpha**3*beta**2 - 53*alpha**3*beta - 6*alpha**3 + 7*alpha**2*beta**3 - 71*alpha**2*beta**2 - 9*alpha**2*beta + 6*alpha**2 - 28*alpha*beta**3 + 8*alpha*beta**2 + 15*alpha*beta + 7*beta**3 + 5*beta**2)/(alpha**6*beta**3 + 6*alpha**5*beta**4 + 3*alpha**5*beta**3 + 15*alpha**4*beta**5 + 15*alpha**4*beta**4 + 3*alpha**4*beta**3 + 20*alpha**3*beta**6 + 30*alpha**3*beta**5 + 12*alpha**3*beta**4 + alpha**3*beta**3 + 15*alpha**2*beta**7 + 30*alpha**2*beta**6 + 18*alpha**2*beta**5 + 3*alpha**2*beta**4 + 6*alpha*beta**8 + 15*alpha*beta**7 + 12*alpha*beta**6 + 3*alpha*beta**5 + beta**9 + 3*beta**8 + 3*beta**7 + beta**6)*(x**4) + (-4*alpha**5 - 22*alpha**4*beta + 2*alpha**4 - 32*alpha**3*beta**2 + 16*alpha**3*beta + 6*alpha**3 - 14*alpha**2*beta**3 + 28*alpha**2*beta**2 + 18*alpha**2*beta + 14*alpha*beta**3 + 10*alpha*beta**2)/(alpha**6*beta**3 + 6*alpha**5*beta**4 + 3*alpha**5*beta**3 + 15*alpha**4*beta**5 + 15*alpha**4*beta**4 + 3*alpha**4*beta**3 + 20*alpha**3*beta**6 + 30*alpha**3*beta**5 + 12*alpha**3*beta**4 + alpha**3*beta**3 + 15*alpha**2*beta**7 + 30*alpha**2*beta**6 + 18*alpha**2*beta**5 + 3*alpha**2*beta**4 + 6*alpha*beta**8 + 15*alpha*beta**7 + 12*alpha*beta**6 + 3*alpha*beta**5 + beta**9 + 3*beta**8 + 3*beta**7 + beta**6)*(x**3) + (2*alpha**4 + 9*alpha**3*beta + 2*alpha**3 + 7*alpha**2*beta**2 + 5*alpha**2*beta)/(alpha**5*beta**3 + 5*alpha**4*beta**4 + 3*alpha**4*beta**3 + 10*alpha**3*beta**5 + 12*alpha**3*beta**4 + 3*alpha**3*beta**3 + 10*alpha**2*beta**6 + 18*alpha**2*beta**5 + 9*alpha**2*beta**4 + alpha**2*beta**3 + 5*alpha*beta**7 + 12*alpha*beta**6 + 9*alpha*beta**5 + 2*alpha*beta**4 + beta**8 + 3*beta**7 + 3*beta**6 + beta**5)*(x**2) + 0*(x) + 0 )
def d1(x, alpha, beta) -> float:
    return ( 1/(alpha**4*beta**2 + 4*alpha**3*beta**3 + 2*alpha**3*beta**2 + 6*alpha**2*beta**4 + 6*alpha**2*beta**3 + alpha**2*beta**2 + 4*alpha*beta**5 + 6*alpha*beta**4 + 2*alpha*beta**3 + beta**6 + 2*beta**5 + beta**4)*(x**7) + (3*alpha + beta - 2)/(alpha**4*beta**2 + 4*alpha**3*beta**3 + 2*alpha**3*beta**2 + 6*alpha**2*beta**4 + 6*alpha**2*beta**3 + alpha**2*beta**2 + 4*alpha*beta**5 + 6*alpha*beta**4 + 2*alpha*beta**3 + beta**6 + 2*beta**5 + beta**4)*(x**6) + (3*alpha**2 + 2*alpha*beta - 6*alpha - 2*beta + 1)/(alpha**4*beta**2 + 4*alpha**3*beta**3 + 2*alpha**3*beta**2 + 6*alpha**2*beta**4 + 6*alpha**2*beta**3 + alpha**2*beta**2 + 4*alpha*beta**5 + 6*alpha*beta**4 + 2*alpha*beta**3 + beta**6 + 2*beta**5 + beta**4)*(x**5) + (alpha**3 + alpha**2*beta - 6*alpha**2 - 4*alpha*beta + 3*alpha + beta)/(alpha**4*beta**2 + 4*alpha**3*beta**3 + 2*alpha**3*beta**2 + 6*alpha**2*beta**4 + 6*alpha**2*beta**3 + alpha**2*beta**2 + 4*alpha*beta**5 + 6*alpha*beta**4 + 2*alpha*beta**3 + beta**6 + 2*beta**5 + beta**4)*(x**4) + (-2*alpha**3 - 2*alpha**2*beta + 3*alpha**2 + 2*alpha*beta)/(alpha**4*beta**2 + 4*alpha**3*beta**3 + 2*alpha**3*beta**2 + 6*alpha**2*beta**4 + 6*alpha**2*beta**3 + alpha**2*beta**2 + 4*alpha*beta**5 + 6*alpha*beta**4 + 2*alpha*beta**3 + beta**6 + 
2*beta**5 + beta**4)*(x**3) + alpha**2/(alpha**3*beta**2 + 3*alpha**2*beta**3 + 2*alpha**2*beta**2 + 3*alpha*beta**4 + 4*alpha*beta**3 + alpha*beta**2 + beta**5 + 2*beta**4 + beta**3)*(x**2) + 0*(x) + 0 )
def d2(x, alpha, beta) -> float:
    return ( (-2*alpha**2 + 4*alpha*beta - 2*alpha + 2*beta)/(alpha**6*beta**3 + 3*alpha**5*beta**3 + 3*alpha**4*beta**3 + alpha**3*beta**3)*(x**7) + (-6*alpha**3 + 9*alpha**2*beta - 2*alpha**2 + 8*alpha*beta**2 - 5*alpha*beta + 4*alpha + 4*beta**2 - 4*beta)/(alpha**6*beta**3 + 3*alpha**5*beta**3 + 3*alpha**4*beta**3 + alpha**3*beta**3)*(x**6) + (-6*alpha**4 + 6*alpha**3*beta + 6*alpha**3 + 16*alpha**2*beta**2 - 18*alpha**2*beta + 10*alpha**2 + 4*alpha*beta**3 - 8*alpha*beta**2 - 2*alpha*beta - 2*alpha + 2*beta**3 - 8*beta**2 + 2*beta)/(alpha**6*beta**3 + 3*alpha**5*beta**3 + 3*alpha**4*beta**3 + alpha**3*beta**3)*(x**5) + (-2*alpha**5 + alpha**4*beta + 10*alpha**4 + 8*alpha**3*beta**2 - 13*alpha**3*beta + 6*alpha**3 + 5*alpha**2*beta**3 - 28*alpha**2*beta**2 + 9*alpha**2*beta - 6*alpha**2 - 5*alpha*beta**3 - 8*alpha*beta**2 + 3*alpha*beta - 4*beta**3 + 4*beta**2)/(alpha**6*beta**3 + 3*alpha**5*beta**3 + 3*alpha**4*beta**3 + alpha**3*beta**3)*(x**4) + (4*alpha**5 - 2*alpha**4*beta - 2*alpha**4 - 16*alpha**3*beta**2 + 8*alpha**3*beta - 6*alpha**3 - 10*alpha**2*beta**3 + 8*alpha**2*beta**2 - 2*alpha*beta**3 + 8*alpha*beta**2 + 2*beta**3)/(alpha**6*beta**3 + 3*alpha**5*beta**3 + 3*alpha**4*beta**3 + alpha**3*beta**3)*(x**3) + (-2*alpha**4 + alpha**3*beta - 
2*alpha**3 + 8*alpha**2*beta**2 - alpha**2*beta + 5*alpha*beta**3 + 4*alpha*beta**2 + 3*beta**3)/(alpha**5*beta**3 + 3*alpha**4*beta**3 + 3*alpha**3*beta**3 + alpha**2*beta**3)*(x**2) + 0*(x) + 0 )
def d3(x, alpha, beta) -> float:
    return ( 1/(alpha**4*beta**2 + 2*alpha**3*beta**2 + alpha**2*beta**2)*(x**7) + (3*alpha + 2*beta - 2)/(alpha**4*beta**2 + 2*alpha**3*beta**2 + alpha**2*beta**2)*(x**6) + (3*alpha**2 + 4*alpha*beta - 6*alpha + beta**2 - 4*beta + 1)/(alpha**4*beta**2 + 2*alpha**3*beta**2 + alpha**2*beta**2)*(x**5) + (alpha**3 + 2*alpha**2*beta - 6*alpha**2 + alpha*beta**2 - 8*alpha*beta + 3*alpha - 2*beta**2 + 2*beta)/(alpha**4*beta**2 + 2*alpha**3*beta**2 + alpha**2*beta**2)*(x**4) + (-2*alpha**3 - 4*alpha**2*beta + 3*alpha**2 - 2*alpha*beta**2 + 4*alpha*beta + beta**2)/(alpha**4*beta**2 + 2*alpha**3*beta**2 + alpha**2*beta**2)*(x**3) + (alpha**2 + 2*alpha*beta + beta**2)/(alpha**3*beta**2 + 2*alpha**2*beta**2 + alpha*beta**2)*(x**2) + 0*(x) + 0 )
def d4(x, alpha, beta) -> float:
    return ( (2*alpha**2 + 2*alpha*beta - 4*alpha - 2*beta)/(alpha**6 + 3*alpha**5*beta + 3*alpha**4*beta**2 + alpha**3*beta**3)*(x**7) + (8*alpha**3 + 12*alpha**2*beta - 19*alpha**2 + 4*alpha*beta**2 - 19*alpha*beta + 8*alpha - 4*beta**2 + 4*beta)/(alpha**6 + 3*alpha**5*beta + 3*alpha**4*beta**2 + alpha**3*beta**3)*(x**6) + (12*alpha**4 + 24*alpha**3*beta - 36*alpha**3 + 14*alpha**2*beta**2 - 54*alpha**2*beta + 32*alpha**2 + 2*alpha*beta**3 - 22*alpha*beta**2 + 32*alpha*beta - 4*alpha - 2*beta**3 + 8*beta**2 - 2*beta)/(alpha**6 + 3*alpha**5*beta + 3*alpha**4*beta**2 + alpha**3*beta**3)*(x**5) + (8*alpha**5 + 20*alpha**4*beta - 34*alpha**4 + 16*alpha**3*beta**2 - 68*alpha**3*beta + 48*alpha**3 + 4*alpha**2*beta**3 - 41*alpha**2*beta**2 + 72*alpha**2*beta - 15*alpha**2 - 7*alpha*beta**3 + 32*alpha*beta**2 - 15*alpha*beta + 4*beta**3 - 4*beta**2)/(alpha**6 + 3*alpha**5*beta + 3*alpha**4*beta**2 + alpha**3*beta**3)*(x**4) + (2*alpha**6 + 6*alpha**5*beta - 16*alpha**5 + 6*alpha**4*beta**2 - 40*alpha**4*beta + 32*alpha**4 + 2*alpha**3*beta**3 - 32*alpha**3*beta**2 + 64*alpha**3*beta - 20*alpha**3 - 8*alpha**2*beta**3 + 40*alpha**2*beta**2 - 30*alpha**2*beta + 8*alpha*beta**3 - 14*alpha*beta**2 - 2*beta**3)/(alpha**6 + 3*alpha**5*beta + 3*alpha**4*beta**2 + alpha**3*beta**3)*(x**3) + (-3*alpha**4 - 6*alpha**3*beta + 8*alpha**3 - 3*alpha**2*beta**2 + 12*alpha**2*beta - 10*alpha**2 + 4*alpha*beta**2 - 10*alpha*beta - 3*beta**2)/(alpha**4 + 
2*alpha**3*beta + alpha**2*beta**2)*(x**2) + 0*(x) + 1 )
def d5(x, alpha, beta) -> float:
    return ( 1/(alpha**4 + 2*alpha**3*beta + alpha**2*beta**2)*(x**7) + (4*alpha + 2*beta - 2)/(alpha**4 + 2*alpha**3*beta + alpha**2*beta**2)*(x**6) + (6*alpha**2 + 6*alpha*beta - 8*alpha + beta**2 - 4*beta + 1)/(alpha**4 + 2*alpha**3*beta + alpha**2*beta**2)*(x**5) + (4*alpha**3 + 6*alpha**2*beta - 12*alpha**2 + 2*alpha*beta**2 - 12*alpha*beta + 4*alpha - 2*beta**2 + 2*beta)/(alpha**4 + 2*alpha**3*beta + alpha**2*beta**2)*(x**4) + (alpha**4 + 2*alpha**3*beta - 8*alpha**3 + alpha**2*beta**2 - 12*alpha**2*beta + 6*alpha**2 - 4*alpha*beta**2 + 6*alpha*beta + beta**2)/(alpha**4 + 2*alpha**3*beta + alpha**2*beta**2)*(x**3) + (-2*alpha**2 - 2*alpha*beta + 4*alpha + 2*beta)/(alpha**2 + alpha*beta)*(x**2) + 1*(x) + 0 )      
def d6(x, alpha, beta) -> float:
    return ( (-2*alpha**2 - 2*alpha*beta - 8*alpha - 4*beta - 6)/(alpha**6 + 3*alpha**5*beta + 6*alpha**5 + 3*alpha**4*beta**2 + 15*alpha**4*beta + 15*alpha**4 + alpha**3*beta**3 + 12*alpha**3*beta**2 + 30*alpha**3*beta + 20*alpha**3 + 3*alpha**2*beta**3 + 18*alpha**2*beta**2 + 30*alpha**2*beta + 15*alpha**2 + 3*alpha*beta**3 + 12*alpha*beta**2 + 15*alpha*beta + 6*alpha + beta**3 + 3*beta**2 + 3*beta + 1)*(x**7) + (-8*alpha**3 - 12*alpha**2*beta - 29*alpha**2 - 4*alpha*beta**2 - 29*alpha*beta - 14*alpha - 8*beta**2 - 7*beta + 7)/(alpha**6 + 3*alpha**5*beta + 6*alpha**5 + 3*alpha**4*beta**2 + 15*alpha**4*beta + 15*alpha**4 + alpha**3*beta**3 + 12*alpha**3*beta**2 + 30*alpha**3*beta + 20*alpha**3 + 3*alpha**2*beta**3 + 18*alpha**2*beta**2 + 30*alpha**2*beta + 15*alpha**2 + 3*alpha*beta**3 + 12*alpha*beta**2 + 15*alpha*beta + 6*alpha + beta**3 + 3*beta**2 + 3*beta + 1)*(x**6) + (-12*alpha**4 - 24*alpha**3*beta - 36*alpha**3 - 14*alpha**2*beta**2 - 54*alpha**2*beta + 4*alpha**2 - 2*alpha*beta**3 - 26*alpha*beta**2 + 4*alpha*beta + 28*alpha - 4*beta**3 + 4*beta**2 + 14*beta)/(alpha**6 + 3*alpha**5*beta + 6*alpha**5 + 3*alpha**4*beta**2 + 15*alpha**4*beta + 15*alpha**4 + alpha**3*beta**3 + 12*alpha**3*beta**2 + 30*alpha**3*beta + 20*alpha**3 + 3*alpha**2*beta**3 + 18*alpha**2*beta**2 + 30*alpha**2*beta + 15*alpha**2 + 3*alpha*beta**3 + 12*alpha*beta**2 + 15*alpha*beta + 6*alpha + beta**3 + 3*beta**2 + 3*beta + 1)*(x**5) + (-8*alpha**5 - 20*alpha**4*beta - 14*alpha**4 - 16*alpha**3*beta**2 - 28*alpha**3*beta + 36*alpha**3 - 4*alpha**2*beta**3 - 19*alpha**2*beta**2 + 54*alpha**2*beta + 42*alpha**2 - 5*alpha*beta**3 + 28*alpha*beta**2 + 42*alpha*beta + 5*beta**3 + 7*beta**2)/(alpha**6 + 3*alpha**5*beta + 6*alpha**5 + 3*alpha**4*beta**2 + 15*alpha**4*beta + 15*alpha**4 + alpha**3*beta**3 + 12*alpha**3*beta**2 + 30*alpha**3*beta + 20*alpha**3 + 3*alpha**2*beta**3 + 18*alpha**2*beta**2 + 30*alpha**2*beta + 15*alpha**2 + 3*alpha*beta**3 + 12*alpha*beta**2 + 15*alpha*beta + 6*alpha + beta**3 + 3*beta**2 + 3*beta + 1)*(x**4) + (-2*alpha**6 - 6*alpha**5*beta + 4*alpha**5 - 6*alpha**4*beta**2 + 10*alpha**4*beta + 34*alpha**4 - 2*alpha**3*beta**3 + 8*alpha**3*beta**2 + 68*alpha**3*beta + 28*alpha**3 + 2*alpha**2*beta**3 + 44*alpha**2*beta**2 + 42*alpha**2*beta + 10*alpha*beta**3 + 14*alpha*beta**2)/(alpha**6 + 3*alpha**5*beta + 6*alpha**5 + 3*alpha**4*beta**2 + 15*alpha**4*beta + 15*alpha**4 + alpha**3*beta**3 + 12*alpha**3*beta**2 + 30*alpha**3*beta + 20*alpha**3 + 3*alpha**2*beta**3 + 18*alpha**2*beta**2 + 30*alpha**2*beta + 15*alpha**2 + 3*alpha*beta**3 + 12*alpha*beta**2 + 15*alpha*beta + 6*alpha + beta**3 + 3*beta**2 + 3*beta + 1)*(x**3) + (3*alpha**6 + 9*alpha**5*beta + 10*alpha**5 + 9*alpha**4*beta**2 + 25*alpha**4*beta + 7*alpha**4 + 3*alpha**3*beta**3 + 20*alpha**3*beta**2 + 14*alpha**3*beta + 5*alpha**2*beta**3 + 7*alpha**2*beta**2)/(alpha**6 + 3*alpha**5*beta + 6*alpha**5 + 3*alpha**4*beta**2 + 15*alpha**4*beta + 15*alpha**4 + alpha**3*beta**3 + 12*alpha**3*beta**2 + 30*alpha**3*beta + 20*alpha**3 + 3*alpha**2*beta**3 + 18*alpha**2*beta**2 + 30*alpha**2*beta + 15*alpha**2 + 3*alpha*beta**3 + 12*alpha*beta**2 + 15*alpha*beta + 6*alpha + beta**3 + 3*beta**2 + 3*beta + 1)*(x**2) + 0*(x) + 0 )
def d7(x, alpha, beta) -> float:
    return ( 1/(alpha**4 + 2*alpha**3*beta + 4*alpha**3 + alpha**2*beta**2 + 6*alpha**2*beta + 6*alpha**2 + 2*alpha*beta**2 + 6*alpha*beta + 4*alpha + beta**2 + 2*beta + 1)*(x**7) + 
(4*alpha + 2*beta - 1)/(alpha**4 + 2*alpha**3*beta + 4*alpha**3 + alpha**2*beta**2 + 6*alpha**2*beta + 6*alpha**2 + 2*alpha*beta**2 + 6*alpha*beta + 4*alpha + beta**2 + 2*beta + 1)*(x**6) + (6*alpha**2 + 6*alpha*beta - 4*alpha + beta**2 - 2*beta)/(alpha**4 + 2*alpha**3*beta + 4*alpha**3 + alpha**2*beta**2 + 6*alpha**2*beta + 6*alpha**2 + 2*alpha*beta**2 + 6*alpha*beta + 4*alpha + beta**2 + 2*beta + 1)*(x**5) + (4*alpha**3 + 6*alpha**2*beta - 6*alpha**2 + 2*alpha*beta**2 - 6*alpha*beta - beta**2)/(alpha**4 + 2*alpha**3*beta + 4*alpha**3 + alpha**2*beta**2 + 6*alpha**2*beta + 6*alpha**2 + 2*alpha*beta**2 + 6*alpha*beta + 4*alpha + beta**2 + 2*beta + 1)*(x**4) + (alpha**4 + 2*alpha**3*beta - 4*alpha**3 + alpha**2*beta**2 
- 6*alpha**2*beta - 2*alpha*beta**2)/(alpha**4 + 2*alpha**3*beta + 4*alpha**3 + alpha**2*beta**2 + 6*alpha**2*beta + 6*alpha**2 + 2*alpha*beta**2 + 6*alpha*beta + 4*alpha + beta**2 + 2*beta + 1)*(x**3) + (-alpha**4 - 2*alpha**3*beta - alpha**2*beta**2)/(alpha**4 + 2*alpha**3*beta + 4*alpha**3 + alpha**2*beta**2 + 6*alpha**2*beta + 6*alpha**2 + 2*alpha*beta**2 + 
6*alpha*beta + 4*alpha + beta**2 + 2*beta + 1)*(x**2) + 0*(x) + 0 )
def d0_prime(x, alpha, beta) -> float:
    return ( 7*(2*alpha**2 + 8*alpha*beta + 2*alpha + 6*beta**2 + 4*beta)/(alpha**6*beta**3 + 6*alpha**5*beta**4 + 3*alpha**5*beta**3 + 15*alpha**4*beta**5 + 15*alpha**4*beta**4 + 3*alpha**4*beta**3 + 20*alpha**3*beta**6 + 30*alpha**3*beta**5 + 12*alpha**3*beta**4 + alpha**3*beta**3 + 15*alpha**2*beta**7 + 30*alpha**2*beta**6 + 18*alpha**2*beta**5 + 3*alpha**2*beta**4 + 6*alpha*beta**8 + 15*alpha*beta**7 + 12*alpha*beta**6 + 3*alpha*beta**5 + beta**9 + 3*beta**8 + 3*beta**7 + beta**6)*(x**6) + 6*(6*alpha**3 + 27*alpha**2*beta + 2*alpha**2 + 28*alpha*beta**2 - alpha*beta - 4*alpha + 7*beta**3 - 7*beta**2 - 8*beta)/(alpha**6*beta**3 + 6*alpha**5*beta**4 + 3*alpha**5*beta**3 + 15*alpha**4*beta**5 + 15*alpha**4*beta**4 + 3*alpha**4*beta**3 + 20*alpha**3*beta**6 + 30*alpha**3*beta**5 + 12*alpha**3*beta**4 + alpha**3*beta**3 + 15*alpha**2*beta**7 + 30*alpha**2*beta**6 + 18*alpha**2*beta**5 + 3*alpha**2*beta**4 + 6*alpha*beta**8 + 15*alpha*beta**7 + 12*alpha*beta**6 + 3*alpha*beta**5 + beta**9 + 3*beta**8 + 3*beta**7 + beta**6)*(x**5) + 5*(6*alpha**4 + 30*alpha**3*beta - 6*alpha**3 
+ 38*alpha**2*beta**2 - 36*alpha**2*beta - 10*alpha**2 + 14*alpha*beta**3 - 46*alpha*beta**2 - 22*alpha*beta + 2*alpha - 14*beta**3 - 4*beta**2 + 4*beta)/(alpha**6*beta**3 + 6*alpha**5*beta**4 + 3*alpha**5*beta**3 + 15*alpha**4*beta**5 + 15*alpha**4*beta**4 + 3*alpha**4*beta**3 + 20*alpha**3*beta**6 + 30*alpha**3*beta**5 + 12*alpha**3*beta**4 + alpha**3*beta**3 
+ 15*alpha**2*beta**7 + 30*alpha**2*beta**6 + 18*alpha**2*beta**5 + 3*alpha**2*beta**4 + 6*alpha*beta**8 + 15*alpha*beta**7 + 12*alpha*beta**6 + 3*alpha*beta**5 + beta**9 + 3*beta**8 + 3*beta**7 + beta**6)*(x**4) + 4*(2*alpha**5 + 11*alpha**4*beta - 10*alpha**4 + 16*alpha**3*beta**2 - 53*alpha**3*beta - 6*alpha**3 + 7*alpha**2*beta**3 - 71*alpha**2*beta**2 - 9*alpha**2*beta + 6*alpha**2 - 28*alpha*beta**3 + 8*alpha*beta**2 + 15*alpha*beta + 7*beta**3 + 5*beta**2)/(alpha**6*beta**3 + 6*alpha**5*beta**4 + 3*alpha**5*beta**3 + 15*alpha**4*beta**5 + 15*alpha**4*beta**4 + 3*alpha**4*beta**3 + 20*alpha**3*beta**6 + 30*alpha**3*beta**5 + 12*alpha**3*beta**4 + alpha**3*beta**3 + 15*alpha**2*beta**7 + 30*alpha**2*beta**6 + 18*alpha**2*beta**5 + 3*alpha**2*beta**4 + 6*alpha*beta**8 + 15*alpha*beta**7 + 12*alpha*beta**6 + 3*alpha*beta**5 + beta**9 + 3*beta**8 + 3*beta**7 + beta**6)*(x**3) + 3*(-4*alpha**5 - 
22*alpha**4*beta + 2*alpha**4 - 32*alpha**3*beta**2 + 16*alpha**3*beta + 6*alpha**3 - 14*alpha**2*beta**3 + 28*alpha**2*beta**2 + 18*alpha**2*beta + 14*alpha*beta**3 + 10*alpha*beta**2)/(alpha**6*beta**3 + 6*alpha**5*beta**4 + 3*alpha**5*beta**3 + 15*alpha**4*beta**5 + 15*alpha**4*beta**4 + 3*alpha**4*beta**3 + 20*alpha**3*beta**6 + 30*alpha**3*beta**5 + 12*alpha**3*beta**4 + alpha**3*beta**3 + 15*alpha**2*beta**7 + 30*alpha**2*beta**6 + 18*alpha**2*beta**5 + 3*alpha**2*beta**4 + 6*alpha*beta**8 + 15*alpha*beta**7 + 12*alpha*beta**6 + 3*alpha*beta**5 + beta**9 + 3*beta**8 + 3*beta**7 + beta**6)*(x**2) + 2*(2*alpha**4 + 9*alpha**3*beta + 2*alpha**3 + 7*alpha**2*beta**2 + 5*alpha**2*beta)/(alpha**5*beta**3 + 5*alpha**4*beta**4 + 3*alpha**4*beta**3 + 10*alpha**3*beta**5 + 12*alpha**3*beta**4 + 3*alpha**3*beta**3 + 10*alpha**2*beta**6 + 18*alpha**2*beta**5 + 9*alpha**2*beta**4 + alpha**2*beta**3 + 5*alpha*beta**7 + 12*alpha*beta**6 + 9*alpha*beta**5 + 2*alpha*beta**4 + beta**8 + 3*beta**7 + 3*beta**6 + beta**5)*x + 0 )
def d1_prime(x, alpha, beta) -> float:
    return ( 7*1/(alpha**4*beta**2 + 4*alpha**3*beta**3 + 2*alpha**3*beta**2 + 6*alpha**2*beta**4 + 6*alpha**2*beta**3 + alpha**2*beta**2 + 4*alpha*beta**5 + 6*alpha*beta**4 + 2*alpha*beta**3 + beta**6 + 2*beta**5 + beta**4)*(x**6) + 6*(3*alpha + beta - 2)/(alpha**4*beta**2 + 4*alpha**3*beta**3 + 2*alpha**3*beta**2 + 6*alpha**2*beta**4 + 6*alpha**2*beta**3 + alpha**2*beta**2 + 4*alpha*beta**5 + 6*alpha*beta**4 + 2*alpha*beta**3 + beta**6 + 2*beta**5 + beta**4)*(x**5) + 5*(3*alpha**2 + 2*alpha*beta - 6*alpha - 2*beta + 1)/(alpha**4*beta**2 + 4*alpha**3*beta**3 + 2*alpha**3*beta**2 + 6*alpha**2*beta**4 + 6*alpha**2*beta**3 + alpha**2*beta**2 + 4*alpha*beta**5 + 6*alpha*beta**4 + 2*alpha*beta**3 + beta**6 + 2*beta**5 + beta**4)*(x**4) + 4*(alpha**3 + alpha**2*beta - 6*alpha**2 - 4*alpha*beta + 3*alpha + beta)/(alpha**4*beta**2 + 4*alpha**3*beta**3 + 2*alpha**3*beta**2 + 6*alpha**2*beta**4 + 6*alpha**2*beta**3 + alpha**2*beta**2 + 4*alpha*beta**5 + 6*alpha*beta**4 + 2*alpha*beta**3 + beta**6 + 2*beta**5 + beta**4)*(x**3) + 3*(-2*alpha**3 - 2*alpha**2*beta + 3*alpha**2 + 2*alpha*beta)/(alpha**4*beta**2 + 4*alpha**3*beta**3 + 2*alpha**3*beta**2 + 6*alpha**2*beta**4 + 6*alpha**2*beta**3 + alpha**2*beta**2 + 4*alpha*beta**5 + 6*alpha*beta**4 + 2*alpha*beta**3 + 
beta**6 + 2*beta**5 + beta**4)*(x**2) + 2*alpha**2/(alpha**3*beta**2 + 3*alpha**2*beta**3 + 2*alpha**2*beta**2 + 3*alpha*beta**4 + 4*alpha*beta**3 + alpha*beta**2 + beta**5 + 2*beta**4 + beta**3)*x + 0 )
def d2_prime(x, alpha, beta) -> float:
    return ( 7*(-2*alpha**2 + 4*alpha*beta - 2*alpha + 2*beta)/(alpha**6*beta**3 + 3*alpha**5*beta**3 + 3*alpha**4*beta**3 + alpha**3*beta**3)*(x**6) + 6*(-6*alpha**3 + 9*alpha**2*beta - 2*alpha**2 + 8*alpha*beta**2 - 5*alpha*beta + 4*alpha + 4*beta**2 - 4*beta)/(alpha**6*beta**3 + 3*alpha**5*beta**3 + 3*alpha**4*beta**3 + alpha**3*beta**3)*(x**5) + 5*(-6*alpha**4 + 6*alpha**3*beta + 6*alpha**3 + 16*alpha**2*beta**2 - 18*alpha**2*beta + 10*alpha**2 + 4*alpha*beta**3 - 8*alpha*beta**2 - 2*alpha*beta - 2*alpha + 2*beta**3 - 8*beta**2 + 2*beta)/(alpha**6*beta**3 + 3*alpha**5*beta**3 + 3*alpha**4*beta**3 + alpha**3*beta**3)*(x**4) + 4*(-2*alpha**5 + alpha**4*beta + 10*alpha**4 + 8*alpha**3*beta**2 - 13*alpha**3*beta + 6*alpha**3 + 5*alpha**2*beta**3 - 28*alpha**2*beta**2 + 9*alpha**2*beta - 6*alpha**2 - 5*alpha*beta**3 - 8*alpha*beta**2 + 3*alpha*beta - 4*beta**3 + 4*beta**2)/(alpha**6*beta**3 + 3*alpha**5*beta**3 + 3*alpha**4*beta**3 + alpha**3*beta**3)*(x**3) + 3*(4*alpha**5 - 2*alpha**4*beta - 2*alpha**4 - 16*alpha**3*beta**2 + 8*alpha**3*beta - 6*alpha**3 - 10*alpha**2*beta**3 + 8*alpha**2*beta**2 - 2*alpha*beta**3 + 8*alpha*beta**2 + 2*beta**3)/(alpha**6*beta**3 + 3*alpha**5*beta**3 + 3*alpha**4*beta**3 + alpha**3*beta**3)*(x**2) + 2*(-2*alpha**4 + alpha**3*beta - 2*alpha**3 + 8*alpha**2*beta**2 - alpha**2*beta + 5*alpha*beta**3 + 4*alpha*beta**2 + 3*beta**3)/(alpha**5*beta**3 + 3*alpha**4*beta**3 + 3*alpha**3*beta**3 + alpha**2*beta**3)*x + 0 )
def d3_prime(x, alpha, beta) -> float:
    return ( 7*1/(alpha**4*beta**2 + 2*alpha**3*beta**2 + alpha**2*beta**2)*(x**6) + 6*(3*alpha + 2*beta - 2)/(alpha**4*beta**2 + 2*alpha**3*beta**2 + alpha**2*beta**2)*(x**5) + 5*(3*alpha**2 + 4*alpha*beta - 6*alpha + beta**2 - 4*beta + 1)/(alpha**4*beta**2 + 2*alpha**3*beta**2 + alpha**2*beta**2)*(x**4) + 4*(alpha**3 + 2*alpha**2*beta - 6*alpha**2 + alpha*beta**2 - 8*alpha*beta + 3*alpha - 2*beta**2 + 2*beta)/(alpha**4*beta**2 + 2*alpha**3*beta**2 + alpha**2*beta**2)*(x**3) + 3*(-2*alpha**3 - 4*alpha**2*beta + 3*alpha**2 - 2*alpha*beta**2 + 4*alpha*beta + beta**2)/(alpha**4*beta**2 + 2*alpha**3*beta**2 + alpha**2*beta**2)*(x**2) + 2*(alpha**2 + 2*alpha*beta + beta**2)/(alpha**3*beta**2 + 2*alpha**2*beta**2 + alpha*beta**2)*x + 0 )
def d4_prime(x, alpha, beta) -> float:
    return ( 7*(2*alpha**2 + 2*alpha*beta - 4*alpha - 2*beta)/(alpha**6 + 3*alpha**5*beta + 3*alpha**4*beta**2 + alpha**3*beta**3)*(x**6) + 6*(8*alpha**3 + 12*alpha**2*beta - 19*alpha**2 + 4*alpha*beta**2 - 19*alpha*beta + 8*alpha - 4*beta**2 + 4*beta)/(alpha**6 + 3*alpha**5*beta + 3*alpha**4*beta**2 + alpha**3*beta**3)*(x**5) + 5*(12*alpha**4 + 24*alpha**3*beta - 36*alpha**3 + 14*alpha**2*beta**2 - 54*alpha**2*beta + 32*alpha**2 + 2*alpha*beta**3 - 22*alpha*beta**2 + 32*alpha*beta - 4*alpha - 2*beta**3 + 8*beta**2 - 2*beta)/(alpha**6 + 3*alpha**5*beta + 3*alpha**4*beta**2 + alpha**3*beta**3)*(x**4) + 4*(8*alpha**5 + 20*alpha**4*beta - 34*alpha**4 + 16*alpha**3*beta**2 - 68*alpha**3*beta + 48*alpha**3 + 4*alpha**2*beta**3 - 41*alpha**2*beta**2 + 72*alpha**2*beta - 15*alpha**2 - 7*alpha*beta**3 + 32*alpha*beta**2 - 15*alpha*beta + 4*beta**3 - 4*beta**2)/(alpha**6 + 3*alpha**5*beta + 3*alpha**4*beta**2 + alpha**3*beta**3)*(x**3) + 3*(2*alpha**6 + 6*alpha**5*beta - 16*alpha**5 + 6*alpha**4*beta**2 - 40*alpha**4*beta + 32*alpha**4 + 2*alpha**3*beta**3 - 32*alpha**3*beta**2 + 64*alpha**3*beta - 20*alpha**3 - 8*alpha**2*beta**3 + 40*alpha**2*beta**2 - 30*alpha**2*beta + 8*alpha*beta**3 - 14*alpha*beta**2 - 2*beta**3)/(alpha**6 + 3*alpha**5*beta + 3*alpha**4*beta**2 + alpha**3*beta**3)*(x**2) + 2*(-3*alpha**4 - 6*alpha**3*beta + 8*alpha**3 - 3*alpha**2*beta**2 + 12*alpha**2*beta - 10*alpha**2 + 4*alpha*beta**2 - 10*alpha*beta - 3*beta**2)/(alpha**4 + 2*alpha**3*beta + alpha**2*beta**2)*x + 0 )
def d5_prime(x, alpha, beta) -> float:
    return ( 7*1/(alpha**4 + 2*alpha**3*beta + alpha**2*beta**2)*(x**6) + 6*(4*alpha + 2*beta - 2)/(alpha**4 + 2*alpha**3*beta + alpha**2*beta**2)*(x**5) + 5*(6*alpha**2 + 6*alpha*beta - 8*alpha + beta**2 - 4*beta + 1)/(alpha**4 + 2*alpha**3*beta + alpha**2*beta**2)*(x**4) + 4*(4*alpha**3 + 6*alpha**2*beta - 12*alpha**2 + 2*alpha*beta**2 - 12*alpha*beta + 4*alpha - 2*beta**2 + 2*beta)/(alpha**4 + 2*alpha**3*beta + alpha**2*beta**2)*(x**3) + 3*(alpha**4 + 2*alpha**3*beta - 8*alpha**3 + alpha**2*beta**2 - 12*alpha**2*beta + 6*alpha**2 - 4*alpha*beta**2 + 6*alpha*beta + beta**2)/(alpha**4 + 2*alpha**3*beta + alpha**2*beta**2)*(x**2) + 2*(-2*alpha**2 - 2*alpha*beta + 4*alpha + 2*beta)/(alpha**2 + alpha*beta)*x + 1 )       
def d6_prime(x, alpha, beta) -> float:
    return ( 7*(-2*alpha**2 - 2*alpha*beta - 8*alpha - 4*beta - 6)/(alpha**6 + 3*alpha**5*beta + 6*alpha**5 + 3*alpha**4*beta**2 + 15*alpha**4*beta + 15*alpha**4 + alpha**3*beta**3 + 12*alpha**3*beta**2 + 30*alpha**3*beta + 20*alpha**3 + 3*alpha**2*beta**3 + 18*alpha**2*beta**2 + 30*alpha**2*beta + 15*alpha**2 + 3*alpha*beta**3 + 12*alpha*beta**2 + 15*alpha*beta + 6*alpha + beta**3 + 3*beta**2 + 3*beta + 1)*(x**6) + 6*(-8*alpha**3 - 12*alpha**2*beta - 29*alpha**2 - 4*alpha*beta**2 - 29*alpha*beta - 14*alpha - 8*beta**2 - 7*beta + 7)/(alpha**6 + 3*alpha**5*beta + 6*alpha**5 + 3*alpha**4*beta**2 + 15*alpha**4*beta + 15*alpha**4 + alpha**3*beta**3 + 12*alpha**3*beta**2 + 30*alpha**3*beta + 20*alpha**3 + 3*alpha**2*beta**3 + 18*alpha**2*beta**2 + 30*alpha**2*beta + 15*alpha**2 + 3*alpha*beta**3 + 12*alpha*beta**2 + 15*alpha*beta + 6*alpha + beta**3 + 3*beta**2 + 3*beta + 1)*(x**5) + 5*(-12*alpha**4 - 
24*alpha**3*beta - 36*alpha**3 - 14*alpha**2*beta**2 - 54*alpha**2*beta + 4*alpha**2 - 2*alpha*beta**3 - 26*alpha*beta**2 + 4*alpha*beta + 28*alpha - 4*beta**3 + 4*beta**2 + 14*beta)/(alpha**6 + 3*alpha**5*beta + 6*alpha**5 + 3*alpha**4*beta**2 + 15*alpha**4*beta + 15*alpha**4 + alpha**3*beta**3 + 12*alpha**3*beta**2 + 30*alpha**3*beta + 20*alpha**3 + 3*alpha**2*beta**3 + 18*alpha**2*beta**2 + 30*alpha**2*beta + 15*alpha**2 + 3*alpha*beta**3 + 12*alpha*beta**2 + 15*alpha*beta + 6*alpha + beta**3 + 3*beta**2 + 3*beta + 1)*(x**4) + 4*(-8*alpha**5 - 20*alpha**4*beta - 14*alpha**4 - 16*alpha**3*beta**2 - 28*alpha**3*beta + 36*alpha**3 - 4*alpha**2*beta**3 - 19*alpha**2*beta**2 + 54*alpha**2*beta + 42*alpha**2 - 5*alpha*beta**3 + 28*alpha*beta**2 + 42*alpha*beta + 5*beta**3 + 7*beta**2)/(alpha**6 + 3*alpha**5*beta + 6*alpha**5 + 3*alpha**4*beta**2 + 15*alpha**4*beta + 15*alpha**4 + alpha**3*beta**3 + 12*alpha**3*beta**2 + 30*alpha**3*beta + 20*alpha**3 + 3*alpha**2*beta**3 + 18*alpha**2*beta**2 + 30*alpha**2*beta + 15*alpha**2 + 3*alpha*beta**3 + 12*alpha*beta**2 + 15*alpha*beta + 6*alpha + beta**3 + 3*beta**2 + 3*beta + 1)*(x**3) + 3*(-2*alpha**6 - 6*alpha**5*beta + 4*alpha**5 - 6*alpha**4*beta**2 + 10*alpha**4*beta + 34*alpha**4 - 2*alpha**3*beta**3 + 8*alpha**3*beta**2 + 68*alpha**3*beta + 28*alpha**3 + 2*alpha**2*beta**3 + 44*alpha**2*beta**2 + 42*alpha**2*beta + 10*alpha*beta**3 + 14*alpha*beta**2)/(alpha**6 + 3*alpha**5*beta + 6*alpha**5 + 3*alpha**4*beta**2 + 15*alpha**4*beta + 15*alpha**4 + alpha**3*beta**3 + 12*alpha**3*beta**2 + 30*alpha**3*beta + 20*alpha**3 + 3*alpha**2*beta**3 + 18*alpha**2*beta**2 + 30*alpha**2*beta + 15*alpha**2 + 3*alpha*beta**3 + 12*alpha*beta**2 + 15*alpha*beta + 6*alpha + beta**3 + 3*beta**2 + 3*beta + 1)*(x**2) + 2*(3*alpha**6 + 9*alpha**5*beta + 10*alpha**5 + 9*alpha**4*beta**2 + 25*alpha**4*beta + 7*alpha**4 + 3*alpha**3*beta**3 + 20*alpha**3*beta**2 + 14*alpha**3*beta + 5*alpha**2*beta**3 + 7*alpha**2*beta**2)/(alpha**6 + 3*alpha**5*beta + 6*alpha**5 + 3*alpha**4*beta**2 + 15*alpha**4*beta + 15*alpha**4 + alpha**3*beta**3 + 12*alpha**3*beta**2 + 30*alpha**3*beta + 20*alpha**3 + 3*alpha**2*beta**3 + 18*alpha**2*beta**2 + 30*alpha**2*beta + 15*alpha**2 + 3*alpha*beta**3 + 12*alpha*beta**2 + 15*alpha*beta + 6*alpha + beta**3 + 3*beta**2 + 3*beta + 1)*x + 0 )
def d7_prime(x, alpha, beta) -> float:
    return ( 7*1/(alpha**4 + 2*alpha**3*beta + 4*alpha**3 + alpha**2*beta**2 + 6*alpha**2*beta + 6*alpha**2 + 2*alpha*beta**2 + 6*alpha*beta + 4*alpha + beta**2 + 2*beta + 1)*(x**6) 
+ 6*(4*alpha + 2*beta - 1)/(alpha**4 + 2*alpha**3*beta + 4*alpha**3 + alpha**2*beta**2 + 6*alpha**2*beta + 6*alpha**2 + 2*alpha*beta**2 + 6*alpha*beta + 4*alpha + beta**2 + 2*beta + 
1)*(x**5) + 5*(6*alpha**2 + 6*alpha*beta - 4*alpha + beta**2 - 2*beta)/(alpha**4 + 2*alpha**3*beta + 4*alpha**3 + alpha**2*beta**2 + 6*alpha**2*beta + 6*alpha**2 + 2*alpha*beta**2 + 
6*alpha*beta + 4*alpha + beta**2 + 2*beta + 1)*(x**4) + 4*(4*alpha**3 + 6*alpha**2*beta - 6*alpha**2 + 2*alpha*beta**2 - 6*alpha*beta - beta**2)/(alpha**4 + 2*alpha**3*beta + 4*alpha**3 + alpha**2*beta**2 + 6*alpha**2*beta + 6*alpha**2 + 2*alpha*beta**2 + 6*alpha*beta + 4*alpha + beta**2 + 2*beta + 1)*(x**3) + 3*(alpha**4 + 2*alpha**3*beta - 4*alpha**3 + alpha**2*beta**2 - 6*alpha**2*beta - 2*alpha*beta**2)/(alpha**4 + 2*alpha**3*beta + 4*alpha**3 + alpha**2*beta**2 + 6*alpha**2*beta + 6*alpha**2 + 2*alpha*beta**2 + 6*alpha*beta + 4*alpha + beta**2 + 2*beta + 1)*(x**2) + 2*(-alpha**4 - 2*alpha**3*beta - alpha**2*beta**2)/(alpha**4 + 2*alpha**3*beta + 4*alpha**3 + alpha**2*beta**2 + 6*alpha**2*beta + 6*alpha**2 + 2*alpha*beta**2 + 6*alpha*beta + 4*alpha + beta**2 + 2*beta + 1)*x + 0 )

#########################################################################################################
def d0_horner(x, alpha, beta) -> float:
    return ( x**2*(alpha**2*(alpha*(2*alpha + 9*beta + 2) + beta*(7*beta + 5))/(alpha**5*beta**3 + 5*alpha**4*beta**4 + 3*alpha**4*beta**3 + 10*alpha**3*beta**5 + 12*alpha**3*beta**4 + 3*alpha**3*beta**3 + 10*alpha**2*beta**6 + 18*alpha**2*beta**5 + 9*alpha**2*beta**4 + alpha**2*beta**3 + 5*alpha*beta**7 + 12*alpha*beta**6 + 9*alpha*beta**5 + 2*alpha*beta**4 + beta**8 + 3*beta**7 + 3*beta**6 + beta**5) + x*(alpha*(alpha*(alpha*(alpha*(-4*alpha - 22*beta + 2) + beta*(16 - 32*beta) + 6) + beta*(beta*(28 - 14*beta) + 18)) + beta**2*(14*beta + 10))/(alpha**6*beta**3 + 6*alpha**5*beta**4 + 3*alpha**5*beta**3 + 15*alpha**4*beta**5 + 15*alpha**4*beta**4 + 3*alpha**4*beta**3 + 20*alpha**3*beta**6 + 30*alpha**3*beta**5 + 12*alpha**3*beta**4 + alpha**3*beta**3 + 15*alpha**2*beta**7 + 30*alpha**2*beta**6 + 18*alpha**2*beta**5 + 3*alpha**2*beta**4 + 6*alpha*beta**8 + 15*alpha*beta**7 + 12*alpha*beta**6 + 3*alpha*beta**5 + beta**9 + 3*beta**8 + 3*beta**7 + beta**6) + x*(x*(x*(x*(alpha*(2*alpha + 8*beta + 2) + beta*(6*beta + 4))/(alpha**6*beta**3 + 6*alpha**5*beta**4 + 3*alpha**5*beta**3 + 15*alpha**4*beta**5 + 15*alpha**4*beta**4 + 3*alpha**4*beta**3 + 20*alpha**3*beta**6 + 30*alpha**3*beta**5 + 12*alpha**3*beta**4 + alpha**3*beta**3 + 15*alpha**2*beta**7 + 30*alpha**2*beta**6 + 18*alpha**2*beta**5 + 3*alpha**2*beta**4 + 
6*alpha*beta**8 + 15*alpha*beta**7 + 12*alpha*beta**6 + 3*alpha*beta**5 + beta**9 + 3*beta**8 + 3*beta**7 + beta**6) + (alpha*(alpha*(6*alpha + 27*beta + 2) + beta*(28*beta - 1) - 4) + beta*(beta*(7*beta - 7) - 8))/(alpha**6*beta**3 + 6*alpha**5*beta**4 + 3*alpha**5*beta**3 + 15*alpha**4*beta**5 + 15*alpha**4*beta**4 + 3*alpha**4*beta**3 + 20*alpha**3*beta**6 + 30*alpha**3*beta**5 + 12*alpha**3*beta**4 + alpha**3*beta**3 + 15*alpha**2*beta**7 + 30*alpha**2*beta**6 + 18*alpha**2*beta**5 + 3*alpha**2*beta**4 + 6*alpha*beta**8 + 
15*alpha*beta**7 + 12*alpha*beta**6 + 3*alpha*beta**5 + beta**9 + 3*beta**8 + 3*beta**7 + beta**6)) + (alpha*(alpha*(alpha*(6*alpha + 30*beta - 6) + beta*(38*beta - 36) - 10) + beta*(beta*(14*beta - 46) - 22) + 2) + beta*(beta*(-14*beta - 4) + 4))/(alpha**6*beta**3 
+ 6*alpha**5*beta**4 + 3*alpha**5*beta**3 + 15*alpha**4*beta**5 + 15*alpha**4*beta**4 + 3*alpha**4*beta**3 + 20*alpha**3*beta**6 + 30*alpha**3*beta**5 + 12*alpha**3*beta**4 + alpha**3*beta**3 + 15*alpha**2*beta**7 + 30*alpha**2*beta**6 + 18*alpha**2*beta**5 + 3*alpha**2*beta**4 + 6*alpha*beta**8 + 15*alpha*beta**7 + 12*alpha*beta**6 + 3*alpha*beta**5 + beta**9 + 3*beta**8 + 3*beta**7 + beta**6)) 
+ (alpha*(alpha*(alpha*(alpha*(2*alpha + 11*beta - 10) + beta*(16*beta - 53) - 6) + beta*(beta*(7*beta - 71) - 9) + 6) + beta*(beta*(8 - 28*beta) + 15)) + beta**2*(7*beta + 5))/(alpha**6*beta**3 + 6*alpha**5*beta**4 + 3*alpha**5*beta**3 + 15*alpha**4*beta**5 + 15*alpha**4*beta**4 + 3*alpha**4*beta**3 + 20*alpha**3*beta**6 + 30*alpha**3*beta**5 + 12*alpha**3*beta**4 + alpha**3*beta**3 + 15*alpha**2*beta**7 + 30*alpha**2*beta**6 + 18*alpha**2*beta**5 + 3*alpha**2*beta**4 + 6*alpha*beta**8 + 15*alpha*beta**7 + 12*alpha*beta**6 + 
3*alpha*beta**5 + beta**9 + 3*beta**8 + 3*beta**7 + beta**6)))) )

def d1_horner(x, alpha, beta) -> float:
    return ( x**2*(alpha**2/(alpha**3*beta**2 + 3*alpha**2*beta**3 + 2*alpha**2*beta**2 + 3*alpha*beta**4 + 4*alpha*beta**3 + alpha*beta**2 + beta**5 + 2*beta**4 + beta**3) + x*(alpha*(alpha*(-2*alpha - 2*beta + 3) + 2*beta)/(alpha**4*beta**2 + 4*alpha**3*beta**3 + 
2*alpha**3*beta**2 + 6*alpha**2*beta**4 + 6*alpha**2*beta**3 + alpha**2*beta**2 + 4*alpha*beta**5 + 6*alpha*beta**4 + 2*alpha*beta**3 + beta**6 + 2*beta**5 + beta**4) + x*(x*(x*(x/(alpha**4*beta**2 + 4*alpha**3*beta**3 + 2*alpha**3*beta**2 + 6*alpha**2*beta**4 + 6*alpha**2*beta**3 + alpha**2*beta**2 + 4*alpha*beta**5 + 6*alpha*beta**4 + 2*alpha*beta**3 + beta**6 + 2*beta**5 + beta**4) + (3*alpha 
+ beta - 2)/(alpha**4*beta**2 + 4*alpha**3*beta**3 + 2*alpha**3*beta**2 + 6*alpha**2*beta**4 + 6*alpha**2*beta**3 + alpha**2*beta**2 
+ 4*alpha*beta**5 + 6*alpha*beta**4 + 2*alpha*beta**3 + beta**6 + 2*beta**5 + beta**4)) + (alpha*(3*alpha + 2*beta - 6) - 2*beta + 1)/(alpha**4*beta**2 + 4*alpha**3*beta**3 + 2*alpha**3*beta**2 + 6*alpha**2*beta**4 + 6*alpha**2*beta**3 + alpha**2*beta**2 + 4*alpha*beta**5 + 6*alpha*beta**4 + 2*alpha*beta**3 + beta**6 + 2*beta**5 + beta**4)) + (alpha*(alpha*(alpha + beta - 6) - 4*beta + 3) + beta)/(alpha**4*beta**2 + 4*alpha**3*beta**3 + 2*alpha**3*beta**2 + 6*alpha**2*beta**4 + 6*alpha**2*beta**3 + alpha**2*beta**2 + 4*alpha*beta**5 + 6*alpha*beta**4 + 2*alpha*beta**3 + beta**6 + 2*beta**5 + beta**4)))) )

def d2_horner(x, alpha, beta) -> float:
    return ( x**2*(x*(x*(x*(x*(x*(alpha*(-2*alpha + 4*beta - 2) + 2*beta)/(alpha**6*beta**3 + 3*alpha**5*beta**3 + 3*alpha**4*beta**3 + alpha**3*beta**3) + (alpha*(alpha*(-6*alpha + 9*beta - 2) + beta*(8*beta - 5) + 4) + beta*(4*beta - 4))/(alpha**6*beta**3 + 3*alpha**5*beta**3 + 3*alpha**4*beta**3 + alpha**3*beta**3)) + (alpha*(alpha*(alpha*(-6*alpha + 6*beta + 6) + beta*(16*beta - 18) + 10) + beta*(beta*(4*beta - 8) - 2) - 2) + beta*(beta*(2*beta - 8) + 2))/(alpha**6*beta**3 + 3*alpha**5*beta**3 + 3*alpha**4*beta**3 + alpha**3*beta**3)) + (alpha*(alpha*(alpha*(alpha*(-2*alpha + beta + 10) + beta*(8*beta - 13) + 6) + beta*(beta*(5*beta - 28) + 9) - 6) + beta*(beta*(-5*beta - 8) + 3)) + beta**2*(4 - 4*beta))/(alpha**6*beta**3 + 3*alpha**5*beta**3 + 3*alpha**4*beta**3 + alpha**3*beta**3)) + (alpha*(alpha*(alpha*(alpha*(4*alpha - 2*beta - 2) + beta*(8 - 16*beta) - 6) + beta**2*(8 - 10*beta)) + beta**2*(8 - 2*beta)) + 2*beta**3)/(alpha**6*beta**3 + 3*alpha**5*beta**3 + 3*alpha**4*beta**3 + alpha**3*beta**3)) + (alpha*(alpha*(alpha*(-2*alpha + beta - 2) + beta*(8*beta - 1)) + beta**2*(5*beta + 4)) + 3*beta**3)/(alpha**5*beta**3 + 3*alpha**4*beta**3 + 3*alpha**3*beta**3 + alpha**2*beta**3)) )

def d3_horner(x, alpha, beta) -> float:
    return ( x**2*(x*(x*(x*(x*(x/(alpha**4*beta**2 + 2*alpha**3*beta**2 + alpha**2*beta**2) + (3*alpha + 2*beta - 2)/(alpha**4*beta**2 + 2*alpha**3*beta**2 + alpha**2*beta**2)) + (alpha*(3*alpha + 4*beta - 6) + beta*(beta - 4) + 1)/(alpha**4*beta**2 + 2*alpha**3*beta**2 + alpha**2*beta**2)) + (alpha*(alpha*(alpha + 2*beta - 6) + beta*(beta - 8) + 3) + beta*(2 - 2*beta))/(alpha**4*beta**2 + 2*alpha**3*beta**2 + alpha**2*beta**2)) + (alpha*(alpha*(-2*alpha - 4*beta + 3) + beta*(4 - 2*beta)) + beta**2)/(alpha**4*beta**2 + 2*alpha**3*beta**2 + alpha**2*beta**2)) + (alpha*(alpha + 2*beta) + beta**2)/(alpha**3*beta**2 + 2*alpha**2*beta**2 + alpha*beta**2)) )     

def d4_horner(x, alpha, beta) -> float:
    return ( x**2*(x*(x*(x*(x*(x*(alpha*(2*alpha + 2*beta - 4) - 2*beta)/(alpha**6 + 3*alpha**5*beta + 3*alpha**4*beta**2 + alpha**3*beta**3) + (alpha*(alpha*(8*alpha + 12*beta - 19) + beta*(4*beta - 19) + 8) + beta*(4 - 4*beta))/(alpha**6 + 3*alpha**5*beta + 3*alpha**4*beta**2 + alpha**3*beta**3)) + (alpha*(alpha*(alpha*(12*alpha + 24*beta - 36) + beta*(14*beta - 54) + 32) + beta*(beta*(2*beta - 22) + 32) - 4) + beta*(beta*(8 - 2*beta) - 2))/(alpha**6 + 3*alpha**5*beta + 3*alpha**4*beta**2 + alpha**3*beta**3)) + (alpha*(alpha*(alpha*(alpha*(8*alpha + 20*beta - 34) + beta*(16*beta - 68) + 48) + beta*(beta*(4*beta - 41) + 72) - 15) + beta*(beta*(32 - 7*beta) - 15)) + beta**2*(4*beta - 4))/(alpha**6 + 3*alpha**5*beta + 3*alpha**4*beta**2 + alpha**3*beta**3)) + (alpha*(alpha*(alpha*(alpha*(alpha*(2*alpha + 6*beta - 16) + beta*(6*beta - 40) + 32) + beta*(beta*(2*beta - 32) + 64) - 20) + beta*(beta*(40 - 8*beta) - 30)) + beta**2*(8*beta - 14)) - 2*beta**3)/(alpha**6 + 3*alpha**5*beta + 3*alpha**4*beta**2 + alpha**3*beta**3)) + (alpha*(alpha*(alpha*(-3*alpha - 6*beta + 8) + beta*(12 - 3*beta) - 10) + beta*(4*beta - 10)) - 3*beta**2)/(alpha**4 + 2*alpha**3*beta + alpha**2*beta**2)) + 1 )

def d5_horner(x, alpha, beta) -> float:
    return ( x*(x*(x*(x*(x*(x*(x/(alpha**4 + 2*alpha**3*beta + alpha**2*beta**2) + (4*alpha + 2*beta - 2)/(alpha**4 + 2*alpha**3*beta + alpha**2*beta**2)) + (alpha*(6*alpha + 6*beta - 8) + beta*(beta - 4) + 1)/(alpha**4 + 2*alpha**3*beta + alpha**2*beta**2)) + (alpha*(alpha*(4*alpha + 6*beta - 12) + beta*(2*beta - 12) + 4) + beta*(2 - 2*beta))/(alpha**4 + 2*alpha**3*beta + alpha**2*beta**2)) + (alpha*(alpha*(alpha*(alpha + 2*beta - 8) + beta*(beta - 12) + 6) + beta*(6 - 4*beta)) + beta**2)/(alpha**4 + 2*alpha**3*beta + alpha**2*beta**2)) + (alpha*(-2*alpha - 2*beta + 4) + 2*beta)/(alpha**2 + alpha*beta)) + 1) )

def d6_horner(x, alpha, beta) -> float:
    return ( x**2*(alpha**2*(alpha*(alpha*(alpha*(3*alpha + 9*beta + 10) + beta*(9*beta + 25) + 7) + beta*(beta*(3*beta + 20) + 14)) 
+ beta**2*(5*beta + 7))/(alpha**6 + 3*alpha**5*beta + 6*alpha**5 + 3*alpha**4*beta**2 + 15*alpha**4*beta + 15*alpha**4 + alpha**3*beta**3 + 12*alpha**3*beta**2 + 30*alpha**3*beta + 20*alpha**3 + 3*alpha**2*beta**3 + 18*alpha**2*beta**2 + 30*alpha**2*beta + 15*alpha**2 + 3*alpha*beta**3 + 12*alpha*beta**2 + 15*alpha*beta + 6*alpha + beta**3 + 3*beta**2 + 3*beta + 1) + x*(alpha*(alpha*(alpha*(alpha*(alpha*(-2*alpha - 6*beta + 4) + beta*(10 - 6*beta) + 34) + beta*(beta*(8 - 2*beta) + 68) + 28) + beta*(beta*(2*beta + 44) + 42)) + 
beta**2*(10*beta + 14))/(alpha**6 + 3*alpha**5*beta + 6*alpha**5 + 3*alpha**4*beta**2 + 15*alpha**4*beta + 15*alpha**4 + alpha**3*beta**3 + 12*alpha**3*beta**2 + 30*alpha**3*beta + 20*alpha**3 + 3*alpha**2*beta**3 + 18*alpha**2*beta**2 + 30*alpha**2*beta + 15*alpha**2 + 3*alpha*beta**3 + 12*alpha*beta**2 + 15*alpha*beta + 6*alpha + beta**3 + 3*beta**2 + 3*beta + 1) + x*(x*(x*(x*(alpha*(-2*alpha - 2*beta - 8) - 4*beta - 6)/(alpha**6 + 3*alpha**5*beta + 6*alpha**5 + 3*alpha**4*beta**2 + 15*alpha**4*beta + 15*alpha**4 + alpha**3*beta**3 + 12*alpha**3*beta**2 + 30*alpha**3*beta + 20*alpha**3 + 3*alpha**2*beta**3 + 18*alpha**2*beta**2 + 30*alpha**2*beta + 15*alpha**2 + 3*alpha*beta**3 + 12*alpha*beta**2 + 15*alpha*beta + 6*alpha + beta**3 + 3*beta**2 + 3*beta + 1) + (alpha*(alpha*(-8*alpha - 
12*beta - 29) + beta*(-4*beta - 29) - 14) + beta*(-8*beta - 7) + 7)/(alpha**6 + 3*alpha**5*beta + 6*alpha**5 + 3*alpha**4*beta**2 + 15*alpha**4*beta + 15*alpha**4 + alpha**3*beta**3 + 12*alpha**3*beta**2 + 30*alpha**3*beta + 20*alpha**3 + 3*alpha**2*beta**3 + 18*alpha**2*beta**2 + 30*alpha**2*beta + 15*alpha**2 + 3*alpha*beta**3 + 12*alpha*beta**2 + 15*alpha*beta + 6*alpha + beta**3 + 3*beta**2 + 3*beta + 1)) + (alpha*(alpha*(alpha*(-12*alpha - 24*beta - 36) + beta*(-14*beta - 54) + 4) + beta*(beta*(-2*beta - 26) + 4) + 28) + 
beta*(beta*(4 - 4*beta) + 14))/(alpha**6 + 3*alpha**5*beta + 6*alpha**5 + 3*alpha**4*beta**2 + 15*alpha**4*beta + 15*alpha**4 + alpha**3*beta**3 + 12*alpha**3*beta**2 + 30*alpha**3*beta + 20*alpha**3 + 3*alpha**2*beta**3 + 18*alpha**2*beta**2 + 30*alpha**2*beta + 15*alpha**2 + 3*alpha*beta**3 + 12*alpha*beta**2 + 15*alpha*beta + 6*alpha + beta**3 + 3*beta**2 + 3*beta + 1)) + (alpha*(alpha*(alpha*(alpha*(-8*alpha - 20*beta - 14) + beta*(-16*beta - 28) + 36) + beta*(beta*(-4*beta - 19) + 54) + 42) + beta*(beta*(28 - 5*beta) + 42)) + beta**2*(5*beta + 7))/(alpha**6 + 3*alpha**5*beta + 6*alpha**5 + 3*alpha**4*beta**2 + 15*alpha**4*beta + 15*alpha**4 + alpha**3*beta**3 + 12*alpha**3*beta**2 + 30*alpha**3*beta + 20*alpha**3 + 3*alpha**2*beta**3 + 18*alpha**2*beta**2 + 30*alpha**2*beta + 15*alpha**2 + 3*alpha*beta**3 + 12*alpha*beta**2 + 15*alpha*beta + 6*alpha + beta**3 + 3*beta**2 + 3*beta + 1)))) )

def d7_horner(x, alpha, beta) -> float:
    return ( x**2*(alpha**2*(alpha*(-alpha - 2*beta) - beta**2)/(alpha**4 + 2*alpha**3*beta + 4*alpha**3 + alpha**2*beta**2 + 6*alpha**2*beta + 6*alpha**2 + 2*alpha*beta**2 + 6*alpha*beta + 4*alpha + beta**2 + 2*beta + 1) + x*(alpha*(alpha*(alpha*(alpha + 2*beta - 4) + beta*(beta - 6)) - 2*beta**2)/(alpha**4 + 2*alpha**3*beta + 4*alpha**3 + alpha**2*beta**2 + 6*alpha**2*beta + 6*alpha**2 + 2*alpha*beta**2 + 6*alpha*beta + 4*alpha + beta**2 + 2*beta + 1) + x*(x*(x*(x/(alpha**4 + 2*alpha**3*beta + 4*alpha**3 + alpha**2*beta**2 + 6*alpha**2*beta + 6*alpha**2 + 2*alpha*beta**2 + 6*alpha*beta + 4*alpha + beta**2 + 2*beta + 1) + (4*alpha + 2*beta - 1)/(alpha**4 + 2*alpha**3*beta + 4*alpha**3 + alpha**2*beta**2 + 6*alpha**2*beta + 6*alpha**2 + 2*alpha*beta**2 + 6*alpha*beta + 4*alpha + beta**2 
+ 2*beta + 1)) + (alpha*(6*alpha + 6*beta - 4) + beta*(beta - 2))/(alpha**4 + 2*alpha**3*beta + 4*alpha**3 + alpha**2*beta**2 + 6*alpha**2*beta + 6*alpha**2 + 2*alpha*beta**2 + 6*alpha*beta + 4*alpha + beta**2 + 2*beta + 1)) + (alpha*(alpha*(4*alpha + 6*beta - 6) + 
beta*(2*beta - 6)) - beta**2)/(alpha**4 + 2*alpha**3*beta + 4*alpha**3 + alpha**2*beta**2 + 6*alpha**2*beta + 6*alpha**2 + 2*alpha*beta**2 + 6*alpha*beta + 4*alpha + beta**2 + 2*beta + 1)))) )

def d0_prime_horner(x, alpha, beta) -> float:
    return ( x*(alpha**2*(alpha*(4*alpha + 18*beta + 4) + beta*(14*beta + 10))/(alpha**5*beta**3 + 5*alpha**4*beta**4 + 3*alpha**4*beta**3 + 10*alpha**3*beta**5 + 12*alpha**3*beta**4 + 3*alpha**3*beta**3 + 10*alpha**2*beta**6 + 18*alpha**2*beta**5 + 9*alpha**2*beta**4 + alpha**2*beta**3 + 5*alpha*beta**7 + 12*alpha*beta**6 + 9*alpha*beta**5 + 2*alpha*beta**4 + beta**8 + 3*beta**7 + 3*beta**6 + beta**5) + x*(alpha*(alpha*(alpha*(alpha*(-12*alpha - 66*beta + 6) + beta*(48 - 96*beta) + 18) + beta*(beta*(84 - 42*beta) + 54)) + beta**2*(42*beta + 30))/(alpha**6*beta**3 + 6*alpha**5*beta**4 + 3*alpha**5*beta**3 + 15*alpha**4*beta**5 + 15*alpha**4*beta**4 + 3*alpha**4*beta**3 + 20*alpha**3*beta**6 + 30*alpha**3*beta**5 + 12*alpha**3*beta**4 + alpha**3*beta**3 + 15*alpha**2*beta**7 + 30*alpha**2*beta**6 + 18*alpha**2*beta**5 + 3*alpha**2*beta**4 + 6*alpha*beta**8 + 15*alpha*beta**7 + 12*alpha*beta**6 + 3*alpha*beta**5 + beta**9 + 3*beta**8 + 3*beta**7 + beta**6) + x*(x*(x*(x*(alpha*(14*alpha + 56*beta + 14) + beta*(42*beta + 28))/(alpha**6*beta**3 + 6*alpha**5*beta**4 + 3*alpha**5*beta**3 + 15*alpha**4*beta**5 + 15*alpha**4*beta**4 + 3*alpha**4*beta**3 + 20*alpha**3*beta**6 + 30*alpha**3*beta**5 + 12*alpha**3*beta**4 + alpha**3*beta**3 + 15*alpha**2*beta**7 + 30*alpha**2*beta**6 + 18*alpha**2*beta**5 + 3*alpha**2*beta**4 + 6*alpha*beta**8 + 15*alpha*beta**7 + 12*alpha*beta**6 + 3*alpha*beta**5 + beta**9 + 3*beta**8 + 3*beta**7 + beta**6) + (alpha*(alpha*(36*alpha + 162*beta + 12) + beta*(168*beta - 6) - 24) + beta*(beta*(42*beta - 42) - 48))/(alpha**6*beta**3 + 6*alpha**5*beta**4 + 3*alpha**5*beta**3 + 15*alpha**4*beta**5 + 15*alpha**4*beta**4 + 3*alpha**4*beta**3 + 20*alpha**3*beta**6 + 30*alpha**3*beta**5 
+ 12*alpha**3*beta**4 + alpha**3*beta**3 + 15*alpha**2*beta**7 + 30*alpha**2*beta**6 + 18*alpha**2*beta**5 + 3*alpha**2*beta**4 + 6*alpha*beta**8 + 15*alpha*beta**7 + 12*alpha*beta**6 + 3*alpha*beta**5 + beta**9 + 3*beta**8 + 3*beta**7 + beta**6)) + (alpha*(alpha*(alpha*(30*alpha + 150*beta - 30) + beta*(190*beta - 180) - 50) + beta*(beta*(70*beta - 230) - 110) + 10) + beta*(beta*(-70*beta - 20) 
+ 20))/(alpha**6*beta**3 + 6*alpha**5*beta**4 + 3*alpha**5*beta**3 + 15*alpha**4*beta**5 + 15*alpha**4*beta**4 + 3*alpha**4*beta**3 + 20*alpha**3*beta**6 + 30*alpha**3*beta**5 + 12*alpha**3*beta**4 + alpha**3*beta**3 + 15*alpha**2*beta**7 + 30*alpha**2*beta**6 + 18*alpha**2*beta**5 + 3*alpha**2*beta**4 + 6*alpha*beta**8 + 15*alpha*beta**7 + 12*alpha*beta**6 + 3*alpha*beta**5 + beta**9 + 3*beta**8 + 3*beta**7 + beta**6)) + (alpha*(alpha*(alpha*(alpha*(8*alpha + 44*beta - 40) + beta*(64*beta - 212) - 24) + beta*(beta*(28*beta - 
284) - 36) + 24) + beta*(beta*(32 - 112*beta) + 60)) + beta**2*(28*beta + 20))/(alpha**6*beta**3 + 6*alpha**5*beta**4 + 3*alpha**5*beta**3 + 15*alpha**4*beta**5 + 15*alpha**4*beta**4 + 3*alpha**4*beta**3 + 20*alpha**3*beta**6 + 30*alpha**3*beta**5 + 12*alpha**3*beta**4 + alpha**3*beta**3 + 15*alpha**2*beta**7 + 30*alpha**2*beta**6 + 18*alpha**2*beta**5 + 3*alpha**2*beta**4 + 6*alpha*beta**8 + 15*alpha*beta**7 + 12*alpha*beta**6 + 3*alpha*beta**5 + beta**9 + 3*beta**8 + 3*beta**7 + beta**6)))) )

def d1_prime_horner(x, alpha, beta) -> float:
    return ( x*(2*alpha**2/(alpha**3*beta**2 + 3*alpha**2*beta**3 + 2*alpha**2*beta**2 + 3*alpha*beta**4 + 4*alpha*beta**3 + alpha*beta**2 + beta**5 + 2*beta**4 + beta**3) + x*(alpha*(alpha*(-6*alpha - 6*beta + 9) + 6*beta)/(alpha**4*beta**2 + 4*alpha**3*beta**3 + 2*alpha**3*beta**2 + 6*alpha**2*beta**4 + 6*alpha**2*beta**3 + alpha**2*beta**2 + 4*alpha*beta**5 + 6*alpha*beta**4 + 2*alpha*beta**3 
+ beta**6 + 2*beta**5 + beta**4) + x*(x*(x*(7*x/(alpha**4*beta**2 + 4*alpha**3*beta**3 + 2*alpha**3*beta**2 + 6*alpha**2*beta**4 + 6*alpha**2*beta**3 + alpha**2*beta**2 + 4*alpha*beta**5 + 6*alpha*beta**4 + 2*alpha*beta**3 + beta**6 + 2*beta**5 + beta**4) + (18*alpha + 6*beta - 12)/(alpha**4*beta**2 + 4*alpha**3*beta**3 + 2*alpha**3*beta**2 + 6*alpha**2*beta**4 + 6*alpha**2*beta**3 + alpha**2*beta**2 + 4*alpha*beta**5 + 6*alpha*beta**4 + 2*alpha*beta**3 + beta**6 + 2*beta**5 + beta**4)) + (alpha*(15*alpha + 10*beta - 30) - 10*beta + 5)/(alpha**4*beta**2 + 4*alpha**3*beta**3 + 2*alpha**3*beta**2 + 6*alpha**2*beta**4 + 6*alpha**2*beta**3 + alpha**2*beta**2 + 
4*alpha*beta**5 + 6*alpha*beta**4 + 2*alpha*beta**3 + beta**6 + 2*beta**5 + beta**4)) + (alpha*(alpha*(4*alpha + 4*beta - 24) - 16*beta + 12) + 4*beta)/(alpha**4*beta**2 + 4*alpha**3*beta**3 + 2*alpha**3*beta**2 + 6*alpha**2*beta**4 + 6*alpha**2*beta**3 + alpha**2*beta**2 + 4*alpha*beta**5 + 6*alpha*beta**4 + 2*alpha*beta**3 + beta**6 + 2*beta**5 + beta**4)))) )

def d2_prime_horner(x, alpha, beta) -> float:
    return ( x*(x*(x*(x*(x*(x*(alpha*(-14*alpha + 28*beta - 14) + 14*beta)/(alpha**6*beta**3 + 3*alpha**5*beta**3 + 3*alpha**4*beta**3 + alpha**3*beta**3) + (alpha*(alpha*(-36*alpha + 54*beta - 12) + beta*(48*beta - 30) + 24) + beta*(24*beta - 24))/(alpha**6*beta**3 + 3*alpha**5*beta**3 + 3*alpha**4*beta**3 + alpha**3*beta**3)) + (alpha*(alpha*(alpha*(-30*alpha + 30*beta + 30) + beta*(80*beta - 90) + 50) + beta*(beta*(20*beta - 40) - 10) - 10) + beta*(beta*(10*beta - 40) + 10))/(alpha**6*beta**3 + 3*alpha**5*beta**3 + 3*alpha**4*beta**3 + alpha**3*beta**3)) + (alpha*(alpha*(alpha*(alpha*(-8*alpha + 4*beta + 40) + beta*(32*beta - 52) + 24) + beta*(beta*(20*beta - 112) + 36) - 24) + beta*(beta*(-20*beta - 32) + 12)) + beta**2*(16 - 16*beta))/(alpha**6*beta**3 + 3*alpha**5*beta**3 + 3*alpha**4*beta**3 + alpha**3*beta**3)) + (alpha*(alpha*(alpha*(alpha*(12*alpha - 6*beta - 6) + beta*(24 - 48*beta) - 18) + beta**2*(24 - 30*beta)) + beta**2*(24 - 6*beta)) + 6*beta**3)/(alpha**6*beta**3 + 3*alpha**5*beta**3 + 3*alpha**4*beta**3 + alpha**3*beta**3)) + (alpha*(alpha*(alpha*(-4*alpha + 2*beta - 4) + beta*(16*beta - 2)) + beta**2*(10*beta + 8)) + 6*beta**3)/(alpha**5*beta**3 + 3*alpha**4*beta**3 + 3*alpha**3*beta**3 + alpha**2*beta**3)) )

def d3_prime_horner(x, alpha, beta) -> float:
    return ( x*(x*(x*(x*(x*(7*x/(alpha**4*beta**2 + 2*alpha**3*beta**2 + alpha**2*beta**2) + (18*alpha + 12*beta - 12)/(alpha**4*beta**2 + 2*alpha**3*beta**2 + alpha**2*beta**2)) + (alpha*(15*alpha + 20*beta - 30) + beta*(5*beta - 20) + 5)/(alpha**4*beta**2 + 2*alpha**3*beta**2 + alpha**2*beta**2)) + (alpha*(alpha*(4*alpha + 8*beta - 24) + beta*(4*beta - 32) + 12) + beta*(8 - 8*beta))/(alpha**4*beta**2 + 2*alpha**3*beta**2 + alpha**2*beta**2)) + (alpha*(alpha*(-6*alpha - 12*beta + 9) + beta*(12 - 6*beta)) + 3*beta**2)/(alpha**4*beta**2 + 2*alpha**3*beta**2 + alpha**2*beta**2)) + (alpha*(2*alpha + 4*beta) + 2*beta**2)/(alpha**3*beta**2 + 2*alpha**2*beta**2 + alpha*beta**2)) )

def d4_prime_horner(x, alpha, beta) -> float:
    return ( x*(x*(x*(x*(x*(x*(alpha*(14*alpha + 14*beta - 28) - 14*beta)/(alpha**6 + 3*alpha**5*beta + 3*alpha**4*beta**2 + alpha**3*beta**3) + (alpha*(alpha*(48*alpha + 72*beta - 114) + beta*(24*beta - 114) + 48) + beta*(24 - 24*beta))/(alpha**6 + 3*alpha**5*beta 
+ 3*alpha**4*beta**2 + alpha**3*beta**3)) + (alpha*(alpha*(alpha*(60*alpha + 120*beta - 180) + beta*(70*beta - 270) + 160) + beta*(beta*(10*beta - 110) + 160) - 20) + beta*(beta*(40 - 10*beta) - 10))/(alpha**6 + 3*alpha**5*beta + 3*alpha**4*beta**2 + alpha**3*beta**3)) + (alpha*(alpha*(alpha*(alpha*(32*alpha + 80*beta - 136) + beta*(64*beta - 272) + 192) + beta*(beta*(16*beta - 164) + 288) - 60) 
+ beta*(beta*(128 - 28*beta) - 60)) + beta**2*(16*beta - 16))/(alpha**6 + 3*alpha**5*beta + 3*alpha**4*beta**2 + alpha**3*beta**3)) + (alpha*(alpha*(alpha*(alpha*(alpha*(6*alpha + 18*beta - 48) + beta*(18*beta - 120) + 96) + beta*(beta*(6*beta - 96) + 192) - 60) + beta*(beta*(120 - 24*beta) - 90)) + beta**2*(24*beta - 42)) - 6*beta**3)/(alpha**6 + 3*alpha**5*beta + 3*alpha**4*beta**2 + alpha**3*beta**3)) + (alpha*(alpha*(alpha*(-6*alpha - 12*beta + 16) + beta*(24 - 6*beta) - 20) + beta*(8*beta - 20)) - 6*beta**2)/(alpha**4 + 2*alpha**3*beta + alpha**2*beta**2)) )

def d5_prime_horner(x, alpha, beta) -> float:
    return ( x*(x*(x*(x*(x*(7*x/(alpha**4 + 2*alpha**3*beta + alpha**2*beta**2) + (24*alpha + 12*beta - 12)/(alpha**4 + 2*alpha**3*beta + alpha**2*beta**2)) + (alpha*(30*alpha + 30*beta - 40) + beta*(5*beta - 20) + 5)/(alpha**4 + 2*alpha**3*beta + alpha**2*beta**2)) + (alpha*(alpha*(16*alpha + 24*beta - 48) + beta*(8*beta - 48) + 16) + beta*(8 - 8*beta))/(alpha**4 + 2*alpha**3*beta + alpha**2*beta**2)) + (alpha*(alpha*(alpha*(3*alpha + 6*beta - 24) + beta*(3*beta - 36) + 18) + beta*(18 - 12*beta)) + 3*beta**2)/(alpha**4 + 2*alpha**3*beta + alpha**2*beta**2)) + (alpha*(-4*alpha - 4*beta + 8) + 4*beta)/(alpha**2 + alpha*beta)) + 1 )

def d6_prime_horner(x, alpha, beta) -> float:
    return ( x*(alpha**2*(alpha*(alpha*(alpha*(6*alpha + 18*beta + 20) + beta*(18*beta + 50) + 14) + beta*(beta*(6*beta + 40) + 28)) 
+ beta**2*(10*beta + 14))/(alpha**6 + 3*alpha**5*beta + 6*alpha**5 + 3*alpha**4*beta**2 + 15*alpha**4*beta + 15*alpha**4 + alpha**3*beta**3 + 12*alpha**3*beta**2 + 30*alpha**3*beta + 20*alpha**3 + 3*alpha**2*beta**3 + 18*alpha**2*beta**2 + 30*alpha**2*beta + 15*alpha**2 + 3*alpha*beta**3 + 12*alpha*beta**2 + 15*alpha*beta + 6*alpha + beta**3 + 3*beta**2 + 3*beta + 1) + x*(alpha*(alpha*(alpha*(alpha*(alpha*(-6*alpha - 18*beta + 12) + beta*(30 - 18*beta) + 102) + beta*(beta*(24 - 6*beta) + 204) + 84) + beta*(beta*(6*beta + 132) 
+ 126)) + beta**2*(30*beta + 42))/(alpha**6 + 3*alpha**5*beta + 6*alpha**5 + 3*alpha**4*beta**2 + 15*alpha**4*beta + 15*alpha**4 + alpha**3*beta**3 + 12*alpha**3*beta**2 + 30*alpha**3*beta + 20*alpha**3 + 3*alpha**2*beta**3 + 18*alpha**2*beta**2 + 30*alpha**2*beta + 15*alpha**2 + 3*alpha*beta**3 + 12*alpha*beta**2 + 15*alpha*beta + 6*alpha + beta**3 + 3*beta**2 + 3*beta + 1) + x*(x*(x*(x*(alpha*(-14*alpha - 14*beta - 56) - 28*beta - 42)/(alpha**6 + 3*alpha**5*beta + 6*alpha**5 + 3*alpha**4*beta**2 + 15*alpha**4*beta + 15*alpha**4 + alpha**3*beta**3 + 12*alpha**3*beta**2 + 30*alpha**3*beta + 20*alpha**3 + 3*alpha**2*beta**3 + 18*alpha**2*beta**2 + 30*alpha**2*beta + 15*alpha**2 + 3*alpha*beta**3 + 12*alpha*beta**2 + 15*alpha*beta + 6*alpha + beta**3 + 3*beta**2 + 3*beta + 1) + (alpha*(alpha*(-48*alpha - 72*beta - 174) + beta*(-24*beta - 174) - 84) + beta*(-48*beta - 42) + 42)/(alpha**6 + 3*alpha**5*beta + 6*alpha**5 + 
3*alpha**4*beta**2 + 15*alpha**4*beta + 15*alpha**4 + alpha**3*beta**3 + 12*alpha**3*beta**2 + 30*alpha**3*beta + 20*alpha**3 + 3*alpha**2*beta**3 + 18*alpha**2*beta**2 + 30*alpha**2*beta + 15*alpha**2 + 3*alpha*beta**3 + 12*alpha*beta**2 + 15*alpha*beta + 6*alpha + beta**3 + 3*beta**2 + 3*beta + 1)) + (alpha*(alpha*(alpha*(-60*alpha - 120*beta - 180) + beta*(-70*beta - 270) + 20) + beta*(beta*(-10*beta - 130) + 20) + 140) + beta*(beta*(20 - 20*beta) + 70))/(alpha**6 + 3*alpha**5*beta + 6*alpha**5 + 3*alpha**4*beta**2 + 15*alpha**4*beta + 15*alpha**4 + alpha**3*beta**3 + 12*alpha**3*beta**2 + 30*alpha**3*beta + 20*alpha**3 + 3*alpha**2*beta**3 + 18*alpha**2*beta**2 + 30*alpha**2*beta + 15*alpha**2 + 3*alpha*beta**3 + 12*alpha*beta**2 + 15*alpha*beta + 6*alpha + beta**3 + 3*beta**2 + 3*beta + 1)) + (alpha*(alpha*(alpha*(alpha*(-32*alpha - 80*beta - 56) + beta*(-64*beta - 112) + 144) + beta*(beta*(-16*beta - 76) + 216) 
+ 168) + beta*(beta*(112 - 20*beta) + 168)) + beta**2*(20*beta + 28))/(alpha**6 + 3*alpha**5*beta + 6*alpha**5 + 3*alpha**4*beta**2 + 15*alpha**4*beta + 15*alpha**4 + alpha**3*beta**3 + 12*alpha**3*beta**2 + 30*alpha**3*beta + 20*alpha**3 + 3*alpha**2*beta**3 + 18*alpha**2*beta**2 + 30*alpha**2*beta + 15*alpha**2 + 3*alpha*beta**3 + 12*alpha*beta**2 + 15*alpha*beta + 6*alpha + beta**3 + 3*beta**2 + 3*beta + 1)))) )

def d7_prime_horner(x, alpha, beta) -> float:
    return ( x*(alpha**2*(alpha*(-2*alpha - 4*beta) - 2*beta**2)/(alpha**4 + 2*alpha**3*beta + 4*alpha**3 + alpha**2*beta**2 + 6*alpha**2*beta + 6*alpha**2 + 2*alpha*beta**2 + 6*alpha*beta + 4*alpha + beta**2 + 2*beta + 1) + x*(alpha*(alpha*(alpha*(3*alpha + 6*beta 
- 12) + beta*(3*beta - 18)) - 6*beta**2)/(alpha**4 + 2*alpha**3*beta + 4*alpha**3 + alpha**2*beta**2 + 6*alpha**2*beta + 6*alpha**2 + 2*alpha*beta**2 + 6*alpha*beta + 4*alpha + beta**2 + 2*beta + 1) + x*(x*(x*(7*x/(alpha**4 + 2*alpha**3*beta + 4*alpha**3 + alpha**2*beta**2 + 6*alpha**2*beta + 6*alpha**2 + 2*alpha*beta**2 + 6*alpha*beta + 4*alpha + beta**2 + 2*beta + 1) + (24*alpha + 12*beta - 6)/(alpha**4 + 2*alpha**3*beta + 4*alpha**3 + alpha**2*beta**2 + 6*alpha**2*beta + 6*alpha**2 + 2*alpha*beta**2 + 6*alpha*beta + 4*alpha + beta**2 + 2*beta + 1)) + (alpha*(30*alpha + 30*beta - 20) + beta*(5*beta - 10))/(alpha**4 + 2*alpha**3*beta + 4*alpha**3 + alpha**2*beta**2 + 6*alpha**2*beta + 6*alpha**2 + 2*alpha*beta**2 + 6*alpha*beta + 4*alpha + beta**2 + 2*beta + 1)) + (alpha*(alpha*(16*alpha + 24*beta - 24) + beta*(8*beta - 24)) - 4*beta**2)/(alpha**4 + 2*alpha**3*beta + 4*alpha**3 + alpha**2*beta**2 + 6*alpha**2*beta + 
6*alpha**2 + 2*alpha*beta**2 + 6*alpha*beta + 4*alpha + beta**2 + 2*beta + 1)))) )

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