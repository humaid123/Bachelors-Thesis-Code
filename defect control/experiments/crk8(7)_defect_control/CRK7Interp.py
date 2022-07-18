# each row is one b_i(theta) 
B_I7 = [
    [1.0, -7.238550783576433, 26.00913483254676, -50.23684777762567, 52.12072084601022, -27.06472451211777, 5.454547288952965], 
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 
    [0.0, 11.153308875889351, -91.76096563989617, 291.70742417220595, -430.4096692910863, 299.45311881989977, -79.78911199784015], 
    [0.0, 2.3487522980730935, -11.672489417201843, -3.3391390765059286, 94.88526224972061, -143.07112658301202, 61.09670974442174], 
    [0.0, -1027.3216753392408, 9198.714323607608, -33189.78048157364, 57750.083134888715, -47698.93315706262, 14951.543653440334], 
    [0.0, 1568.546608927282, -13995.388525416005, 50256.21246981024, -86974.512803622, 71494.79770959976, -22324.571394333743], 
    [0.0, -2000.882061921042, 17864.363803476917, -64205.190751556285, 111224.84899303781, -91509.33921021303, 28594.46085938938], 
    [0.0, 1496.6204006934463, -13397.55405171476, 48323.56021994375, -84051.4283423393, 69399.85821115709, -21748.118154466232], 
    [0.0, -16.413207755609335, 147.60970454070025, -535.7199637147321, 938.2862470778207, -779.4383096393493, 245.43939702786273], 
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 
    [0.0, -4.296724431782465, 38.64447461116781, -140.3503471762809, 246.39546696975026, -205.83416869641673, 65.44129872356201], 
    [0.0, -20.416280692948217, 153.52132325248365, -436.550261021122, 598.2146442626508, -398.7823950071291, 104.01296920606484], 
    [0.0, 16.530071842642716, -96.6861433615782, 268.9599342195317, -428.6819097889647, 354.5782311524334, -114.70018406406496], 
    [0.0, -18.630641713134295, 164.1994112280183, -579.2722562495404, 980.1982557088668, -786.2241790155139, 239.7294100413036]
]
n_stages_for_interp =  17
n_degree_polys = 7

def sigma_prod(arr1, arr2, start, end):
    res = 0
    for i in range(start, end):
        res += (  arr1[i] * arr2[i]  )
    return res

class CRK7Interp:
    def __init__(self, k, y_i, x_i, x_i_plus_1):
        self.k = k[:] # make a copy of the stages
        self.x_i = x_i
        self.x_i_plus_1 = x_i_plus_1
        self.h_i = self.x_i_plus_1 - self.x_i
        self.y_i = y_i
    """
    the interpolant is of the form
    u(t + theta * h) = b1(theta)f1 + b2(theta)f2 + ... + b12(theta)f12 
    """
    def eval(self, x):
        theta = (x - self.x_i) / (self.h_i)
        b_evals = [self._b_eval(i, theta) for i in range(n_stages_for_interp)]
        return self.y_i + self.h_i * sigma_prod(self.k, b_evals, 0, n_stages_for_interp)
    
    def _b_eval(self, i, theta):
        # each row of B_I8 is an interp and is of length 8
        b_i = B_I7[i]
        theta_powers = [theta**pow for pow in range(1, n_degree_polys + 1)] 
        return sigma_prod(b_i, theta_powers, 0, n_degree_polys)

    def prime(self, x):
        theta = (x - self.x_i) / (self.h_i)
        b_primes = [self._b_prime(i, theta) for i in range(n_stages_for_interp)]
        return sigma_prod(self.k, b_primes, 0, n_stages_for_interp)

    def _b_prime(self, i, theta):
        b_i = B_I7[i]
        b_i_prime_coeff = [coeff * bij for (coeff, bij) in zip(range(1, n_degree_polys + 1), b_i)]
        theta_powers = [theta**pow for pow in range(n_degree_polys)] # go from 0 to n_degree_polys - 1
        res = sigma_prod(b_i_prime_coeff, theta_powers, 0, n_degree_polys)
        return res

class CRK7ContinuousSolution:
    def __init__(self):
        self.interps = []

    def append(self, interp):
        self.interps.append(interp)

    def extend(self, interps):
        self.interps.extend(interps)
    
    def eval(self, x) -> float:
        for interp in self.interps:
            if (interp.x_i <= x <= interp.x_i_plus_1):
                return interp.eval(x)

        print(f"ERROR: {x} is outside of the solution range: {self.interps[0].x_i} <= x <= {self.interps[-1].x_i_plus_1}")
        return -1
    def prime(self, x) -> float:
        for interp in self.interps:
            if (interp.x_i <= x <= interp.x_i_plus_1):
                return interp.prime(x)

        print(f"ERROR: {x} is outside of the solution range: {self.interps[0].x_i} <= x <= {self.interps[-1].x_i_plus_1}")
        return -1
        
    def create_error_sampling(self):
        res = []
        for interp in self.interps:
            res.append(
                (interp.x_i, interp.x_i_plus_1, interp)
            )
        return res

    def create_defect_sampling(self):
        res = []
        for interp in self.interps:
            res.append( (interp.x_i, interp.x_i_plus_1, interp) )
        return res
