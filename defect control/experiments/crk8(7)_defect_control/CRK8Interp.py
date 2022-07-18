# each row is one b_i(theta) 
B_I8 = [
    [1.0, -10.039154650554519, 53.79210495862331, -165.0579057235472, 298.026456543461, -311.91254487079004, 174.60598526911716, -40.37066163211959], 
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 
    [0.0, 158.1976739121776, -1543.96141721949, 6241.39874782878, -13136.516156406109, 15106.948493169599, -8996.489626298231, 2170.776389952444], 
    [0.0, 110.78115200797782, -1081.1905145356177, 4370.666940459977, -9199.113723922197, 10578.949209629855, -6299.975594978841, 1520.1305005543413], 
    [0.0, -7011.442038211314, 68429.55220744078, -276623.5714822198, 582220.4545548494, -669551.5244611246, 398731.3087623333, -96210.47174510667], 
    [0.0, 11206.397569848148, -109371.04854950662, 442127.8393698155, -930563.7629864562, 1070145.1335855902, -637292.8058429047, 153773.3309185794], 
    [0.0, -14179.231640455684, 138385.00931963572, -559415.549024087, 1177423.7946992505, -1354033.3227908213, 806353.893882505, -194566.3328138133], 
    [0.0, 10247.761767921746, -100015.05326375231, 404306.62401434296, -850959.9711689702, 978601.0462088685, -582776.4729907749, 140619.0037156383], 
    [0.0, -105.49303976850968, 1029.5801395803103, -4162.034181876453, 8759.996193602336, -10073.965556886049, 5999.247741473951, -1447.5674285888924], 
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 
    [0.0, -14.863613373267432, 145.76359364894867, -587.6557063401914, 1227.3721512545558, -1394.4931057405536, 816.8562950730669, -192.97961452255882], 
    [0.0, 14.349685752905462, -150.29493444816657, 629.481242570029, -1352.5182073090607, 1575.8969337088804, -946.7876580472948, 229.87293777270722], 
    [0.0, -102.54524701110401, 1074.0326612646807, -4498.377917100411, 9665.320624003281, -11261.62224831288, 6765.902468760784, -1642.7103416043497], 
    [0.0, -38.13206313286474, 399.3854658292329, -1672.7487204919717, 3594.1072548585666, -4187.7015568029265, 2515.9412806490636, -610.8516609091005], 
    [0.0, -66.38279583069588, 595.8297683881103, -2188.7370600929717, 4213.839795282853, -4484.035731929197, 2500.6482514253466, -571.1622272434449], 
    [0.0, -90.4188757317306, 931.9503884048154, -3962.898377713156, 8733.31742002555, -10445.908189887661, 6426.218942917599, -1592.261308015418], 
    [0.0, -59.738843630388715, 544.8870146891725, -2090.4303749263127, 4194.418982707227, -4603.369436819628, 2619.2014135592976, -604.9687555793671], 
    [0.0, -59.20053764683937, 571.7660156218088, -2308.9495644453605, 4881.2341106861395, -5660.118807771202, 3408.7066890374217, -833.4379054819676]
]
n_stages_for_interp =  21
n_degree_polys = 8


def sigma_prod(arr1, arr2, start, end):
    res = 0
    for i in range(start, end):
        res += (  arr1[i] * arr2[i]  )
    return res

class CRK8Interp:
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
        b_i = B_I8[i]
        theta_powers = [theta**pow for pow in range(1, n_degree_polys + 1)] 
        return sigma_prod(b_i, theta_powers, 0, n_degree_polys)

    def prime(self, x):
        theta = (x - self.x_i) / (self.h_i)
        b_primes = [self._b_prime(i, theta) for i in range(n_stages_for_interp)]
        return sigma_prod(self.k, b_primes, 0, n_stages_for_interp)

    def _b_prime(self, i, theta):
        b_i = B_I8[i]
        b_i_prime_coeff = [coeff * bij for (coeff, bij) in zip(range(1, n_degree_polys + 1), b_i)]
        theta_powers = [theta**pow for pow in range(n_degree_polys)] # go from 0 to n_degree_polys - 1
        res = sigma_prod(b_i_prime_coeff, theta_powers, 0, n_degree_polys)
        return res

class CRK8ContinuousSolution:
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
