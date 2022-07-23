
# each row is one of the bi(theta)
B_I5 = [
    [1, -5.308169607103577, 10.18168044895868, -7.520036991611715, 0.9340485368631161, 0.746867191577065], 
    [0, 0, 0, 0, 0, 0], 
    [0, 0, 0, 0, 0, 0], 
    [0, 6.272050253212501, -16.02618147467746, 12.844356324519618, -1.1487945044767591, -1.6831681430145498], 
    [0, 6.876491702846304, -24.635767260846333, 33.21078648379717, -17.49461528263644, 2.4640414758066496], 
    [0, -35.5444517105996, 165.7016170190242, -385.4635395491143, 442.43241370157017, -182.7206429912112], 
    [0, 1918.6548566980114, -9268.121508966042, 20858.33702877255, -22645.82767158481, 8960.474176055992], 
    [0, -1883.0698021327182, 9101.025187200634, -20473.188551959534, 22209.765551256532, -8782.1682509635], 
    [0, 0.11902479635123643, -0.12502696705039376, 1.7799569193949991, -4.660932123043763, 2.886977374347921], 
    [0, -8, 32, -40, 16, 0]
]        
n_stages_for_interp =  10
n_degree_polys = 6

def sigma_prod(arr1, arr2, start, end):
    res = 0
    for i in range(start, end):
        res += (  arr1[i] * arr2[i]  )
    return res

#################################################################################
class CRK5Interp:
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
        # each row of B_I6 is an interp and is of length 6
        b_i = B_I5[i]
        # b_i taken correctly print(f"b{i} = ", b_i)
        theta_powers = [theta**pow for pow in range(1, n_degree_polys + 1)] 
        # theta powers are correct
        # print("theta powers", theta_powers)
        # print("actual powers", theta, theta**2, theta**3, theta**4, theta**5, theta**6)
    
        res = sigma_prod(b_i, theta_powers, 0, n_degree_polys)
        # b_eval is as expected....
        # actual_eval = (   
        #           b_i[0] * theta 
        #         + b_i[1] * theta ** 2 
        #         + b_i[2] * theta ** 3
        #         + b_i[3] * theta ** 4
        #         + b_i[4] * theta ** 5
        #         + b_i[5] * theta ** 6  
        # )
        # print("actual eval", actual_eval)
        # print("res", res)
        return res

class CRK5ContinuousSolution:
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
