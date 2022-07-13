
# each row is one of the bi(theta)
B_I6 = [
    [1, -7.778593856495575, 27.0524385722671, -45.780190114576975, 36.72377741043638, -11.183042432947357], 
    [0, 0, 0, 0, 0, 0], 
    [0, 0, 0, 0, 0, 0], 
    [0, 16.632102138279762, -86.25583404770623, 171.7330546182696, -149.67744091315947, 47.826380659879696], 
    [0, 27.10835046149758, -140.58676162962996, 279.9044757968917, -243.95644583707966, 77.95131832728772], 
    [0, 283.70753264670356, -1471.3371557366656, 2929.3928569314394, -2553.17199842168, 815.8141610498723], 
    [0, -11365.512865164834, 58942.74718938947, -117353.43045697975, 102281.77209230464, -32682.059078573824], 
    [0, 11100.250191051131, -57567.067013355576, 114614.48808378985, -99894.591091309, 31919.283963225014], 
    [0, -3.0022825150732126, 14.946122435958785, -27.826954732510288, 21.824672217437076, -5.941557405812358], 
    [0, -19.610347376201034, 93.13370014508226, -165.3493635542416, 129.73901617804057, -37.91300539268019], 
    [0, -18.23029074639409, 96.74593449012313, -199.08634973839895, 180.85605899200485, -60.285352997334954], 
    [0, -13.563796638614157, 90.62137973668116, -204.04515601697273, 190.48135937835858, -63.493786459452856]
]
n_degree_polys = 6
n_stages = 12

def sigma_prod(arr1, arr2, start, end):
    res = 0
    for i in range(start, end):
        res += (  arr1[i] * arr2[i]  )
    return res

#################################################################################
class CRK6Interp:
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
        b_evals = [self._b_eval(i, theta) for i in range(n_stages)]
        return self.y_i + self.h_i * sigma_prod(self.k, b_evals, 0, n_stages)
    
    def _b_eval(self, i, theta):
        # each row of B_I6 is an interp and is of length 6
        b_i = B_I6[i]
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

class CRK6ContinuousSolution:
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
