
# these formulae were derived using sympy: see get_hb_coefficients.py and test_hb_coefficients.py for the derivations and the tests
# we create a class for each of these so that the class abstracts a Hermite Birkhoff interpolant
# we then need to instantiate one object of this class for each interval and we get to do Hermite Birkhoff interpolation

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

# I use a class to represent the Hermite Birkhoff interpolant
# we will have an instance of this class on each step
class HB:
    def __init__(   
        self, 
        x_i_minus_1, x_i, x_i_plus_1,
        y_i_minus_1, f_i_minus_1,
        y_i, f_i,
        y_i_plus_1, f_i_plus_1 
    ):
        h_i = x_i_plus_1 - x_i
        h_i_minus_1 = x_i - x_i_minus_1
        
        self.alpha = h_i_minus_1 / h_i

        if (self.alpha != 1):
            print("alpha not 1", self.alpha)
        # print("alpha", self.alpha)

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