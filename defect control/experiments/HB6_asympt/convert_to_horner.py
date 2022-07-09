

from sympy.polys.polyfuncs import horner
from sympy.abc import x
from sympy import symbols

alpha = symbols("alpha")
beta = symbols("beta")


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
print("def d0_horner(x, alpha, beta) -> float:\n    return (", horner(d0(x, alpha, beta)), ")\n")
print("def d1_horner(x, alpha, beta) -> float:\n    return (", horner(d1(x, alpha, beta)), ")\n")
print("def d2_horner(x, alpha, beta) -> float:\n    return (", horner(d2(x, alpha, beta)), ")\n")
print("def d3_horner(x, alpha, beta) -> float:\n    return (", horner(d3(x, alpha, beta)), ")\n")
print("def d4_horner(x, alpha, beta) -> float:\n    return (", horner(d4(x, alpha, beta)), ")\n")
print("def d5_horner(x, alpha, beta) -> float:\n    return (", horner(d5(x, alpha, beta)), ")\n")


print("def d0_prime_horner(x, alpha, beta) -> float:\n    return (", horner(d0_prime(x, alpha, beta)), ")\n")
print("def d1_prime_horner(x, alpha, beta) -> float:\n    return (", horner(d1_prime(x, alpha, beta)), ")\n")
print("def d2_prime_horner(x, alpha, beta) -> float:\n    return (", horner(d2_prime(x, alpha, beta)), ")\n")
print("def d3_prime_horner(x, alpha, beta) -> float:\n    return (", horner(d3_prime(x, alpha, beta)), ")\n")
print("def d4_prime_horner(x, alpha, beta) -> float:\n    return (", horner(d4_prime(x, alpha, beta)), ")\n")
print("def d5_prime_horner(x, alpha, beta) -> float:\n    return (", horner(d5_prime(x, alpha, beta)), ")\n")
