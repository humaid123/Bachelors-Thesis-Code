

from sympy.polys.polyfuncs import horner
from sympy.abc import x
from sympy import symbols

alpha = symbols("alpha")


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


print("def d0_horner(x, alpha) -> float:\n    return (", horner(d0(x, alpha)), ")\n")
print("def d1_horner(x, alpha) -> float:\n    return (", horner(d1(x, alpha)), ")\n")
print("def d2_horner(x, alpha) -> float:\n    return (", horner(d2(x, alpha)), ")\n")
print("def d3_horner(x, alpha) -> float:\n    return (", horner(d3(x, alpha)), ")\n")
print("def d4_horner(x, alpha) -> float:\n    return (", horner(d4(x, alpha)), ")\n")
print("def d5_horner(x, alpha) -> float:\n    return (", horner(d5(x, alpha)), ")\n")

print("def d0_prime_horner(x, alpha) -> float:\n    return (", horner(d0_prime(x, alpha)), ")\n")
print("def d1_prime_horner(x, alpha) -> float:\n    return (", horner(d1_prime(x, alpha)), ")\n")
print("def d2_prime_horner(x, alpha) -> float:\n    return (", horner(d2_prime(x, alpha)), ")\n")
print("def d3_prime_horner(x, alpha) -> float:\n    return (", horner(d3_prime(x, alpha)), ")\n")
print("def d4_prime_horner(x, alpha) -> float:\n    return (", horner(d4_prime(x, alpha)), ")\n")
print("def d5_prime_horner(x, alpha) -> float:\n    return (", horner(d5_prime(x, alpha)), ")\n")