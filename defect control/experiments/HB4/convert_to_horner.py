
from sympy.polys.polyfuncs import horner
from sympy.abc import x

def h00(t):
    return 2*(t**3) - 3*(t**2) + 1

def h10(t):
    return t**3 - 2*(t**2) + t

def h01(t):
    return -2*(t**3) + 3*(t**2)

def h11(t):
    return t**3 - t**2

def h00_prime(t):
    return 6*(t**2) - 6*t

def h10_prime(t):
    return 3*(t**2) - 4*t + 1

def h01_prime(t):
    return -6*(t**2) + 6*t

def h11_prime(t):
    return 3*(t**2) - 2*t

print("def h00_horner(x):\n    return (", horner(h00(x)), ")")
print("def h10_horner(x):\n    return (", horner(h10(x)), ")")
print("def h01_horner(x):\n    return (", horner(h01(x)), ")")
print("def h11_horner(x):\n    return (", horner(h11(x)), ")")

print("def h00_prime_horner(x):\n    return (", horner(h00_prime(x)), ")")
print("def h10_prime_horner(x):\n    return (", horner(h10_prime(x)), ")")
print("def h01_prime_horner(x):\n    return (", horner(h01_prime(x)), ")")
print("def h11_prime_horner(x):\n    return (", horner(h11_prime(x)), ")")