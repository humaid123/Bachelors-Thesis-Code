from sympy import Matrix, solve_linear_system, Indexed, symbols, init_printing


init_printing()

# we create the coefficients for a 3-step interpolant where the steps have sizes: [alpha * h, h, beta*h] 
# we want formulae in alpha and beta for our septics d0, d1, d2, d3, d4, d5, d6, d7
beta = symbols("beta")
gamma = - 1 - symbols("alpha")

# the first three rows are the evaluations of a*(x**7) + b*(x**6) + c*(x**5) + d*(x**4) + e*(x**3) + f(x**2) + gx + h at 0, beta, -1, -(1 + alpha)
# the next three rows are the evaluations of 7a*(x**6) + 6b*(x**5) + 5c*(x**4) + 4d*(x**3) + 3e*(x**2) + 2fx + g at 0, beta, -1, -(1 + alpha)
# A contains the RHS of the system, x will contain the coefficients [a, b, c, d, e, f, g, h]
A = [
    [  0             ,  0            ,  0            ,       0  ,     0  ,     1],  # quintic at 0
    [  -1            ,  1            ,  -1           , 1        ,    -1  ,     1],  # quintic at -1
    [  0             ,  0            ,  0            , 0        ,  1     ,     0],  # quintic prime at 0
    [  5*(beta**4)   ,  4*(beta**3)  ,  3*(beta**2)  , 2*beta   ,  1     ,     0],  # quintic prime at beta
    [  5             ,  -4           ,  3            , -2       ,  1     ,     0],  # quintic prime at -1
    [  5*(gamma**4)  ,  4*(gamma**3) ,  3*(gamma**2) , 2*gamma  ,  1     ,     0],  # quintic prime at gamma = - (1 + alpha)
]

# we need to define the LHS, 'b' for the system for each of the quintics for the hermite birkhoff d0, d1, d2, d3, d4, d5, d6
# such that 
    # def hb(x):
    #       pheta = x - x_i / (x_i - x_i_minus_1)  # h is between x_i and x_i_minus_1 now
    #       return (  
    #                                        h_i * d0(pheta) * f_i_minus_2
    #            + d1(pheta) * y_i_minus_1 + h_i * d2(pheta) * f_i_minus_1
    #            + d3(pheta) * y_i         + h_i * d4(pheta) * f_i 
    #                                      + h_i * d5(pheta) * f_i_plus_1
    #       )

#[alpha * h, h, beta * h]
# thus we get that 
# at pheta == gamma, -1, 0, alpha, ONLY one of d0, d2, d4 evaluate to 1, everything else evaluate to 0
                                # only one of d1_prime, d3_prime, d5_prime evaluate to 1, everythin else evaluate to 1
# thus each of b_for_d_i is an 6-length vector that is 0 for each except for one value where it is 1
    # as each b_for_d_i is 0 at each eval of the matrix A except for 1

"""
when pheta == 0:
    d3(0) = 1           [1 0        0 0 0 0]
    d4_prime(0) = 1     [0 0        1 0 0 0]

when pheta == beta
    d5_prime(beta) = 1  [0 0        0 1 0 0]

when pheta == -1
    d1(-1) = 1          [0 1        0 0 0 0]
    d2_prime(-1) = 1    [0 0        0 0 1 0]

when pheta == gamma
    d0_prime(gamma) = 1 [0 0        0 0 0 1]

"""

# the following matrix is such that
# b_for_d[i] is the RHS 'b', for d_i    
b_for_d = [
    [0, 0,        0, 0, 0, 1], # d0_prime(gamma) == 1
    [0, 1,        0, 0, 0, 0], # d1(-1) == 1
    [0, 0,        0, 0, 1, 0], # d2_prime(-1) == 1
    [1, 0,        0, 0, 0, 0], # d3(0) == 1
    [0, 0,        1, 0, 0, 0], # d4_prime(0) == 1
    [0, 0,        0, 1, 0, 0], # d5_prime(beta) == 1
]

# sympy Matrix expects the system as the augmented matrix
# to solve ax + by = c and px + qy = r
# we need to provide it a matrix as such:
# matrix = [
#    [a, b, c]
#    [p, q, r]
# ]
# reference: https://stackoverflow.com/questions/31547657/how-can-i-solve-system-of-linear-equations-in-sympy
def create_system(A, b):
    this_A = [row[:] for row in A] # need to make a deep copy of A as we change it here
    for i in range(len(this_A)):
        this_A[i].append(b[i])
    # for row in this_A:
    #    print(row)
    return this_A



ans = [0] * len(b_for_d)        # create an array to store the results of solving the system
# solve the systems
for i in range(len(b_for_d)):
    system = Matrix(create_system(A, b_for_d[i]))
    res = solve_linear_system(system, Indexed('a', i), Indexed('b', i), Indexed('c', i), Indexed('d', i), Indexed('e', i), Indexed('f',i))
    ans[i] = res

# printing the solution as python functions
# ans is an array of dictionaries
#       Each element of ans is a dictionary that contains the coefficients for each that prime
#       the dictionary contains keys a_i, b_i, c_i, d_i, e_i, f_i and the values at those keys are STRINGs in alpha representing the formula to get the value of this coefficent at a given alpha
#       to get a_i, we need to use Indexed('a', i)
def get_quintic_as_string(a, b, c, d, e, f):
    return (f"{a}*(x**5) + {b}*(x**4) + {c}*(x**3) + {d}*(x**2) + {e}*x + {f}")
def get_quintic_primes_as_string(a, b, c, d, e, f):
    return f"5*{a}*(x**4) + 4*{b}*(x**3) + 3*{c}*(x**2) + 2*{d}*x + {e}"

for i, res in enumerate(ans):
    print(f"def d{i}(x, alpha, beta):\n    return (", get_quintic_as_string(res[Indexed('a', i)], res[Indexed('b', i)], res[Indexed('c', i)], res[Indexed('d', i)], res[Indexed('e', i)], res[Indexed('f', i)]), ")")

for i, res in enumerate(ans):
    print(f"def d{i}_prime(x, alpha, beta):\n    return (", get_quintic_primes_as_string(res[Indexed('a', i)], res[Indexed('b', i)], res[Indexed('c', i)], res[Indexed('d', i)], res[Indexed('e', i)], res[Indexed('f', i)]), ")")