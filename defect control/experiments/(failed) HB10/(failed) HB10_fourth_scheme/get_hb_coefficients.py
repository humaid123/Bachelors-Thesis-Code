from sympy import Matrix, solve_linear_system, Indexed, symbols, init_printing


init_printing()

# we create the coefficients for a 4-step interpolant where the steps have sizes: [alpha*h, h, h, beta*h] 
# we are going to cut the step from x_i-1 to x_i in half and use x_i-0.5 to x_i as the base step

# NOTE:
# THE DIFFERENCE BETWEEN THIS SCHEME AND the first scheme is that we center around x_i_minus_0_5
# so pheta = (x - x_i_minus_0_5) / h        where h = x_i_minus_0_5 - x_i_minus_1

# we want formulae in alpha and beta for our 9 nonics d0, d1, d2, d3, d4, d5, d6, d7, d8, d9
beta = symbols("beta")
alpha = symbols("alpha")

gamma = 1 + beta 
tau = -1 - alpha

# the first three rows are the evaluations of a*(x**7) + b*(x**6) + c*(x**5) + d*(x**4) + e*(x**3) + f(x**2) + gx + h at 0, beta, -1, -(1 + alpha)
# the next three rows are the evaluations of 7a*(x**6) + 6b*(x**5) + 5c*(x**4) + 4d*(x**3) + 3e*(x**2) + 2fx + g at 0, beta, -1, -(1 + alpha)
# A contains the RHS of the system, x will contain the coefficients [a, b, c, d, e, f, g, h]
A = [
    [0              , 0             , 0              , 0             , 0              ,  0            ,  0            ,       0  ,     0  ,     1],  # nonic at 0
    [1              , 1             , 1              ,  1            ,  1             ,  1            ,  1            , 1        ,     1  ,     1],  # nonic at 1
    [gamma**9       , gamma**8      , gamma**7       ,  gamma**6     ,  gamma**5      ,  gamma**4     ,  gamma**3     , gamma**2 ,  gamma ,     1],  # nonic at gamma = 1 + beta
    [-1             , 1             , -1             ,  1            ,  -1            ,  1            ,  -1           , 1        ,    -1  ,     1],  # nonic at -1
    [tau**9         , tau**8        , tau**7         ,  tau**6       ,  tau**5        ,  tau**4       ,  tau**3       , tau**2   ,  tau   ,     1],  # nonic at tau = - 1 - alpha

    [0              , 0             , 0              ,  0            ,  0             ,  0            ,  0            , 0        ,  1     ,     0],  # nonic prime at 0
    [ 9             , 8             , 7              ,  6            ,  5             ,  4            ,  3            , 2        ,  1     ,     0],  # nonic prime at 1
    [ 9*(gamma**8)  , 8*(gamma**7)  , 7*(gamma**6)   ,  6*(gamma**5) ,  5*(gamma**4)  ,  4*(gamma**3) ,  3*(gamma**2) , 2*gamma  ,  1     ,     0],  # nonic prime at gamma = 1 + beta
    [ 9             , -8            , 7              ,  -6           ,  5             ,  -4           ,  3            , -2       ,  1     ,     0],  # nonic prime at -1
    [ 9*(tau**8)    , 8*(tau**7)    , 7*(tau**6)     ,  6*(tau**5)   ,  5*(tau**4)    ,  4*(tau**3)   ,  3*(tau**2)   , 2*tau    ,  1     ,     0],  # nonic prime at tau = - 1 - alpha
]

# we need to define the LHS, 'b' for the system for each of the quintics for the hermite birkhoff d0, d1, d2, d3, d4, d5, d6
# such that 
    # def hb(x):
    #       pheta = x - x_i / (x_i - x_i_minus_1)  # h is between x_i and x_i_minus_1 now
    #       return (  
    #              d0(pheta) * y_i_minus_2   + h_i * d1(pheta) * f_i_minus_2
    #            + d2(pheta) * y_i_minus_1   + h_i * d3(pheta) * f_i_minus_1
    #            + d4(pheta) * y_i_minus_0.5 + h_i * d5(pheta) * f_i_minus_0.5
    #            + d6(pheta) * y_i           + h_i * d7(pheta) * f_i 
    #            + d8(pheta) * y_i_plus_1    + h_i * d9(pheta) * f_i_plus_1
    #       )

#[alpha * h, h, h, beta * h]
# thus we get that 
"""
when pheta == 0
    d4(0) == 1              [1 0 0 0 0     0 0 0 0 0]
    d5_prime(0) == 1        [0 0 0 0 0     1 0 0 0 0]

when pheta == 1
    d6(1) == 1              [0 1 0 0 0     0 0 0 0 0]
    d7_prime(1) == 1        [0 0 0 0 0     0 1 0 0 0]

when pheta == gamma == 1 + beta
    d8(gamma) == 1          [0 0 1 0 0     0 0 0 0 0]
    d9_prime(gamma) == 1    [0 0 0 0 0     0 0 1 0 0]

when pheta == -1
    d2(-1) == 1             [0 0 0 1 0     0 0 0 0 0]
    d3_prime(-1) == 1       [0 0 0 0 0     0 0 0 1 0]

when pheta == -1 - alpha == tau
    d0(tau) == 1            [0 0 0 0 1     0 0 0 0 0]
    d1(tau) == 1            [0 0 0 0 0     0 0 0 0 1]

"""



# the following matrix is such that
# b_for_d[i] is the RHS 'b', for d_i    
b_for_d = [
    [0, 0, 0, 0, 1,     0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0,     0, 0, 0, 0, 1],

    [0, 0, 0, 1, 0,     0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0,     0, 0, 0, 1, 0],

    [1, 0, 0, 0, 0,     0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0,     1, 0, 0, 0, 0],

    [0, 1, 0, 0, 0,     0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0,     0, 1, 0, 0, 0],
    
    [0, 0, 1, 0, 0,     0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0,     0, 0, 1, 0, 0],
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
    res = solve_linear_system(system, Indexed('a', i), Indexed('b', i), Indexed('c', i), Indexed('d', i), Indexed('e', i), Indexed('f',i), Indexed("g", i), Indexed("h", i), Indexed("i", i), Indexed("j", i))
    ans[i] = res

# printing the solution as python functions
# ans is an array of dictionaries
#       Each element of ans is a dictionary that contains the coefficients for each that prime
#       the dictionary contains keys a_i, b_i, c_i, d_i, e_i, f_i and the values at those keys are STRINGs in alpha representing the formula to get the value of this coefficent at a given alpha
#       to get a_i, we need to use Indexed('a', i)
def get_nonic_as_string(a, b, c, d, e, f, g, h, i, j):
    return (f"{a}*(x**9) + {b}*(x**8) + {c}*(x**7) + {d}*(x**6) + {e}*(x**5) + {f}*(x**4) + {g}*(x**3) + {h}*(x**2) + {i}*(x) + {j}")
def get_nonic_primes_as_string(a, b, c, d, e, f, g, h, i, j):
    return f"9*{a}*(x**8) + 8*{b}*(x**7) + 7*{c}*(x**6) + 6*{d}*(x**5) + 5*{e}*(x**4) + 4*{f}*(x**3) + 3*{g}*(x**2) + 2*{h}*(x) + {i}"

for i, res in enumerate(ans):
    print(f"def d{i}(x, alpha, beta):\n    return (", get_nonic_as_string(res[Indexed('a', i)], res[Indexed('b', i)], res[Indexed('c', i)], res[Indexed('d', i)], res[Indexed('e', i)], res[Indexed('f', i)], res[Indexed("g", i)], res[Indexed("h", i)], res[Indexed("i", i)], res[Indexed("j", i)]   ), ")")

for i, res in enumerate(ans):
    print(f"def d{i}_prime(x, alpha, beta):\n    return (", get_nonic_primes_as_string(res[Indexed('a', i)], res[Indexed('b', i)], res[Indexed('c', i)], res[Indexed('d', i)], res[Indexed('e', i)], res[Indexed('f', i)], res[Indexed('g', i)], res[Indexed('h', i)], res[Indexed("i", i)], res[Indexed("j", i)] ), ")")