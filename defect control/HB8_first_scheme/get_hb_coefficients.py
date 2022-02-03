from sympy import Matrix, solve_linear_system, Indexed, symbols, init_printing


init_printing()

# we create the coefficients for a 3-step interpolant where the steps have sizes: [alpha * h, beta * h, h] 
# we want formulae in alpha and beta for our septics d0, d1, d2, d3, d4, d5, d6, d7
beta = symbols("beta")
alpha = symbols("alpha")
gamma = - alpha - beta 

# the first three rows are the evaluations of a*(x**7) + b*(x**6) + c*(x**5) + d*(x**4) + e*(x**3) + f(x**2) + gx + h at 1, 0, -alpha, -alpha-beta
# the next three rows are the evaluations of 7a*(x**6) + 6b*(x**5) + 5c*(x**4) + 4d*(x**3) + 3e*(x**2) + 2fx + g at 1, 0, -alpha, -alpha-beta
# A contains the RHS of the system, x will contain the coefficients [a, b, c, d, e, f, g, h]
A = [
    [1                 ,  1                ,  1                ,  1                ,  1                , 1            ,  1       ,     1],  # septic at 1
    [0                 ,  0                ,  0                ,  0                ,  0                , 0            ,  0       ,     1],  # septic at 0
    [(-alpha)**7       ,  (-alpha)**6      ,  (-alpha)**5      ,  (-alpha)**4      ,  (-alpha)**3      , (-alpha)**2  ,  -alpha  ,     1],  # septic at -alpha
    [gamma**7          ,  gamma**6         ,  gamma**5         ,  gamma**4         ,  gamma**3         , gamma**2     ,  gamma   ,     1],  # septic at gamma = -alpha-beta
    [ 7                ,  6                ,  5                ,  4                ,  3                , 2            ,  1       ,     0],  # septic prime at 1
    [ 0                ,  0                ,  0                ,  0                ,  0                , 0            ,  1       ,     0],  # septic prime at 0
    [ 7*((-alpha)**6)  ,  6*((-alpha)**5)  ,  5*((-alpha)**4)  ,  4*((-alpha)**3)  ,  3*((-alpha)**2)  , 2*(-alpha)   ,  1       ,     0],  # septic prime at -alpha
    [ 7*(gamma**6)     ,  6*(gamma**5)     ,  5*(gamma**4)     ,  4*(gamma**3)     ,  3*(gamma**2)     , 2*gamma      ,  1       ,     0],  # septic prime at gamma = -alpha-beta
]

# we need to define the LHS, 'b' for the system for each of the quintics for the hermite birkhoff d0, d1, d2, d3, d4, d5, d6
# such that 
    # def hb(x):
    #       pheta = (x - x_i) / (x_i_plus_1 - x_i) 
    #       return (  
    #              d0(pheta) * y_i_minus_2 + h_i * d1(pheta) * f_i_minus_2
    #            + d2(pheta) * y_i_minus_1 + h_i * d3(pheta) * f_i_minus_1
    #            + d4(pheta) * y_i         + h_i * d5(pheta) * f_i 
    #            + d6(pheta) * y_i_plus_1  + h_i * d7(pheta) * f_i_plus_1
    #       )

#[beta * h, alpha*h, h]
# thus we get that 
# at each of pheta == gamma, -alpha, 0, 1, ONLY one of d0, d2, d4, d6 evaluate to 1, everything else evaluate to 0
                                # only one of d1_prime, d3_prime, d5_prime, d7_prime evaluate to 1, everything else evaluate to 1
# thus each of b_for_d_i is an 8-length vector that is 0 for each except for one value where it is 1
    # as each b_for_d_i is 0 at each eval of the matrix A except for 1

# when pheta == 0
    # as u(t_i) === y_i                         =>    d4(pheta) == 1 and all the others are 0 
        # b_for_d4 equal 0 for every evaluation in A except for when septic at 0 => [0, 1, 0, 0, 0, 0, 0, 0]
    # as u_prime(t_i) === f_i                   =>    d5_prime(pheta) == 1, and all the others are 0 
        # b_for_d3 equal 0 for every evaluation in A except for when septic prime at 0 => [0, 0, 0, 0, 0, 1, 0, 0]

# when pheta == 1
    # as u(t_i_plus_1) === y_i_plus_1           => d6(pheta) == 1 and all others are 0
        # b_for_d6 equal 0 for every evaluation in A except for when septic at 1 => [1, 0, 0 , 0, 0, 0, 0, 0]
    # as u_prime(t_i_plus_1) === f_i_plus_1     => d7(pheta) == 1 and all others are 0
        # b_for_d7 equal 0 for every evaluation in A except for when septic prime at 1 => [0, 0, 0, 0, 1, 0, 0, 0]

# when pheta == -alpha
    # as u(t_i_minus_1) === y_i_minus_1           =>    d2(pheta) == 1 and all the others are 0 
        # b_for_d2 equal 0 for every evaluation in A except for when septic at -alpha => [0, 0, 1, 0, 0, 0, 0, 0]
    # as u_prime(t_i_plus_1) === f_i_plus_1     =>    d3_prime(pheta) == 1 and all the others are 0 
        # b_for_d3 equal 0 for every evaluation in A except for when septic prime at -alpha => [0, 0, 0, 0, 0, 0, 1, 0]

# when pheta == -alpha-gamma
    # as u(t_i_minus_2) === y_i_minus_2           =>   d0(pheta) == 1 and all others are 0
        # b_for_d0 equal 0 for every evaluation in A except when septic at -alpha-gamma => [0, 0, 0, 1, 0, 0, 0, 0]
    # as u_prime(t_i_minus_2) === f_i_minus_2     =>   d1(pheta) == 1 and all others are 0
        # b_for_d1 equal 0 for every evaluation in A except when septic prime at -alpha-gamma => [0, 0, 0, 0, 0, 0, 0, 1]

# the following matrix is such that
# b_for_d[i] is the RHS 'b', for d_i    
b_for_d = [
    [0, 0, 0, 1, 0, 0, 0, 0], # b_for_d0 => d0(gamma) == 1
    [0, 0, 0, 0, 0, 0, 0, 1], # b_for_d1 => d1_prime(gamma) == 1
    [0, 0, 1, 0, 0, 0, 0, 0], # b_for_d2 => d2(-alpha) == 1
    [0, 0, 0, 0, 0, 0, 1, 0], # b_for_d3 => d3_prime(-alpha) == 1
    [0, 1, 0, 0, 0, 0, 0, 0], # b_for_d4 => d4(0) == 1
    [0, 0, 0, 0, 0, 1, 0, 0], # b_for_d5 => d5_prime(0) == 1
    [1, 0, 0, 0, 0, 0, 0, 0], # b_for_d6 => d6(1) == 1
    [0, 0, 0, 0, 1, 0, 0, 0], # b_for_d7 => d7_prime(1) == 1
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
    res = solve_linear_system(system, Indexed('a', i), Indexed('b', i), Indexed('c', i), Indexed('d', i), Indexed('e', i), Indexed('f',i), Indexed("g", i), Indexed("h", i))
    ans[i] = res

# printing the solution as python functions
# ans is an array of dictionaries
#       Each element of ans is a dictionary that contains the coefficients for each that prime
#       the dictionary contains keys a_i, b_i, c_i, d_i, e_i, f_i and the values at those keys are STRINGs in alpha representing the formula to get the value of this coefficent at a given alpha
#       to get a_i, we need to use Indexed('a', i)
def get_septic_as_string(a, b, c, d, e, f, g, h):
    return (f"{a}*(x**7) + {b}*(x**6) + {c}*(x**5) + {d}*(x**4) + {e}*(x**3) + {f}*(x**2) + {g}*(x) + {h}")
def get_septic_primes_as_string(a, b, c, d, e, f, g, h):
    return f"7*{a}*(x**6) + 6*{b}*(x**5) + 5*{c}*(x**4) + 4*{d}*(x**3) + 3*{e}*(x**2) + 2*{f}*x + {g}"

for i, res in enumerate(ans):
    print(f"def d{i}(x, alpha, beta) -> float:\n    return (", get_septic_as_string(res[Indexed('a', i)], res[Indexed('b', i)], res[Indexed('c', i)], res[Indexed('d', i)], res[Indexed('e', i)], res[Indexed('f', i)], res[Indexed("g", i)], res[Indexed("h", i)]), ")")

for i, res in enumerate(ans):
    print(f"def d{i}_prime(x, alpha, beta) -> float:\n    return (", get_septic_primes_as_string(res[Indexed('a', i)], res[Indexed('b', i)], res[Indexed('c', i)], res[Indexed('d', i)], res[Indexed('e', i)], res[Indexed('f', i)], res[Indexed('g', i)], res[Indexed('h', i)]), ")")