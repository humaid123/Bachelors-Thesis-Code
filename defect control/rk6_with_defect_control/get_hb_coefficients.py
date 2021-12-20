from sympy import Matrix, solve_linear_system, Indexed, symbols, init_printing


init_printing()

# we create the coefficients for a 2-step interpolant where the steps have sizes: [alpha * h, h] 
# => the first step has size alpha * h, the second step has size h
# we would like to get formulae where we can plug in 'alpha' and the formula will return
# the appropriate Hermite Birkhoff coeffficients
alpha = symbols("alpha")

# the first three rows are the evaluations of a*(x**5) + b*(x**4) + c*(x**3) + d*(x**2) + e*x + f at -alpha, 0, 1
# the next three rows are the evaluation of the derivative 5*a*(x**4) + 4*b*(x**3) + 3*c*(x**2) + 2*d*x + e at -alpha, 0, 1
# A contains the RHS of the system, x will contain the coefficients [a, b, c, d, e, d]
A = [
    [(-alpha)**5    , (-alpha)**4   , (-alpha)**3    ,  (-alpha)**2  , -alpha    , 1],  # quintic at -alpha
    [ 0             ,  0            ,  0             ,  0            ,  0        , 1],  # quintic at 0
    [ 1             ,  1            ,  1             ,  1            ,  1        , 1],  # quintic at 1
    [ 5*(-alpha)**4 , 4*(-alpha)**3 ,  3*(-alpha)**2 ,  2*(-alpha)   ,  1        , 0],  # quintic prime at -alpha
    [ 0             ,  0            ,  0             ,  0            ,  1        , 0],  # quintic prime at 0
    [ 5             ,  4            ,  3             ,  2            ,  1        , 0]   # quintic prime at 1
]

# we need to define the LHS, 'b' for the system for each of the quintics for the hermite birkhoff d0, d1, d2, d3, d4, d5, d6
# such that 
    # def hb(x):
    #       pheta = x - x_i / (x_i_plus_1 - x_i)
    #       return (  
    #              d0(pheta) * y_i_minus_1 + h_i * d1(pheta) * f_i_minus_1
    #            + d2(pheta) * y_i         + h_i * d3(pheta) * f_i 
    #            + d4(pheta) * y_i_plus_1  + h_i * d5(pheta) * f_i_plus_1
    #       )

# thus we get that 
# when pheta == -alpha
    # as u(t_i_minus_1) === y_i_minus_1,        =>    d0(pheta) == 1 and all the others are 0 
        # b_for_d0 equal 0 for every evaluation except for when quintic at -alpha => [1, 0, 0, 0, 0, 0]
    # as u_prime(t_i_minus_1) === f_i_minus_1   =>    d1_prime(pheta) == 1 and all the others are 0 
        # b_for_d1 equal 0 for every evaluation except for when quintic_prime at -alpha => [0, 0, 0, 1, 0, 0]

# when pheta == 0
    # as u(t_i) === y_i                         =>    d2(pheta) == 1 and all the others are 0 
        # b_for_d2 equal 0 for every evaluation except for when quintic at 0 => [0, 1, 0, 0, 0, 0]
    # as u_prime(t_i) === f_i                   =>    d3_prime(pheta) == 1, and all the others are 0 
        # b_for_d3 equal 0 for every evaluation except for when quintic prime at 0 => [0, 0, 0, 0, 1, 0]

# when pheta == 1
    # as u(t_i_plus_1) === y_i_plus_1           =>    d4(pheta) == 1 and all the others are 0 
        # b_for_d4 equal 0 for every evaluation except for when quintic at 1 => [0, 0, 1, 0, 0, 0]
    # as u_prime(t_i_plus_1) === f_i_plus_1     =>    d5_prime(pheta) == 1 and all the others are 0 
        # b_for_d5 equal 0 for every evaluation except for when quintic prime at 1 => [0, 0, 0, 0, 0, 1]

# the following matrix is such that
# b_for_d[i] is the RHS 'b', for d_i    
b_for_d = [
    [1, 0, 0, 0, 0, 0], # b_for_d0 => d_0(-alpha) == 1 only
    [0, 0, 0, 1, 0, 0], # b_for_d1 => d_1_prime(-alpha) == 1 only
    [0, 1, 0, 0, 0, 0], # b_for_d2 => d_2(0) == 1 only
    [0, 0, 0, 0, 1, 0], # b_for_d3 => d_3_prime(0) == 1 only
    [0, 0, 1, 0, 0, 0], # b_for_d4 => d_4(1) == 1 only
    [0, 0, 0, 0, 0, 1]  # b_for_d5 => d_5_prime(1) == 1 only
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
    print(f"d{i} has coefficients: ", res)
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
    print(f"def d{i}(x):\n    return", get_quintic_as_string(res[Indexed('a', i)], res[Indexed('b', i)], res[Indexed('c', i)], res[Indexed('d', i)], res[Indexed('e', i)], res[Indexed('f', i)]))

for i, res in enumerate(ans):
    print(f"def d{i}_prime(x):\n    return", get_quintic_primes_as_string(res[Indexed('a', i)], res[Indexed('b', i)], res[Indexed('c', i)], res[Indexed('d', i)], res[Indexed('e', i)], res[Indexed('f', i)]))