# the following is the quintics we got from solving the system with sympy
# we can now look to test these by seeing if they satisfy our conditions
def compare(computed_value, expected_value, for_which):
    passed = True
    for computed, expected in zip(computed_value, expected_value):
        if (abs(computed - expected) > 1e-12):
            print(f"FAIL: {for_which} by => {abs(computed-expected)}")
            passed = False
    if passed:
        print(f"PASS: {for_which}")
    else:    
        print(computed_value)
        print(expected_value)

# these are the expected values of the computations
b_for_d = [
    [1, 0, 0, 0, 0, 0], # b_for_d0 => d_0(-alpha) == 1 only
    [0, 0, 0, 1, 0, 0], # b_for_d1 => d_1_prime(-alpha) == 1 only
    [0, 1, 0, 0, 0, 0], # b_for_d2 => d_2(0) == 1 only
    [0, 0, 0, 0, 1, 0], # b_for_d3 => d_3_prime(0) == 1 only
    [0, 0, 1, 0, 0, 0], # b_for_d4 => d_4(1) == 1 only
    [0, 0, 0, 0, 0, 1]  # b_for_d5 => d_5_prime(1) == 1 only
]

for alpha in [1, 1/2, 2, 1/4, 4, 1/8, 8, 1/(2**7), 2**7, 1/(2**10), 2**10]:
    print(f"TESTING for alpha={alpha}")
    print("==============================================================================================")

    def d0(x):
        return (4*alpha + 2)/(alpha**6 + 3*alpha**5 + 3*alpha**4 + alpha**3)*(x**5) + (5*alpha**2 - 5*alpha - 4)/(alpha**6 + 3*alpha**5 + 3*alpha**4 + alpha**3)*(x**4) + (-10*alpha**2 - 2*alpha + 2)/(alpha**6 + 3*alpha**5 + 3*alpha**4 + alpha**3)*(x**3) + (5*alpha + 3)/(alpha**5 + 3*alpha**4 + 3*alpha**3 + alpha**2)*(x**2) + 0*x + 0
    def d1(x):
        return 1/(alpha**4 + 2*alpha**3 + alpha**2)*(x**5) + (alpha - 2)/(alpha**4 + 2*alpha**3 + alpha**2)*(x**4) + (1 - 2*alpha)/(alpha**4 + 2*alpha**3 + alpha**2)*(x**3) + 1/(alpha**3 + 2*alpha**2 + alpha)*(x**2) + 0*x + 0
    def d2(x):
        return (2*alpha - 2)/alpha**3*(x**5) + (4*alpha**2 - 7*alpha + 4)/alpha**3*(x**4) + (2*alpha**3 - 8*alpha**2 + 8*alpha - 2)/alpha**3*(x**3) + (-3*alpha**2 + 4*alpha - 3)/alpha**2*(x**2) + 0*x + 1
    def d3(x):
        return alpha**(-2)*(x**5) + (2*alpha - 2)/alpha**2*(x**4) + (alpha**2 - 4*alpha + 1)/alpha**2*(x**3) + (2 - 2*alpha)/alpha*(x**2) + 1*x + 0
    def d4(x):
        return (-2*alpha - 4)/(alpha**3 + 3*alpha**2 + 3*alpha + 1)*(x**5) + (-4*alpha**2 - 5*alpha + 5)/(alpha**3 + 3*alpha**2 + 3*alpha + 1)*(x**4) + (-2*alpha**3 + 2*alpha**2 + 10*alpha)/(alpha**3 + 3*alpha**2 + 3*alpha + 1)*(x**3) + (3*alpha**3 + 5*alpha**2)/(alpha**3 + 3*alpha**2 + 3*alpha + 1)*(x**2) + 0*x + 0
    def d5(x):
        return 1/(alpha**2 + 2*alpha + 1)*(x**5) + (2*alpha - 1)/(alpha**2 + 2*alpha + 1)*(x**4) + (alpha**2 - 2*alpha)/(alpha**2 + 2*alpha + 1)*(x**3) + -alpha**2/(alpha**2 + 2*alpha + 1)*(x**2) + 0*x + 0

    def d0_prime(x):
        return 5*(4*alpha + 2)/(alpha**6 + 3*alpha**5 + 3*alpha**4 + alpha**3)*(x**4) + 4*(5*alpha**2 - 5*alpha - 4)/(alpha**6 + 3*alpha**5 + 3*alpha**4 + alpha**3)*(x**3) + 3*(-10*alpha**2 - 2*alpha + 2)/(alpha**6 + 3*alpha**5 + 3*alpha**4 + alpha**3)*(x**2) + 2*(5*alpha + 3)/(alpha**5 + 3*alpha**4 + 3*alpha**3 + alpha**2)*x + 0
    def d1_prime(x):
        return 5*1/(alpha**4 + 2*alpha**3 + alpha**2)*(x**4) + 4*(alpha - 2)/(alpha**4 + 2*alpha**3 + alpha**2)*(x**3) + 3*(1 - 2*alpha)/(alpha**4 + 2*alpha**3 + alpha**2)*(x**2) + 2*1/(alpha**3 + 2*alpha**2 + alpha)*x + 0
    def d2_prime(x):
        return 5*(2*alpha - 2)/alpha**3*(x**4) + 4*(4*alpha**2 - 7*alpha + 4)/alpha**3*(x**3) + 3*(2*alpha**3 - 8*alpha**2 + 8*alpha - 2)/alpha**3*(x**2) + 2*(-3*alpha**2 + 4*alpha - 3)/alpha**2*x + 0
    def d3_prime(x):
        return 5*alpha**(-2)*(x**4) + 4*(2*alpha - 2)/alpha**2*(x**3) + 3*(alpha**2 - 4*alpha + 1)/alpha**2*(x**2) + 2*(2 - 2*alpha)/alpha*x + 1
    def d4_prime(x):
        return 5*(-2*alpha - 4)/(alpha**3 + 3*alpha**2 + 3*alpha + 1)*(x**4) + 4*(-4*alpha**2 - 5*alpha + 5)/(alpha**3 + 3*alpha**2 + 3*alpha + 1)*(x**3) + 3*(-2*alpha**3 + 2*alpha**2 + 10*alpha)/(alpha**3 + 3*alpha**2 + 3*alpha + 1)*(x**2) + 2*(3*alpha**3 + 5*alpha**2)/(alpha**3 + 3*alpha**2 + 3*alpha + 1)*x + 0
    def d5_prime(x):
        return 5*1/(alpha**2 + 2*alpha + 1)*(x**4) + 4*(2*alpha - 1)/(alpha**2 + 2*alpha + 1)*(x**3) + 3*(alpha**2 - 2*alpha)/(alpha**2 + 2*alpha + 1)*(x**2) + 2*-alpha**2/(alpha**2 + 2*alpha + 1)*x + 0

    ds = [d0, d1, d2, d3, d4, d5]
    d_primes = [d0_prime, d1_prime, d2_prime, d3_prime, d4_prime, d5_prime]

    # we try to get recompute b_for_d[i] to check if the system for each quintic is satisfied 
    # the first three evaluations are the evalutions of the quintic itself at -alpha, 0, 1
    # the next three evaluations are evaluations of its derivatives at -alpha, 0, 1
    computed_values = []
    for i in range(6):
        res = [
            ds[i](-alpha),
            ds[i](0),
            ds[i](1),
            d_primes[i](-alpha),
            d_primes[i](0),
            d_primes[i](1),
        ]
        computed_values.append(res)

    for i in range(len(b_for_d)):
        expected_value = b_for_d[i]
        computed_value = computed_values[i]
        compare(computed_value, expected_value, f"for b_for_d[{i}] at alpha={alpha}")
    print("\n\n\n")


