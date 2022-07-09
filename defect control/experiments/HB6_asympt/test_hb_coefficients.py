from HB6_asympt import d0, d1, d2, d3, d4, d5, d0_prime, d1_prime, d2_prime, d3_prime, d4_prime, d5_prime


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
    [0, 0,        0, 0, 0, 1], # d0_prime(gamma) == 1
    [0, 1,        0, 0, 0, 0], # d1(-1) == 1
    [0, 0,        0, 0, 1, 0], # d2_prime(-1) == 1
    [1, 0,        0, 0, 0, 0], # d3(0) == 1
    [0, 0,        1, 0, 0, 0], # d4_prime(0) == 1
    [0, 0,        0, 1, 0, 0], # d5_prime(beta) == 1
]

# values = [1, 1/2, 2, 1/4, 4, 1/8, 8, 1/(2**7), 2**7, 1/(2**10), 2**10]
values = [1, 1/2, 2, 1/4, 4]
for alpha in values:
    for beta in values:
        print(f"TESTING for alpha={alpha}, beta={beta}")
        print("==============================================================================================")

        ds = [d0, d1, d2, d3, d4, d5]
        d_primes = [d0_prime, d1_prime, d2_prime, d3_prime, d4_prime, d5_prime]

        gamma = - (1 + alpha)
        # we try to get recompute b_for_d[i] to check if the system for each quintic is satisfied 
        # the first three evaluations are the evalutions of the quintic itself at -alpha, 0, 1
        # the next three evaluations are evaluations of its derivatives at -alpha, 0, 1
        computed_values = []
        for i in range(len(ds)):
            res = [
                ds[i](0, alpha, beta),
                ds[i](-1, alpha, beta),

                d_primes[i](0, alpha, beta),
                d_primes[i](beta, alpha, beta),
                d_primes[i](-1, alpha, beta),
                d_primes[i](gamma, alpha, beta),
            ]
            computed_values.append(res)

        for i in range(len(b_for_d)):
            expected_value = b_for_d[i]
            computed_value = computed_values[i]
            compare(computed_value, expected_value, f"for b_for_d[{i}] at alpha={alpha}")
        print("\n\n\n")


