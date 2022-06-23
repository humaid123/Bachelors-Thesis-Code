from rk6 import rk_error_control_static_alpha_beta
from math import exp, floor
import matplotlib.pyplot as plt

def create_t_eval(start, end, num_points = 100):
    res = [start]
    h = (end - start)/ (num_points - 1)

    for _ in range(num_points - 1):
        res.append(
            res[-1] + h
        )
    return res

t_span_3 = [0, 10]
y0_3 = [1]

def model3(t, y):
    return [(1/4)*y*(1 - y/20)]

def solution3(t):
    return [20 / ( 1 + 19 * exp(-x/4) ) for x in t]

t_span_2 = [0, 10]
y0_2 = [1]

def model2(t, y):
    return [-2*t*y**2]

def solution2(t):
    return [1/(1+x**2) for x in t]


def experiment(model, solution, y0, t_span):
    t_eval = create_t_eval(t_span[0], t_span[1])
    actual_solutions = solution(t_eval)

    tols = [1e-8, 1e-9, 1e-10]
    for tol in tols:
        res, sol, first_deriv, error_samplings = rk_error_control_static_alpha_beta(model, t_span, y0[0], tol, solution)
        computed_solutions = [sol(x) for x in t_eval]
         
        plt.figure()
        plt.plot(t_eval, actual_solutions, label="actual solution")
        plt.plot(t_eval, computed_solutions, label="computed solution")
        #for (x, y) in res:
        #    plt.axvline(x=x)
        plt.title(f"solutions for tol={tol}")
        plt.legend()
        plt.show()


        error = [abs(actual_solution - computed_solution) for actual_solution, computed_solution in zip(actual_solutions, computed_solutions)]
        plt.figure()
        plt.plot(t_eval, error)
        plt.title(f"error for tol={tol}")
        plt.show()

        # shape of errors - graphs
        plt.figure()
        for (x_i, x_i_plus_1, hb) in error_samplings:
            num_points = 100
            pts_to_sample = create_t_eval(x_i, x_i_plus_1, num_points)
            errors = []
            for i, pt in enumerate(pts_to_sample):
                y = solution([pt])[0]
                hb_eval = hb.eval(pt)
                error = abs(hb_eval - y)
                errors.append( error )

                # print the error at the extremities
                if i == 0:
                    interpolation_error = hb_eval - hb.y_i
                    print("error=", error, "interpolation_error=", interpolation_error)

                if i == len(pts_to_sample) - 1:
                    interpolation_error = hb_eval - hb.y_i_plus_1
                    print("error=", error, "interpolation_error=", interpolation_error)

            maximum_error = max(errors)
            scaled_errors = [error / (maximum_error) for error in errors]

            # str_x_i = "{:.3f}".format(x_i)
            # str_x_i_plus_1 = "{:.3f}".format(x_i_plus_1)
            x_axis = [i/(num_points - 1) for i in range(num_points)]
            plt.plot(x_axis, scaled_errors, label=f"x_{str(x_i)}_{str(x_i_plus_1)}")
        plt.title("plot of shape of errors")
        plt.xlabel("step scaled between 0 and 1")
        plt.ylabel('scaled error')
        # plt.legend()
        plt.show()


        # actual_f_evals = [model(x, solution([x])[0] )[0] for x in t_eval]
        # hb_prime_evals = [first_deriv(x) for x in t_eval]
        # plt.figure()
        # plt.plot(t_eval, actual_f_evals)
        # plt.plot(t_eval, hb_prime_evals)
        # plt.title(f"first derivative for tol={tol}")
        # plt.show()

        # defects = [abs(actual_f_eval - hb_prime_eval) for (actual_f_eval, hb_prime_eval) in zip(actual_f_evals, hb_prime_evals)]
        # plt.figure()
        # plt.plot(t_eval, defects)
        # plt.title(f"global defect for tol={tol}")
        # plt.show()
        
        # # defect graphs
        # # we only pick 10 of the derivs
        # pick_every = (len(derivs) // 5) + 1
        # # print("pick every", pick_every)
        # plotted_derivs = []
        # for i in range(0, len(derivs), pick_every):
        #     plotted_derivs.append(
        #         derivs[i]
        #     )

        # plt.figure()
        # for (x_i_minus_1, x_i_plus_1, hb) in plotted_derivs:
        #     num_points = 100
        #     pts_to_sample = create_t_eval(x_i_minus_1, x_i_plus_1, num_points)
        #     defects = []
        #     for pt in pts_to_sample:
        #         y = solution([pt])[0]
        #         f_eval  = model(pt, y)[0]
        #         hb_prime_eval = hb.prime(pt)
        #         defects.append( abs(hb_prime_eval - f_eval) )
        #     maximum_defect = max(defects)
        #     minimum_defect = min(defects)
        #     plot_vals = [(defect - minimum_defect) / (maximum_defect - minimum_defect) for defect in defects]
        #     #plt.plot(xs, defects, label=f"x_{str(x_i_minus_1)}_{str(x_i_plus_1)}")

        #     x_axis = [i/(num_points - 1) for i in range(num_points)]
        #     str_x_i_minus_1 = "{:.3f}".format(x_i_minus_1)
        #     str_x_i_plus_1 = "{:.3f}".format(x_i_plus_1)
        #     plt.plot(x_axis, plot_vals, label=f"x_{str_x_i_minus_1}_{str_x_i_plus_1}")
        # plt.title("plot of defects")
        # plt.xlabel("step scaled between 0 and 2")
        # plt.ylabel('defect')
        # plt.legend()
        # plt.show()

# experiment(model3, solution3, y0_3, t_span_3)
experiment(model2, solution2, y0_2, t_span_2)