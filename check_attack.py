"""
An attack on a key exchange protocol based on max-times and min-times
algebras from [M. I. Durcheva, An application of different dioids in public
key cryptography. In AIP Conference Proceedings, vol. 1631, pp. 336-343,
AIP, 2014].

I. Buchinskiy, M. Kotov, A. Treier, 2022
"""

import argparse
from attack import attack
from generate_instance import generate_random_instance


def check_attack(count, n, c_bound, d_bound, p_bound, t_bound):
    """
    Generates and runs instances to check the attack.
    c_bound is the upper bound for coefficients of matrices and polynomials.
    d_bound is the apper bound for degrees of polynomials.
    p_bound and t_bound are the bound to search polynomials p' and t' respectively (also, q' and r').
    """
    failed = 0
    incorrect = 0

    for i in range(count):
        print(i, end=" ")
        inst = generate_random_instance(n, c_bound, d_bound)

        k1 = attack(inst.M, inst.N, inst.X, inst.A, inst.B, p_bound, t_bound)

        if not k1:
            print("FAILED")
            print("M =", inst.M)
            print("N =", inst.N)
            print("X =", inst.X)
            print("p =", inst.p)
            print("t =", inst.t)
            print("q =", inst.q)
            print("r =", inst.r)
            failed += 1
        elif k1 != inst.kA:
            print("INCORRECT")
            print("M =", inst.M)
            print("N =", inst.N)
            print("X =", inst.X)
            print("p =", inst.p)
            print("t =", inst.t)
            print("q =", inst.q)
            print("r =", inst.r)
            incorrect += 1
        else:
            print("OK")

    print("failed =", failed, "incorrect =", incorrect,
          "success rate =", (count - failed - incorrect) / count)


def get_arguments_parser():
    """
    Creates arguments parser with necessary options.
    """
    parser = argparse.ArgumentParser(
        description="""
        The script to check the attack.
        """,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "--count",
        help="Number of tests",
        required=True,
        type=int
    )
    parser.add_argument(
        "--size",
        help="Size of matrices",
        required=True,
        type=int
    )
    parser.add_argument(
        "--d_bound",
        help="Bound for degrees of polynomials",
        required=True,
        type=int
    )
    parser.add_argument(
        "--c_bound",
        help="Upper bound for coefficients",
        required=True,
        type=int
    )
    parser.add_argument(
        "--p_bound",
        help="Upper bound for degree of p",
        required=True,
        type=int
    )
    parser.add_argument(
        "--t_bound",
        help="Upper bound for degree of t",
        required=True,
        type=int
    )

    return parser


if __name__ == "__main__":
    args = get_arguments_parser().parse_args()

    check_attack(args.count, args.size, args.c_bound,
                 args.d_bound, args.p_bound, args.t_bound)
