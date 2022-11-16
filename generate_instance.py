"""
An attack on a key exchange protocol based on max-times and min-times
algebras from [M. I. Durcheva, An application of different dioids in public
key cryptography. In AIP Conference Proceedings, vol. 1631, pp. 336-343,
AIP, 2014].

I. Buchinskiy, M. Kotov, A. Treier, 2022
"""

import random
import tropical_algebra as ta
from matrix_utils import calc_triple_product


def generate_random_matrix(n, l, u):
    """
    Generates a random integer matrix of size n, and the entries are in [l, u].
    """
    return [[random.randint(l, u) for j in range(n)] for i in range(n)]


def generate_random_max_poly(d_bound, l, u, sparse_rate):
    """
    Generates a random polinomial of degree in [1, d_bound] over max-times, and the coefficients are in [l, u].
    sparse_rate defines how many coefficients will be zero.
    """
    d = random.randint(1, d_bound)

    result = [random.randint(l, u) for i in range(d + 1)]

    for i in random.sample(range(1, d + 1), int(d * sparse_rate)):
        result[i] = ta.zero_max_times()

    while result[0] == ta.zero_max_times():
        result[0] = random.randint(l, u)

    return result


def generate_random_min_poly(d_bound, l, u, sparse_rate):
    """
    Generates a random polinomial of degree in [1, d_bound]] over min-times, and the coefficients are in [l, u].
    sparse_rate defines how many coefficients will be zero over min-times.
    """
    d = random.randint(1, d_bound)

    result = [random.randint(l, u) for i in range(d + 1)]

    for i in random.sample(range(1, d + 1), int(d * sparse_rate)):
        result[i] = ta.zero_min_times()

    while result[0] == ta.zero_min_times():
        result[0] = random.randint(l, u)

    return result


class Instance:
    pass


def generate_random_instance(n, u, d):
    """
    Generates a random instance of the protocol.
    """
    while True:
        result = Instance()
        result.M = generate_random_matrix(n, 1, u)
        result.N = generate_random_matrix(n, 1, u)
        result.X = generate_random_matrix(n, 1, u)
        result.p = generate_random_max_poly(d, 1, u, 0.5)
        result.t = generate_random_min_poly(d, 1, u, 0.5)
        result.q = generate_random_max_poly(d, 1, u, 0.5)
        result.r = generate_random_min_poly(d, 1, u, 0.5)
        result.A = calc_triple_product(
            result.M, result.N, result.X, result.p, result.t)
        result.B = calc_triple_product(
            result.M, result.N, result.X, result.q, result.r)
        result.kA = calc_triple_product(
            result.M, result.N, result.B, result.p, result.t)
        result.kB = calc_triple_product(
            result.M, result.N, result.A, result.q, result.r)

        if result.kA == result.kB:
            return result
