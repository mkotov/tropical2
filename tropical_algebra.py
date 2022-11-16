"""
An attack on a key exchange protocol based on max-times and min-times
algebras from [M. I. Durcheva, An application of different dioids in public
key cryptography. In AIP Conference Proceedings, vol. 1631, pp. 336-343,
AIP, 2014].

I. Buchinskiy, M. Kotov, A. Treier, 2022
"""

import sys

INFTY = "infty"
"""This constant represent +infinity."""


def zero_max_times():
    """
    Returns the zero element of R_max-times.
    """
    return 0


def one_max_times():
    """
    Returns the unit element of R_max-times.
    """
    return 1


def sum_max_times(a, b):
    """
    Returns the sum of two elements in R_max-times.
    """
    if a == INFTY:
        return INFTY
    if b == INFTY:
        return INFTY
    return max(a, b)


def mul_max_times(a, b):
    """
    Returns the product of two elements in R_max-times.
    """
    if a == INFTY:
        if b == 0:
            return 0
        return INFTY
    if b == INFTY:
        if a == 0:
            return 0
        return INFTY
    return a * b


def zero_min_times():
    """
    Returns the zero element of R_min-times.
    """
    return INFTY


def one_min_times():
    """
    Returns the unit element of R_min-times.
    """
    return 1


def sum_min_times(a, b):
    """
    Returns the sum of two elements in R_min-times.
    """
    if a == INFTY:
        return b
    if b == INFTY:
        return a
    return min(a, b)


def mul_min_times(a, b):
    """
    Returns the product of two elements in R_min-times.
    """
    if a == INFTY:
        return INFTY
    if b == INFTY:
        return INFTY
    return a * b


def sum_matrices_semiring(A, B, sum_elements):
    """
    Returns the sum of two matrices over a semiring.
    """
    n = len(A)
    return [[sum_elements(A[i][j], B[i][j]) for j in range(n)] for i in range(n)]


def sum_matrices_max_times(A, B):
    """
    Returns the sum of two matrices over R_max-times.
    """
    return sum_matrices_semiring(A, B, sum_max_times)


def sum_matrices_min_times(A, B):
    """
    Returns the sum of two matrices over R_min-times.
    """
    return sum_matrices_semiring(A, B, sum_min_times)


def zero_matrix_semiring(n, zero_element):
    """
    Returns the zero matrix of size n over a semiring.
    """
    return [[zero_element() for row in range(n)] for col in range(n)]


def zero_matrix_max_times(n):
    """
    Returns the zero matrix of size n over R_max-times.
    """
    return zero_matrix_semiring(n, zero_max_times)


def zero_matrix_min_times(n):
    """
    Returns the zero matrix of size n over R_min-times.
    """
    return zero_matrix_semiring(n, zero_min_times)


def mul_matrices_semiring(A, B, sum_elements, mul_elements, zero_element):
    """
    Returns the product of two matrices over a semiring.
    """
    n = len(A)
    C = zero_matrix_semiring(n, zero_element)

    for i in range(n):
        for j in range(n):
            for k in range(n):
                C[i][j] = sum_elements(C[i][j], mul_elements(A[i][k], B[k][j]))
    return C


def mul_matrices_max_times(A, B):
    """
    Returns the product of two matrices over R_max-times.
    """
    return mul_matrices_semiring(A, B, sum_max_times, mul_max_times, zero_max_times)


def mul_matrices_min_times(A, B):
    """
    Returns the product of two matrices over R_min-times.
    """
    return mul_matrices_semiring(A, B, sum_min_times, mul_min_times, zero_min_times)


def mul_matrix_by_coef_semiring(A, coef, mul_elements):
    """
    Returns the product of an element of a semiring and a matrix over the semiring.
    """
    n = len(A)
    return [[mul_elements(A[i][j], coef) for j in range(n)] for i in range(n)]


def mul_matrix_by_coef_max_times(A, coef):
    """
    Returns the product of an element of R_max-times and a matrix over this structure.
    """
    return mul_matrix_by_coef_semiring(A, coef, mul_max_times)


def mul_matrix_by_coef_min_times(A, coef):
    """
    Returns the product of an element of R_min-times and a matrix over this structure.
    """
    return mul_matrix_by_coef_semiring(A, coef, mul_min_times)


def one_matrix_semiring(n, zero_element, one_element):
    """
    Returns the unit matrix of size n over a semiring.
    """
    return [[zero_element() if row != col else one_element() for col in range(n)] for row in range(n)]


def one_matrix_max_times(n):
    """
    Returns the unit matrix of size n over R_max_times.
    """
    return one_matrix_semiring(n, zero_max_times, one_max_times)


def one_matrix_min_times(n):
    """
    Returns the unit matrix of size n over R_min_times.
    """
    return one_matrix_semiring(n, zero_min_times, one_min_times)


def pwr_matrix_semiring(A, m, sum_elements, mul_elements, zero_element, one_element):
    """
    Returns a matrix raised to the power m over a semiring.
    """
    n = len(A)
    if m == 0:
        return one_matrix_semiring(n, zero_element, one_element)
    if m % 2 == 0:
        return pwr_matrix_semiring(mul_matrices_semiring(A, A, sum_elements, mul_elements, zero_element), m // 2,
                                   sum_elements, mul_elements, zero_element, one_element)
    else:
        return mul_matrices_semiring(A, pwr_matrix_semiring(A, m - 1, sum_elements, mul_elements, zero_element,
                                                            one_element), sum_elements, mul_elements, zero_element)


def pwr_matrix_max_times(A, m):
    """
    Returns a matrix raised to the power m over R_max-times.
    """
    return pwr_matrix_semiring(A, m, sum_max_times, mul_max_times, zero_max_times, one_max_times)


def pwr_matrix_min_times(A, m):
    """
    Returns a matrix raised to the power m over R_min-times.
    """
    return pwr_matrix_semiring(A, m, sum_min_times, mul_min_times, zero_min_times, one_min_times)


def calc_poly_matrix_semiring(A, p, sum_elements, mul_elements, zero_element, one_element):
    """
    Given a matrix A and a polynomial p over a semiring. Returns p(A).
    """
    n = len(A)
    d = len(p) - 1
    C = zero_matrix_semiring(n, zero_element)
    D = one_matrix_semiring(n, zero_element, one_element)
    for i in range(d + 1):
        C = sum_matrices_semiring(C, mul_matrix_by_coef_semiring(
            D, p[d - i], mul_elements), sum_elements)
        if i != d:
            D = mul_matrices_semiring(
                D, A, sum_elements, mul_elements, zero_element)

    return C


def calc_poly_matrix_max_times(A, p):
    """
    Given a matrix A and a polynomial p over R_max-times. Returns p(A).
    """
    return calc_poly_matrix_semiring(A, p, sum_max_times, mul_max_times, zero_max_times, one_max_times)


def calc_poly_matrix_min_times(A, p):
    """
    Given a matrix A and a polynomial p over R_min-times. Returns p(A).
    """
    return calc_poly_matrix_semiring(A, p, sum_min_times, mul_min_times, zero_min_times, one_min_times)
