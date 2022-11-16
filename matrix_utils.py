"""
An attack on a key exchange protocol based on max-times and min-times
algebras from [M. I. Durcheva, An application of different dioids in public
key cryptography. In AIP Conference Proceedings, vol. 1631, pp. 336-343,
AIP, 2014].

I. Buchinskiy, M. Kotov, A. Treier, 2022
"""

import tropical_algebra as ta


def calc_min(A):
    """
    Computes the minimum of the elements of the matrix.
    """
    n = len(A)
    m = A[0][0]
    for i in range(n):
        for j in range(n):
            if m == 0:
                return m
            if m == ta.INFTY:
                m = A[i][j]
            if A[i][j] == ta.INFTY:
                continue
            if A[i][j] < m:
                m = A[i][j]
    return m


def calc_max(A):
    """
    Computes the maximum of the elements of the matrix.
    """
    n = len(A)
    m = A[0][0]
    for i in range(n):
        for j in range(n):
            if m == ta.INFTY or A[i][j] == ta.INFTY:
                return ta.INFTY
            if A[i][j] > m:
                m = A[i][j]
    return m


def calc_triple_product(M, N, X, p, t):
    """
    Returns (p(M) boxtimes X) otimes t(N).
    """
    return ta.mul_matrices_min_times(
        ta.mul_matrices_max_times(
            ta.calc_poly_matrix_max_times(M, p),
            X),
        ta.calc_poly_matrix_min_times(N, t))
