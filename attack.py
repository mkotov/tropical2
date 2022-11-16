"""
An attack on a key exchange protocol based on max-times and min-times
algebras from [M. I. Durcheva, An application of different dioids in public
key cryptography. In AIP Conference Proceedings, vol. 1631, pp. 336-343,
AIP, 2014].

I. Buchinskiy, M. Kotov, A. Treier, 2022
"""

import tropical_algebra as ta
from matrix_utils import calc_min, calc_max, calc_triple_product


def find_t_coeff(A, B):
    """
    Returns max(A / B) if this number is an integer, infty otherwise.
    This function doesn't use float-point arithmetic, so instead of comparing a/b vs c/d, we compare ad vs cb.
    """
    n = len(A)
    a = A[0][0]
    b = B[0][0]
    for i in range(n):
        for j in range(n):
            if A[i][j] * b >= B[i][j] * a:
                a = A[i][j]
                b = B[i][j]
    if a % b == 0:
        return a // b
    return ta.INFTY


def is_matrix_div_matrix_const(A, B):
    """
    If A / B is a const, then returns this const. Returns None otherwise.
    This function doesn't use float-point arithmetic, so instead of comparing a/b vs c/d, we compare ad vs cb and return a pair of numbers.
    """
    n = len(A)
    ra = None
    rb = None
    ok = False
    # Find a pair of non-zero finite elements of matrices, to compute a const.
    for i in range(n):
        for j in range(n):
            if A[i][j] != 0 and A[i][j] != ta.INFTY and B[i][j] != 0 and B[i][j] != ta.INFTY:
                ra = A[i][j]
                rb = B[i][j]
                ok = True
                break
        if ok:
            break

    if not ra or not rb:
        return None

    # Skip pairs 0 and 0, infty and infty.
    for i in range(n):
        for j in range(n):
            if A[i][j] == 0 and B[i][j] == 0:
                continue
            if A[i][j] == ta.INFTY and B[i][j] == ta.INFTY:
                continue
            if A[i][j] == ta.INFTY or B[i][j] == ta.INFTY:
                return None

            if A[i][j] * rb != B[i][j] * ra:
                return None

    return ra, rb


def is_matrix_repeated(As):
    """
    Returns True iff the last matrix is const * As[i] for some i < len(As) - 1.
    """
    for i in range(len(As) - 1):
        if is_matrix_div_matrix_const(As[-1], As[i]):
            return True

    return False


def find_polys(n, M, N, X, A, p_bound, t_bound):
    """
    Given the public matrices M, N, X, Alice's matrix A. Returns p' and t'.
    """
    minA = calc_min(A)
    maxA = calc_max(A)

    p = [1]

    Mi = [ta.one_matrix_max_times(n)]
    Nj = [ta.one_matrix_min_times(n)]
    Npd = None

    for i in range(p_bound + 1):
        if is_matrix_repeated(Mi):
            return None, None

        MiX = ta.mul_matrices_max_times(Mi[-1], X)
        minMiX = calc_min(MiX)

        if minMiX > minA:
            return None, None

        t = []
        tN = ta.zero_matrix_min_times(n)

        for j in range(t_bound + 1):
            if Npd:
                if j == Npd:
                    break
            else:
                if is_matrix_repeated(Nj):
                    Npd = j
                    break

            if calc_min(ta.mul_matrices_min_times(MiX, Nj[j])) > maxA:
                break

            t.insert(0, find_t_coeff(A, ta.mul_matrices_min_times(MiX, Nj[j])))
            tN = ta.sum_matrices_min_times(
                tN, ta.mul_matrix_by_coef_min_times(Nj[j], t[0]))

            if ta.mul_matrices_min_times(MiX, tN) == A:
                return p, t

            if len(Nj) == j + 1 and j < t_bound:
                Nj.append(ta.mul_matrices_min_times(Nj[-1], N))

        p.append(0)

        if i < p_bound:
            Mi.append(ta.mul_matrices_max_times(Mi[-1], M))

    return None, None


def attack(M, N, X, A, B, p_bound, t_bound):
    """
    The implementation of our attack on the protocol.
    """
    n = len(M)
    p1, t1 = find_polys(n, M, N, X, A, p_bound, t_bound)
    if p1 is None or t1 is None:
        return None

    q1, r1 = find_polys(n, M, N, X, B, p_bound, t_bound)
    if q1 is None or r1 is None:
        return None

    k1 = calc_triple_product(M, N, B, p1, t1)
    k2 = calc_triple_product(M, N, A, q1, r1)

    if k1 == k2:
        return k1

    return None
