"""
An attack on a key exchange protocol based on max-times and min-times
algebras from [M. I. Durcheva, An application of different dioids in public
key cryptography. In AIP Conference Proceedings, vol. 1631, pp. 336-343,
AIP, 2014].

I. Buchinskiy, M. Kotov, A. Treier, 2022
"""

import sys
import unittest
import tropical_algebra
import random


class TestTropicalAlgebra(unittest.TestCase):
    def test_sum_max_times(self):
        for i in range(100):
            for j in range(100):
                self.assertEqual(
                    max(i, j), tropical_algebra.sum_max_times(i, j))

    def test_mul_max_times(self):
        for i in range(100):
            for j in range(100):
                self.assertEqual(i * j, tropical_algebra.mul_max_times(i, j))

    def test_sum_matrices_max_times(self):
        A = [[1, 2, 3],
             [4, 5, 6],
             [7, 8, 9]]
        B = [[2, 1, 4],
             [3, 6, 5],
             [8, 7, 10]]
        C = [[2, 2, 4],
             [4, 6, 6],
             [8, 8, 10]]
        self.assertEqual(C, tropical_algebra.sum_matrices_max_times(A, B))

    def test_mul_matrices_max_times(self):
        A = [[1, 2, 3],
             [4, 5, 6],
             [7, 8, 9]]
        B = [[2, 1, 4],
             [3, 6, 5],
             [8, 7, 10]]
        C = [[24, 21, 30],
             [48, 42, 60],
             [72, 63, 90]]
        self.assertEqual(C, tropical_algebra.mul_matrices_max_times(A, B))

    def test_mul_matrix_by_coef_max_times(self):
        for i in range(20):
            for j in range(i):
                for coef in range(i):
                    A = [[i + j for j in range(10)] for i in range(10)]
                    expected = [
                        [A[i][j] * coef for j in range(10)] for i in range(10)]
                    self.assertEqual(
                        expected, tropical_algebra.mul_matrix_by_coef_max_times(A, coef))

    def test_pwr_matrix_max_times_1(self):
        A = [[1, 2, 3],
             [4, 5, 6],
             [7, 8, 9]]
        pwr = 3
        expected = [[189, 216, 243],
                    [378, 432, 486],
                    [567, 648, 729]]
        self.assertEqual(
            expected, tropical_algebra.pwr_matrix_max_times(A, pwr))

    def test_pwr_matrix_max_times_2(self):
        A = [[1, 2, 3],
             [4, 5, 6],
             [7, 8, 9]]
        pwr = 4
        expected = [[1701, 1944, 2187],
                    [3402, 3888, 4374],
                    [5103, 5832, 6561]]
        self.assertEqual(
            expected, tropical_algebra.pwr_matrix_max_times(A, pwr))

    def test_calc_poly_matrix_max_times(self):
        A = [[5, 7, 1], [4, 2, 3], [2, 5, 6]]
        poly = [1, 5, 10, 0]
        expected = [[140, 196, 126], [112, 140, 108], [120, 180, 216]]
        self.assertEqual(
            expected, tropical_algebra.calc_poly_matrix_max_times(A, poly))

    def test_sum_min_times(self):
        for i in range(100):
            for j in range(100):
                self.assertEqual(
                    min(i, j), tropical_algebra.sum_min_times(i, j))

    def test_mul_min_times(self):
        for i in range(100):
            for j in range(100):
                self.assertEqual(i * j, tropical_algebra.mul_min_times(i, j))

    def test_sum_matrices_min_times(self):
        A = [[1, 2, 3],
             [4, 5, 6],
             [7, 8, 9]]
        B = [[2, 1, 4],
             [3, 6, 5],
             [8, 7, 10]]
        C = [[1, 1, 3],
             [3, 5, 5],
             [7, 7, 9]]
        self.assertEqual(C, tropical_algebra.sum_matrices_min_times(A, B))

    def test_mul_matrices_min_times(self):
        A = [[1, 2, 3],
             [4, 5, 6],
             [7, 8, 9]]
        B = [[2, 1, 4],
             [3, 6, 5],
             [8, 7, 10]]
        C = [[2, 1, 4],
             [8, 4, 16],
             [14, 7, 28]]
        self.assertEqual(C, tropical_algebra.mul_matrices_min_times(A, B))

    def test_mul_matrix_by_coef_min_times(self):
        for i in range(20):
            for j in range(i):
                for coef in range(i):
                    A = [[i + j for j in range(10)] for i in range(10)]
                    expected = [
                        [A[i][j] * coef for j in range(10)] for i in range(10)]
                    self.assertEqual(
                        expected, tropical_algebra.mul_matrix_by_coef_min_times(A, coef))

    def test_pwr_matrix_min_times_1(self):
        A = [[1, 2, 3],
             [4, 5, 6],
             [7, 8, 9]]
        pwr = 2
        C = [[1, 2, 3],
             [4, 8, 12],
             [7, 14, 21]]
        self.assertEqual(C, tropical_algebra.pwr_matrix_min_times(A, pwr))

    def test_pwr_matrix_min_times_2(self):
        A = [[3, 25, 37],
             [44, 52, 64],
             [71, 83, 95]]
        pwr = 4
        C = [[81, 675, 999],
             [1188, 9900, 14652],
             [1917, 15975, 23643]]
        self.assertEqual(C, tropical_algebra.pwr_matrix_min_times(A, pwr))

    def test_calc_poly_matrix_min_times(self):
        A = [[1, 2, 3],
             [4, 5, 6],
             [7, 8, 9]]
        poly = [2, 3, 4, 5]
        expected = [[2, 4, 6],
                    [8, 5, 24],
                    [14, 28, 5]]
        self.assertEqual(
            expected, tropical_algebra.calc_poly_matrix_min_times(A, poly))

    def test_from_paper_1(self):
        # see M. I. Durcheva, An application of different dioids in public key cryptography, 2014.
        M = [[5, 7, 1], [4, 2, 3], [2, 5, 6]]
        N = [[2, 1, 3], [7, 5, 4], [3, 1, 9]]
        X = [[5, 2, 8], [6, 7, 4], [3, 1, 5]]
        p = [1, 5, 10, 0]
        t = [3, 1, tropical_algebra.INFTY]
        q = [1, 5, 0]
        r = [10, tropical_algebra.INFTY, 1,
             tropical_algebra.INFTY, tropical_algebra.INFTY]
        pM = [[140, 196, 126], [112, 140, 108], [120, 180, 216]]
        tN = [[2, 1, 3], [7, 5, 4], [3, 1, 9]]
        A = [[2352, 1120, 3528], [1680, 840, 2520], [2160, 1080, 3240]]
        qM = [[28, 35, 21], [20, 28, 18], [20, 30, 36]]
        rN = [[4, 2, 4], [12, 4, 20], [6, 3, 4]]
        B = [[840, 420, 840], [672, 336, 640], [720, 360, 720]]
        K = [[263424, 125440, 263424], [
            188160, 94080, 188160], [311040, 155520, 311040]]
        self.assertEqual(pM, tropical_algebra.calc_poly_matrix_max_times(M, p))
        self.assertEqual(qM, tropical_algebra.calc_poly_matrix_max_times(M, q))
        self.assertEqual(tN, tropical_algebra.calc_poly_matrix_min_times(N, t))
        self.assertEqual(rN, tropical_algebra.calc_poly_matrix_min_times(N, r))
        self.assertEqual(A, tropical_algebra.mul_matrices_min_times(
            tropical_algebra.mul_matrices_max_times(pM, X), tN))
        self.assertEqual(B, tropical_algebra.mul_matrices_min_times(
            tropical_algebra.mul_matrices_max_times(qM, X), rN))
        self.assertEqual(K, tropical_algebra.mul_matrices_min_times(
            tropical_algebra.mul_matrices_max_times(pM, B), tN))
        self.assertEqual(K, tropical_algebra.mul_matrices_min_times(
            tropical_algebra.mul_matrices_max_times(qM, A), rN))

    def test_zero_and_one_1(self):
        a = 2
        self.assertEqual(a, tropical_algebra.mul_max_times(
            a, tropical_algebra.one_max_times()))
        self.assertEqual(tropical_algebra.zero_max_times(),
                         tropical_algebra.mul_max_times(a, tropical_algebra.zero_max_times()))
        self.assertEqual(a, tropical_algebra.mul_max_times(
            tropical_algebra.one_max_times(), a))
        self.assertEqual(tropical_algebra.zero_max_times(),
                         tropical_algebra.mul_max_times(tropical_algebra.zero_max_times(), a))
        self.assertEqual(a, tropical_algebra.mul_min_times(
            a, tropical_algebra.one_min_times()))
        self.assertEqual(tropical_algebra.zero_min_times(),
                         tropical_algebra.mul_min_times(a, tropical_algebra.zero_min_times()))
        self.assertEqual(a, tropical_algebra.mul_min_times(
            tropical_algebra.one_min_times(), a))
        self.assertEqual(tropical_algebra.zero_min_times(),
                         tropical_algebra.mul_min_times(tropical_algebra.zero_min_times(), a))
        b = tropical_algebra.INFTY
        self.assertEqual(b, tropical_algebra.mul_max_times(
            b, tropical_algebra.one_max_times()))
        self.assertEqual(tropical_algebra.zero_max_times(),
                         tropical_algebra.mul_max_times(b, tropical_algebra.zero_max_times()))
        self.assertEqual(b, tropical_algebra.mul_max_times(
            tropical_algebra.one_max_times(), b))
        self.assertEqual(tropical_algebra.zero_max_times(),
                         tropical_algebra.mul_max_times(tropical_algebra.zero_max_times(), b))
        self.assertEqual(b, tropical_algebra.mul_min_times(
            b, tropical_algebra.one_min_times()))
        self.assertEqual(tropical_algebra.zero_min_times(),
                         tropical_algebra.mul_min_times(b, tropical_algebra.zero_min_times()))
        self.assertEqual(b, tropical_algebra.mul_min_times(
            tropical_algebra.one_min_times(), b))
        self.assertEqual(tropical_algebra.zero_min_times(),
                         tropical_algebra.mul_min_times(tropical_algebra.zero_min_times(), b))
        c = 0
        self.assertEqual(c, tropical_algebra.mul_max_times(
            c, tropical_algebra.one_max_times()))
        self.assertEqual(tropical_algebra.zero_max_times(),
                         tropical_algebra.mul_max_times(c, tropical_algebra.zero_max_times()))
        self.assertEqual(c, tropical_algebra.mul_max_times(
            tropical_algebra.one_max_times(), c))
        self.assertEqual(tropical_algebra.zero_max_times(),
                         tropical_algebra.mul_max_times(tropical_algebra.zero_max_times(), c))
        self.assertEqual(c, tropical_algebra.mul_min_times(
            c, tropical_algebra.one_min_times()))
        self.assertEqual(tropical_algebra.zero_min_times(),
                         tropical_algebra.mul_min_times(c, tropical_algebra.zero_min_times()))
        self.assertEqual(c, tropical_algebra.mul_min_times(
            tropical_algebra.one_min_times(), c))
        self.assertEqual(tropical_algebra.zero_min_times(),
                         tropical_algebra.mul_min_times(tropical_algebra.zero_min_times(), c))

    def test_powers_1(self):
        A = [[1, 2], [3, 4]]
        self.assertEqual(tropical_algebra.one_matrix_max_times(
            2), tropical_algebra.pwr_matrix_max_times(A, 0))
        self.assertEqual(A, tropical_algebra.pwr_matrix_max_times(A, 1))
        self.assertEqual(tropical_algebra.mul_matrices_max_times(
            A, A), tropical_algebra.pwr_matrix_max_times(A, 2))
        self.assertEqual(tropical_algebra.one_matrix_min_times(
            2), tropical_algebra.pwr_matrix_min_times(A, 0))
        self.assertEqual(A, tropical_algebra.pwr_matrix_min_times(A, 1))
        self.assertEqual(tropical_algebra.mul_matrices_min_times(
            A, A), tropical_algebra.pwr_matrix_min_times(A, 2))

    def test_distributivity(self):
        a = 0
        b = tropical_algebra.INFTY
        c = 1
        self.assertEqual(tropical_algebra.mul_max_times(a, tropical_algebra.sum_max_times(b, c)),
                         tropical_algebra.sum_max_times(tropical_algebra.mul_max_times(a, b), tropical_algebra.mul_max_times(a, c)))
        self.assertEqual(tropical_algebra.mul_max_times(b, tropical_algebra.sum_max_times(a, c)),
                         tropical_algebra.sum_max_times(tropical_algebra.mul_max_times(b, a), tropical_algebra.mul_max_times(b, c)))
        self.assertEqual(tropical_algebra.mul_max_times(b, tropical_algebra.sum_max_times(a, a)),
                         tropical_algebra.sum_max_times(tropical_algebra.mul_max_times(b, a), tropical_algebra.mul_max_times(b, a)))
        self.assertEqual(tropical_algebra.mul_max_times(a, tropical_algebra.sum_max_times(b, b)),
                         tropical_algebra.sum_max_times(tropical_algebra.mul_max_times(a, b), tropical_algebra.mul_max_times(a, b)))
        self.assertEqual(tropical_algebra.mul_min_times(a, tropical_algebra.sum_min_times(b, c)),
                         tropical_algebra.sum_min_times(tropical_algebra.mul_min_times(a, b), tropical_algebra.mul_min_times(a, c)))
        self.assertEqual(tropical_algebra.mul_min_times(b, tropical_algebra.sum_min_times(a, c)),
                         tropical_algebra.sum_min_times(tropical_algebra.mul_min_times(b, a), tropical_algebra.mul_min_times(b, c)))
        self.assertEqual(tropical_algebra.mul_min_times(b, tropical_algebra.sum_min_times(a, a)),
                         tropical_algebra.sum_min_times(tropical_algebra.mul_min_times(b, a), tropical_algebra.mul_min_times(b, a)))
        self.assertEqual(tropical_algebra.mul_min_times(a, tropical_algebra.sum_min_times(b, b)),
                         tropical_algebra.sum_min_times(tropical_algebra.mul_min_times(a, b), tropical_algebra.mul_min_times(a, b)))


if __name__ == "__main__":
    unittest.main()
