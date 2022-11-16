"""
An attack on a key exchange protocol based on max-times and min-times
algebras from [M. I. Durcheva, An application of different dioids in public
key cryptography. In AIP Conference Proceedings, vol. 1631, pp. 336-343,
AIP, 2014].

I. Buchinskiy, M. Kotov, A. Treier, 2022
"""

import tropical_algebra as ta
import unittest
import attack
import generate_instance
import matrix_utils


class TestAttack(unittest.TestCase):
    def test_instance_from_paper(self):
        i = generate_instance.Instance()
        i.M = [[5, 7, 1], [4, 2, 3], [2, 5, 6]]
        i.N = [[2, 1, 3], [7, 5, 4], [3, 1, 9]]
        i.X = [[5, 2, 8], [6, 7, 4], [3, 1, 5]]
        i.p = [1, 5, 10, 0]
        i.t = [3, 1, ta.INFTY]
        i.q = [1, 5, 0]
        i.r = [10, ta.INFTY, 1, ta.INFTY, ta.INFTY]
        i.A = matrix_utils.calc_triple_product(i.M, i.N, i.X, i.p, i.t)
        i.B = matrix_utils.calc_triple_product(i.M, i.N, i.X, i.q, i.r)
        i.kA = matrix_utils.calc_triple_product(i.M, i.N, i.B, i.p, i.t)
        i.kB = matrix_utils.calc_triple_product(i.M, i.N, i.A, i.q, i.r)

        self.assertEqual(i.kA, i.kB)

        k1 = attack.attack(i.M, i.N, i.X, i.A, i.B, 100, 100)

        self.assertEqual(k1, i.kA)

        K = [[263424, 125440, 263424], [
            188160, 94080, 188160], [311040, 155520, 311040]]

        self.assertEqual(K, i.kA)


if __name__ == "__main__":
    unittest.main()
