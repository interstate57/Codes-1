import unittest

import numpy as np
from numpy.polynomial.polynomial import *

from GF import Poly, F2, euclid, value

data = [
    [1, 1, 0, 0, 1],    # x^4 + x + 1
    [1, 0, 1],          # x^2 + 1
    [1, 1],             # x + 1
    [0, 0, 1],          # x^2
    [1],                # 1
    [0],
    [1, 1, 0, 1],
    [0, 1]
]

data_poly = list(map(
    lambda x: Poly(x),
    map(
        lambda ar: [F2(x) for x in ar],
        data
    )
))

a,b,c,d,e,f,g,h = data_poly


class TestEverything(unittest.TestCase):
    def test_euclid(self):
        gcd, k, l = euclid(b, c)
        self.assertEqual(gcd, c)
        self.assertEqual(k, f)
        self.assertEqual(l, e)

        gcd, k, l = euclid(a, b)
        self.assertEqual(gcd, e)
        self.assertEqual(k, h)
        self.assertEqual(l, g)
        #print ("k:", k, "\n", "l:", l)

    def test_f2(self):
        zero = F2(0)
        one = F2(1)
        self.assertEqual(value(0), 0)
        self.assertEqual(value(1), 1)
        self.assertEqual(value(F2(1)), 1)
        self.assertEqual(value(F2(0)), 0)
        self.assertNotEqual(F2(1), 0)
        self.assertEqual(one * one, one)
        self.assertEqual(one - one, zero)
        self.assertEqual(one / one * one + zero * one - one * one * one, zero)
        self.assertEqual(zero / one, zero)
        with self.assertRaises(RuntimeError):
            one / zero

    def test_f2q(self):
        pass


if __name__ == "__main__":
    unittest.main()
