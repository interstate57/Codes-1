import unittest

import numpy as np
from numpy.polynomial.polynomial import *

from GF import Poly, F2, F2Q, euclid, value, to_poly, to_dec, gen_pow_matrix#, min_poly

data = [
    [1, 1, 0, 0, 1],    # x^4 + x + 1
    [1, 0, 1],          # x^2 + 1
    [1, 1],             # x + 1
    [0, 0, 1],          # x^2
    [1],                # 1
    [0],
    [1, 1, 0, 1],
    [0, 1],
    [1, 0, 0, 1]
]

# FIXME: use toF2
data_poly = list(map(
    lambda x: Poly(x),
    map(
        lambda ar: [F2(x) for x in ar],
        data
    )
))

a,b,c,d,e,f,g,h,i = data_poly


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
        module = a
        one = F2Q(e, module)
        zero = F2Q(f, module)
        x = F2Q(h, module)
        another_x = F2Q(h, module)
        x3plus1 = F2Q(i, module)
        self.assertEqual(x, another_x)
        self.assertEqual(x * x3plus1, one)
        self.assertEqual(x3plus1 * x, one)
        self.assertEqual(one / x3plus1, x)
        self.assertEqual(one / x, x3plus1)
        self.assertEqual(x + x, x - x)
        self.assertEqual(x3plus1 + x3plus1, zero)

    def test_pow_matrix(self):
        self.assertEqual(to_poly(19), a)
        self.assertEqual(to_dec(a), 19)
        result = gen_pow_matrix(11)
        print (result)

    def test_min_poly(self):
        primitive = F2Q(h, g)
        pass
        #min_poly([primitive, primitive * primitive, primitive*primitive*primitive])
        #aa = Poly([F2Q.one(g)], domain=F2Q)
        #bb = Poly([F2Q.zero(g), F2Q.one(g)], domain=F2Q)
        #cc = bb * aa


if __name__ == "__main__":
    unittest.main()
