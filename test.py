import unittest

import numpy as np
from copy import copy
from itertools import combinations

from GF import plus, minus, mult, div, power
from GF import polyprod, min_poly, shift, polyadd, polysub, polydiv, polyval, set_degree, euclid, is_zero
from GF import Poly, F2, F2Q, euclid_poly, value, to_poly, to_dec, gen_pow_matrix

from solve import linsolve

from bch import BCH, hamming_distance, choose_poly


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



class TestBase(unittest.TestCase):
    def assertEqualArrays(self, first, second):
        self.assertEqual(list(first), list(second))

    def assertEqualTables(self, first, second):
        self.assertEqual(len(first), len(second))
        for (a, b) in zip(first, second):
            self.assertEqualArrays(a, b)


class TestGF(TestBase):

    def test_euclid_poly(self):
        gcd, k, l = euclid_poly(b, c)
        self.assertEqual(gcd, c)
        self.assertEqual(k, f)
        self.assertEqual(l, e)

        gcd, k, l = euclid_poly(a, b)
        self.assertEqual(gcd, e)
        self.assertEqual(k, h)
        self.assertEqual(l, g)

    def test_euclid(self):
        # with fake max_deg
        pm = gen_pow_matrix(11)
        gcd, k, l = euclid(np.array([1, 0, 1]), np.array([1, 1]), pm)
        self.assertEqualArrays(gcd, [1, 1])
        self.assertEqualArrays(k, [0])
        self.assertEqualArrays(l, [1])

        gcd, k, l = euclid(np.array([1, 0, 0, 1, 1]), np.array([1, 0, 1]), pm)
        self.assertEqualArrays(gcd, [1])
        self.assertEqualArrays(k, [1, 0])
        self.assertEqualArrays(l, [1, 0, 1, 1])

        # with active max_deg
        pm = gen_pow_matrix(19)
        gcd, k, l = euclid(np.array([1, 0, 0, 0, 0, 0, 0, 0]), np.array([2, 1, 3, 5, 4, 2, 1]), pm, 3)
        self.assertEqualArrays(gcd, [6, 0, 13])
        self.assertEqualArrays(l, [2, 6, 9, 13])

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

    # FIXME improve
    def test_quick_f2q(self):
        pm = gen_pow_matrix(11)
        self.assertEqual(div(0, 4, pm), 0)

    def test_pow_matrix(self):
        self.assertEqual(to_poly(19), a)
        self.assertEqual(to_dec(a), 19)
        self.assertEqualTables(
            gen_pow_matrix(11),
            [
                [7, 2],
                [1, 4],
                [3, 3],
                [2, 6],
                [6, 7],
                [4, 5],
                [5, 1]
            ]
        )

    def test_polyval(self):
        pm = gen_pow_matrix(11)
        self.assertEqual(polyval([1, 0, 1, 1], 1, pm), 1)
        self.assertEqual(polyval([1, 0, 1, 1], 2, pm), 0)
        self.assertEqual(polyval([1, 0, 1, 1], 0, pm), 1)

    def test_polyprod(self):
        pm = gen_pow_matrix(11)
        self.assertEqualArrays(
            polyprod([0, 0, 0], [0, 0], pm),
            [0]
        )

    def test_poly_operations(self):
        self.assertEqualArrays(polyadd([1, 0, 0], [1]), [1, 0, 1])
        a = np.array([3, 2, 1])
        self.assertEqualArrays(set_degree(a, 4), np.array([0, 0, 3, 2, 1]))
        b = np.array([1])
        pm = gen_pow_matrix(11)
        self.assertEqualArrays(polyprod(a, b, pm), a)
        b = np.array([2, 1])
        self.assertEqualArrays(polyprod(a, b, pm), [6, 7, 0, 1])
        b = np.array([0])
        self.assertEqualArrays(polyprod(b, a, pm), [0])
        self.assertEqualArrays(polyprod(a, b, pm), [0])
        c = np.array([3, 2, 1])
        d = np.array([3, 0, 1])
        self.assertEqualArrays(polyadd(c, d), [2, 0])
        self.assertEqualArrays(polysub(c, d), [2, 0])
        e = [1, 0, 0, 0, 1]
        f = [1, 0]
        q, r = polydiv(e, f, pm)
        self.assertEqualArrays(q, [1, 0, 0, 0])
        self.assertEqualArrays(r, [1])
        q, r = polydiv(f, e, pm)
        self.assertEqualArrays(q, [0])
        self.assertEqualArrays(r, f)
        g = [0, 1]
        q, r = polydiv(f, g, pm)
        self.assertEqualArrays(q, f)
        self.assertEqualArrays(r, [0])


    def test_min_poly(self):
        # for x^3 + x + 1
        pm = gen_pow_matrix(11)
        m_p, roots = min_poly(np.array([2, 3, 4]), pm)
        self.assertEqualArrays(m_p, [1, 1, 1, 1, 1, 1, 1])
        self.assertEqual(set(roots), set([2, 3, 4, 5, 6, 7]))
        # for x^4 + x + 1
        pm = gen_pow_matrix(19)
        m_p, roots = min_poly(np.array([2, 4, 8, 3]), pm)
        self.assertEqualArrays(m_p, [1, 1, 1, 0, 1, 0, 0, 0, 1])
        self.assertEqual(set(roots), set([2, 3, 4, 5, 8, 10, 12, 15]))

    def test_other(self):
        ar = np.array([1, 2, 3, 0, 0])
        self.assertEqualArrays(shift(ar, 1), [0, 1, 2, 3, 0])
        self.assertEqualArrays(shift(ar, 3), [3, 0, 0, 1, 2])
        self.assertEqualArrays(shift(ar, 5), [1, 2, 3, 0, 0])
        self.assertEqual(is_zero([0, 0, 0]), True)
        self.assertEqual(is_zero([0, 4]), False)


class TestBCH(TestBase):
    def test_hamming_distance(self):
        self.assertEqual(hamming_distance(57, 9), 2)

    def test_init(self):
        bch = BCH(7, 1)
        self.assertEqualArrays(bch.g, [1, 0, 1, 1])
        self.assertEqual(set(bch.R), set([2, 4, 6]))
        self.assertEqual(bch.dist(), 3)

        bch = BCH(7, 2)
        self.assertEqualArrays(bch.g, [1, 1, 1, 1, 1, 1, 1])
        self.assertEqual(set(bch.R), set([2, 4, 6, 3, 5, 7]))
        self.assertEqual(bch.dist(), 7)

        # cases of extra distance: 31 & 6; 63 & 14; 127 & 30
        #bch = BCH(127, 30)
        #print(bch.dist())

    def test_encode_1(self):
        bch = BCH(7, 1)
        U = np.array([[1, 1, 0, 1], [1, 0, 0, 0]])
        result = bch.encode(U)
        self.assertEqualTables(result, [[1, 1, 0, 1, 0, 0, 1], [1, 0, 0, 0, 1, 0, 1]])
        # Check that the code is systematic
        for msg in bch._all_messages():
            encoded = bch._encode_one(msg)
            self.assertEqualArrays(msg, encoded[:bch.k])

    def test_encode_2(self):
        bch = BCH(15, 3)
        msg = [0, 1, 1, 1, 1]
        self.assertEqualArrays(bch._encode_one(msg), [0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0])


    def check_decode_codewords(self, n, t):
        bch = BCH(n, t)
        # Check that the code is systematic
        for msg in bch._all_messages():
            encoded = bch._encode_one(msg)
            decoded = bch._decode_one(encoded)
            self.assertEqualArrays(msg, decoded)

    def check_decode_with_errors(self, n, t, method):
        bch = BCH(15, 3)

        #msg = [0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1] # no mistakes
        #msg.reverse()
        #msg = [1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1] # 3 mistakes
        for message in bch._all_codewords():
            result = bch._decode_one(message, method=method)
            self.assertEqualArrays(message[:bch.k], result)
            for i in [1, 2, 3]:
                for comb in combinations(range(bch.n), i):
                    new_message = copy(message)
                    for k in comb:
                        new_message[k] = 1 - new_message[k]
                    result = bch._decode_one(new_message, method=method)
                    self.assertEqualArrays(result, message[:bch.k])


    def test_compare_methods_internal(self):
        bch = BCH(15, 3)
        s = [13, 7, 3, 9, 4, 2, 1]
        by_euclid, er1 = bch._euclid(s)
        by_pgz, er2 = bch._pgz(s)
        self.assertEqual(er1, False)
        self.assertEqual(er2, False)
        self.assertEqualArrays(by_euclid, by_pgz)

    def test_compare_methods(self):
        bch = BCH(15, 3)
        word = [0, 0, 1, 1, 1]
        encoded = bch._encode_one(word)
        # 1 mistake
        encoded[2] = 1 - encoded[2]
        self.assertEqualArrays(word, bch._decode_one(encoded, method='euclid'))
        self.assertEqualArrays(word, bch._decode_one(encoded, method='pgz'))

        # 2 mistakes
        encoded[5] = 1 - encoded[5]
        self.assertEqualArrays(word, bch._decode_one(encoded, method='euclid'))
        self.assertEqualArrays(word, bch._decode_one(encoded, method='pgz'))

        # 3 mistakes
        encoded[8] = 1 - encoded[8]
        self.assertEqualArrays(word, bch._decode_one(encoded, method='euclid'))
        self.assertEqualArrays(word, bch._decode_one(encoded, method='pgz'))

    def _test_compare_methods_big(self):
        bch = BCH(15, 3)
        msg = [0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1]
        #for message in bch._all_codewords():
        for message in [msg]:
            result = bch._decode_one(message)
            self.assertEqualArrays(message[:bch.k], result)
            print ("Pizda")
            for i in [1]:
                for comb in combinations(range(bch.n), i):
                    new_message = copy(message)
                    for k in comb:
                        new_message[k] = 1 - new_message[k]
                    result = bch._decode_one(new_message, 'euclid')
                    self.assertEqualArrays(result, message[:bch.k])



    def test_decode(self):
        print ("Caution: long test")
        self.check_decode_codewords(7, 1)
        self.check_decode_codewords(7, 2)
        self.check_decode_codewords(15, 2)
        self.check_decode_codewords(15, 3)
        self.check_decode_with_errors(7, 1, 'euclid')
        self.check_decode_with_errors(15, 3, 'pgz')


if __name__ == "__main__":
    unittest.main()
