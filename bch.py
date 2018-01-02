import math
import numpy as np
from copy import copy
from itertools import combinations

from GF import gen_pow_matrix, min_poly, F2Q, power, shift, Poly, makeF2, to_rev_poly, to_dec, set_degree, from_rev_poly, degree
from GF import polyval, euclid, is_zero
from solve import linsolve


def hamming_distance(a, b):
    return sum(map(int, str(bin(a^b))[2:]))


def choose_poly(degree):
    with open("primpoly.txt", "r") as f:
        line = f.readlines()[0]
        numbers = map(int, line.strip().split(', '))
        # pick one number greater than 2^degree
        return list(filter(lambda n: n > 2**degree, numbers))[0]


class BCH:
    def __init__(self, n, t):
        if 2 * t + 1 >= n:
            raise RuntimeError("Code of length {0} can't correct as many mistakes as {1}".format(n, t))
        self.n = n
        self.t = t
        self.q = math.floor(math.log(n + 1, 2))
        d_c = 2 * t + 1

        # 1. Choose a primitive poly of this degree
        module = choose_poly(self.q)
        self.pm = gen_pow_matrix(module)

        # 2. Get its roots
        necessary_roots = [power(2, i, self.pm) for i in range(1, d_c)]
        a = min_poly(necessary_roots, self.pm)
        self.g, self.R = a[0], a[1]
        self.k = self.n - degree(self.g)

    def encode(self, U):
        # Encode with g(x) and return a matrix
        return np.array([self._encode_one(u) for u in U])

    def decode(self, W, method='euclid'):
        return np.array([self._decode_one(w, method) for w in W])
        # Decode with the given method
        pass

    def dist(self):
        # Return code distance. Exhaustive search.
        lengths = sorted(
            map(
                lambda x: hamming_distance(x[0], x[1]),
                combinations(
                    map(from_rev_poly, self._all_codewords()),
                    2
                )
            )
        )
        return list(lengths)[0]

    def _all_messages(self):
        return [list(set_degree(to_rev_poly(dec), self.k - 1)) for dec in range(2**self.k)]

    def _all_codewords(self):
        return self.encode(self._all_messages())

    def _encode_one(self, u):
        # v(x) = x^m u(x) + (x^m u(x) % g(x))
        def as_poly(coeff):
            poly_coeff = list(coeff)
            poly_coeff.reverse()
            return Poly(makeF2(poly_coeff))

        g_poly, u_poly = map(as_poly, (self.g, u))
        m = g_poly.degree()
        x_m = as_poly([1] + [0] * m)
        v_poly = x_m * u_poly + (x_m * u_poly % g_poly)
        v = list(v_poly.coef)
        v.reverse()
        return set_degree(np.array(v), self.n - 1)

    # Calculates error locator by syndrome poly. Euclidean way.
    def _euclid(self, s):
        a = np.array([1] + [0] * (2*self.t + 1))
        b = s
        nod, k, l = euclid(a, b, self.pm, self.t)
        # result: sigma(x) polynomial. False positive cases will be checked outside, in _decode
        return l, False

    def _decode_one(self, _w, method='euclid'):
        w = copy(_w)
        roots = [power(2, i, self.pm) for i in range(1, 2*self.t + 1)]
        roots.reverse()
        # syndrome:
        s = [polyval(w, alpha, self.pm) for alpha in roots]
        # We'l correct mistakes if there are
        if not is_zero(s):
            s += [1]
            # find error locators
            if method == 'euclid':
                error_locator, error_happened = self._euclid(s)
            elif method == 'pgz':
                error_locator, error_happened = self._pgz(s)
            else:
                raise LogicError("Unknown decode method")
            if error_happened:
                return np.array([np.nan for i in range(self.n)])

            # find roots
            locator_roots = list(filter(lambda x: polyval(error_locator, power(2, x, self.pm), self.pm) == 0, range(1, 2**self.q)))
            if len(locator_roots) != degree(error_locator):
                return np.array([np.nan for i in range(self.n)])
            # correct the word
            for root in locator_roots:
                w[root - 1] = 1 - w[root - 1]

        w = list(w)
        # the code is systematic and the message is in k beginning bits, so
        return w[:self.k]

    # fix interface
    def _pgz(self, _s):
        # here syndrome comes in descending order, from max to min degree
        s = copy(list(_s))
        s.reverse()
        for r in range(self.t, 0, -1):
            # build matrix
            A = np.array([s[i:i+r] for i in range(1, r+1)])
            b = np.array(s[r+1:2*r+1])
            result = linsolve(A, b, self.pm)
            if not np.isnan(result[0]):
                return np.array(list(result) + [1]), False
        return None, True
