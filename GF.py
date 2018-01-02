from copy import copy
import numpy as np
from numpy.polynomial.polynomial import Polynomial

"""
Pretty good thing to remember: numpy polynomials are [a_0, a_1, ..., a_n] <!>
"""

"""
Realization of F2 class:
    
"""
def value(obj):
    return int(obj)


def euclid_poly(A, B):
    a = A
    b = B
    c1 = Poly(F2(1))
    c2 = Poly(F2(0))
    d1 = Poly(F2(0))
    d2 = Poly(F2(1))

    while b.notNull():
        q = a // b
        r = a - b*q
        c1, c2 = c2, c1 - q * c2
        d1, d2 = d2, d1 - q * d2
        a = b
        b = r
    NOD = a
    #a = c1 * A + d1 * B
    return (NOD, c1, d1)


class F2:
    # consists of : value (0, 1) + operations

    def __init__(self, value):
        self.value = value % 2

    def __neg__(self):
        return F2(-self.value)

    def __int__(self):
        return self.value

    def __repr__(self):
        return str(self.value)

    def __eq__(self, other):
        if value(self) == value(other):
            return True
        return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def __add__(self, other : int):
        return F2((value(self) + value(other)) % 2)

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        return F2((value(self) - value(other)) % 2)

    def __rsub__(self, other):
        return self.__sub__(other)

    def __mul__(self, other):
        return F2((value(self) * value(other)))

    def __rmul__(self, other):
        return __mul__(self, other)

    def __rtruediv__(self, other):
        return self.__truediv(other)

    def __truediv__(self, other):
        if other == 0:
            raise RuntimeError("Don't divide by zero, please")
        return self



class Poly(Polynomial):
    def notNull(self):
        if self.degree() == 0 and self.coef[0] == 0:
            return False
        return True

    def __repr__(self):
        return str(list(self.coef))

    def __str__(self):
        return str(list(self.coef))


def poly_one():
    return Poly(makeF2([1]))


def poly_zero():
    return Poly(makeF2([0]))


def poly_x():
    return Poly(makeF2([0, 1]))



def makeF2(array):
    return [F2(int(x)) for x in array]


def assertSameModule(lhs, rhs):
    if not lhs.P == rhs.P:
        raise RuntimeError("Different modules detected")

class F2Q:
    # ``module'' should be an irreducible F2[x]
    def __init__(self, poly, module):
        self.polynom = poly % module
        self.P = module

    def __hash__(self):
        as_str = '0b' + ''.join(map(str, self.polynom.coef)) + ''.join(map(str, self.P.coef))
        return int(as_str, 2)

    @staticmethod
    def one(module):
        return F2Q(Poly(makeF2([1])), module)

    @staticmethod
    def zero(module):
        return F2Q(Poly(makeF2([0])), module)

    @staticmethod
    def x(module):
        return F2Q(Poly(makeF2([0, 1])), module)

    def __eq__(self, other):
        assertSameModule(self, other)
        return self.polynom == other.polynom

    def __mul__(self, other):
        assertSameModule(self, other)
        return F2Q((self.polynom * other.polynom) % self.P, self.P)

    def __truediv__(self, other):
        assertSameModule(self, other)
        gcd, k, l = euclid_poly(other.polynom, self.P)
        inverse = F2Q(k, self.P)
        return self * inverse

    def __add__(self, other):
        assertSameModule(self, other)
        return F2Q((self.polynom + other.polynom) % self.P, self.P)

    def __sub__(self, other):
        assertSameModule(self, other)
        return F2Q((self.polynom - other.polynom) % self.P, self.P)

    def __str__(self):
        return self.polynom.__repr__() + " mod " + self.P.__repr__()

    def __repr__(self):
        return self.polynom.__repr__() + " mod " + self.P.__repr__()


def to_poly(dec):
    coeff_list = list(map(int, str(bin(dec))[2:]))
    coeff_list.reverse()
    return Poly(makeF2(coeff_list))


def to_rev_poly(dec):
    return list(map(int, str(bin(dec))[2:]))


def from_rev_poly(poly):
    rev = list(poly)
    rev.reverse()
    as_poly = Poly(makeF2(map(int, rev)))
    return to_dec(as_poly)

def to_dec(poly):
    coeff_list = list(poly.coef)
    coeff_list.reverse()
    return int('0b'+''.join(map(str, coeff_list)), 2)


def min_poly(rs, pm):

    def cycle_poly(cycle, pm):
        result = np.array([1])
        for root in cycle:
            result = polyprod(result, np.array([1, root]), pm)
        return result

    cycles = []
    remaining = set(rs)
    while len(remaining) > 0:
        element = list(remaining)[0]
        cycle = [element]
        next = mult(element, element, pm)
        while next != element:
            cycle.append(next)
            next = mult(next, next, pm)
        cycles.append(cycle)
        remaining -= set(cycle)

    # result is the product of all minimal polynomials for each cycle
    result = np.array([1])
    for cycle in cycles:
        result = polyprod(result, cycle_poly(cycle, pm), pm)

    assert all([a in [0, 1] for a in result]), "Internal error: wrong minimal polynomial"
    return result, np.array(sum(cycles, []))


def gen_pow_matrix(primpoly):
    module = to_poly(primpoly)
    column_1 = [0] * (2 ** module.degree())
    column_2 = [0] * (2 ** module.degree())
    i = 1
    x = F2Q(Poly(makeF2([0, 1])), module)
    deg_i = x
    while 1:
        column_2[i] = to_dec(deg_i.polynom)
        column_1[to_dec(deg_i.polynom)] = i
        if deg_i == F2Q.one(module):
            break
        deg_i *= x
        i += 1

    if (i != len(column_1) - 1):
        raise RuntimeError("Not a primitive polynomial. " + str(i))

    return np.array(list(zip(column_1[1:], column_2[1:])))


# Quick f2q operations
def mult(lhs, rhs, pm):
    if lhs == 0 or rhs == 0:
        return 0
    j1 = pm[lhs - 1][0]
    j2 = pm[rhs - 1][0]
    module = len(pm)
    return pm[((j1 + j2) % module - 1) % module][1]

def div(lhs, rhs, pm):
    if rhs == 0:
        raise RuntimeError("F2Q error: division by zero")
    if lhs == 0:
        return 0
    j1 = pm[lhs - 1][0]
    j2 = pm[rhs - 1][0]
    module = len(pm)
    return pm[((j1 - j2) % module - 1) % module][1]

def plus(lhs, rhs):
    return lhs ^ rhs

def minus(lhs, rhs):
    return plus(lhs, rhs)

def power(base, deg, pm):
    result = 1
    for i in range(deg):
        result = mult(result, base, pm)
    return result

# -------- other -------------------------

def shift(ar, step):
    lst = list(ar)
    return np.array(lst[-step:] + lst[:-step])

def set_length(p, l):
    assert l >= len(p), "Can't shorten a polynomial from {0} to {1}".format(len(p), l)
    return np.array([0] * (l - len(p)) + list(p))

def set_degree(p, deg):
    return set_length(p, deg + 1)

def is_zero(p):
    return all([item == 0 for item in p])

# -------- Functions for F2q[x] ----------

def polyval(p, x, pm):
    result = 0
    deg = 1
    for i in range(len(p)):
        result = plus(result, mult(deg, p[-i-1], pm))
        deg = mult(deg, x, pm)
    return result

def strip_leading_zeroes(p):
    while len(p) > 1 and p[0] == 0:
        p = p[1:]
    return p

def make_equilong(p1, p2):
    common_length = max(len(p1), len(p2))
    return map(lambda x: set_length(x, common_length), (p1, p2))

def polyadd(p1, p2):
    pp1, pp2 = make_equilong(p1, p2)
    result = [plus(a, b) for (a, b) in zip(pp1, pp2)]
    result = strip_leading_zeroes(result)
    return np.array(result)

def polysub(p1, p2):
    return polyadd(p1, p2) # F_2^q

def polyprod(p1, p2, pm):
    n = len(p1) - 1
    k = len(p2) - 1
    result = np.zeros(n + k + 1, dtype=int)
    for m in range(n + k + 1):
        i_max = min(m, n)
        j_max = min(m, k)
        i_min = m - j_max
        j_min = m - i_max
        result[m] = 0
        for i in range (i_min, i_max + 1):
            j = m - i
            result[m] = plus(result[m], mult(p1[i], p2[j], pm))
    result = strip_leading_zeroes(result)
    return result

def degree(ar):
    if len(ar) == 1 and ar[0] == 0:
        return -1
    return len(ar) - 1

# Add leading zeroes to a polynomial
def deg_resize(p, deg):
    diff = deg + 1 - len(p)
    return np.array([0] * diff + list(p))

def polydiv(p1, p2, pm):

    a, b = p1, p2

    if is_zero(b):
        raise RuntimeError("Polynomial error: division by zero")
    b = strip_leading_zeroes(b)

    # separate, non-convenient case
    if list(b) == [1]:
        return copy(a), [0]

    q = []
    next_deg = degree(a)
    while degree(a) >= degree(b):
        next_deg -= 1
        q_i = div(a[0], b[0], pm)
        m = degree(a) - degree(b)
        a = polysub(a, polyprod(np.array([q_i] + [0] * m), b, pm))
        a = deg_resize(a, next_deg)
        q += [q_i]

    if not q:
        q = [0]
    return strip_leading_zeroes(q), strip_leading_zeroes(a)

def euclid(A, B, pm, max_deg=0):
    # dirty hack
    if max_deg == 0:
        max_deg = -2
    a, b = A, B
    c1 = np.array([1])
    c2 = np.array([0])
    d1 = np.array([0])
    d2 = np.array([1])

    while not is_zero(b):
        q, r = polydiv(a, b, pm)
        c_new = polysub(c1, polyprod(q, c2, pm))
        d_new = polysub(d1, polyprod(q, d2, pm))
        if degree(r) <= max_deg:
            return r, c_new, d_new
        a, b = b, r
        c1, c2 = c2, c_new
        d1, d2 = d2, d_new

    NOD = a
    # a = c1 * A + d1 * B
    return (NOD, c1, d1)
