from copy import copy
import numpy as np
from numpy.polynomial.polynomial import *

"""
Pretty good thing to remember: numpy polynomials are [a_0, a_1, ..., a_n] <!>
"""

"""
Realization of F2 class:
    
"""
def value(obj):
    return int(obj)

def euclid(A, B):
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
    return [F2(x) for x in array]


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
        gcd, k, l = euclid(other.polynom, self.P)
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

def to_dec(poly):
    coeff_list = list(poly.coef)
    coeff_list.reverse()
    return int('0b'+''.join(map(str, coeff_list)), 2)

"""
def min_poly(rs):
    def convert(coeff):
        as_poly_coeff = coeff.polynom.coef
        if (len(as_poly_coef) != 1):
            raise RuntimeError("Internal error: length " + str(len(as_poly_coef)))
        return as_poly_coeff[0]

    def min_poly(cycle):
        module = cycle[0].P

        if (cycle[-1].polynom == Poly(makeF2([0, 1]))):
            return module

        #result = F2Q.one(module)
        #result = poly_one()
        #F2Q_zero = F2Q.zero(module)
        one = Poly([F2Q.one(module)])
        x = Poly([F2Q.zero(module), F2Q.one(module)])

        result = Poly([F2Q.one(module)])
        for element in cycle:
            #result = result * (F2Q.x(module) - element).polynom
            result = result * (x + Poly([element]))

        as_F2_poly = Poly(map(convert, result.coef))
        return as_F2_poly
    
    cycles = []
    remaining = set(rs)
    while len(remaining) > 0:
        element = list(remaining)[0]
        cycle = [element]
        next = element * element
        while next != element:
            cycle.append(next)
            next = next * next
        cycles.append(cycle)
        remaining -= set(cycle)
    # now we have cycles
    # primitive poly for the first cycle is definitely P
    result = Poly(makeF2([1]))
    for cycle in cycles:
        result *= min_poly(cycle)
    return result
"""

# maybe viet's theorem is more convenient than repairing polynomial
def min_poly(rs, pm):
    # FIXME
    return 0, 1


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
    j1 = pm[lhs - 1][0]
    j2 = pm[rhs - 1][0]
    module = len(pm)
    return pm[((j1 + j2) % module - 1) % module][1]

def div(lhs, rhs, pm):
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

