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

    @staticmethod
    def one(module):
        return F2Q(Poly(makeF2([1])), module)

    def zero(module):
        return F2Q(Poly(makeF2([0])), module)

    def __eq__(self, other):
        assertSameModule(self, other)
        return self.polynom == other.polynom

    def __mul__(self, other):
        assertSameModule(self, other)
        return F2Q((self.polynom * other.polynom) % self.P, self.P)

    def __truediv__(self, other):
        assertSameModule(self, other)
        print (type(other.polynom), type(self.P))
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
        #return self.polynom.__repr__()
        return self.polynom.__repr__() + " mod " + self.P.__repr__()

    def __repr__(self):
        return self.polynom.__repr__() + " mod " + self.P.__repr__()

