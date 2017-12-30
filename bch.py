import math
import numpy
from itertools import combinations

from GF import gen_pow_matrix, min_poly, F2Q, power


def hamming_distance(a, b):
    return sum(map(int, str(bin(a^b))[2:]))


def choose_poly(degree):
    with open("primpoly.txt", "r") as f:
        line = f.readlines()[0]
        numbers = map(int, line.strip().split(', '))
        return list(filter(lambda n: n > 2**degree, numbers))[0]
        # pick one number greater than 2^degree


class BCH:
    def __init__(self, n, t):
        q = math.floor(math.log(n + 1, 2))
        d_c = 2 * t + 1

        # 1. Choose a primitive poly of this degree
        module = choose_poly(q)
        self.pm = gen_pow_matrix(module)

        # 2. Get its roots
        necessary_roots = [power(2, i, self.pm) for i in range(1, d_c)]
        a = min_poly(necessary_roots, self.pm)
        # FIXME repair min_poly; extend roots with zero
        self.R, self.g = a[0], a[1]

        # Matrix
        line = numpy.concatenate(
            [0] * q,
            self.g.coef,
            [0] * (q-1)
        )
        self.generator = np.array([line.rotate(i) for i in range(q)])

    def encode(self, U):
        # Encode with g(x) and return a matrix
        return U * self.generator

    def decode(self, W, method='euclid'):
        # Decode with the given method
        pass

    def dist(self):
        # Return code distance. Exhaustive search.
        lengths = sorted(map(hamming_distance, combinations(self.R, 2)))
        # FIXME without zero there may be only one code word, as in trivial 1->7 code
        return list(lengths)[0]


if __name__ == "__main__":
    #b = BCH(7, 2)
    print (hamming_distance(8, 11))
