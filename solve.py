import numpy as np
from functools import reduce

from GF import plus, minus, mult, div, gen_pow_matrix


# --- row operations ---
def row_plus(a, b):
    assert len(a) == len(b), "Can't add vectors of different length"
    return list(map(lambda x: plus(x[0], x[1]), zip(a, b)))

def row_minus(a, b):
    return row_plus(a, b)

def row_mult(a, k, pm):
    return [mult(c, k, pm) for c in a]

def scalar_mult(a, b, pm):
    assert len(a) == len(b), "Can't multiply vectors of different length"
    result = 0
    for x, y in zip(a, b):
        result = plus(result, mult(x, y, pm))
    return result

# ----------------------

def choose_nonzero_row(A, _k):
    k = _k
    while k < A.shape[0] and A[k, _k] == 0:
        k += 1
    if (k >= A.shape[0]):
        return np.nan
    return k

# Gauss
def linsolve(A, b, pm):
    assert A.shape[0] == A.shape[1], "Not a square matrix. No gauss."
    assert A.shape[0] == b.shape[0], "Invalid system of equations"
    A = np.hstack((A, np.array([b]).T))

    n = A.shape[0]
    for k in range(n):
        # Maybe A[k, k] = 0. Swap line k with some better underlying one.
        nonzero_row = choose_nonzero_row(A, k)
        if np.isnan(nonzero_row): # no suitable rows. FIXME
            return np.array([np.nan] * b.shape[0])
        if nonzero_row != k:
            A[k, :], A[nonzero_row, :] = A[nonzero_row, :], np.copy(A[k, :])

        # Divide row by diagonal element
        if A[k, k] != 1:
            inverse = div(1, A[k, k], pm)
            A[k] = row_mult(A[k],inverse, pm)

        # Zeroize lower part of column k.
        for row in range(k + 1, n):
            A[row] = row_minus(A[row], row_mult(A[k], A[row, k], pm))

    x = np.array([0 for i in range(n)])
    for k in range(n - 1, -1, -1):
        x[k] = div(
            minus(
                A[k, -1],
                scalar_mult(
                    A[k, k:n],
                    x[k:n],
                    pm
                )
            ),
            A[k, k],
            pm
        )
    return x
