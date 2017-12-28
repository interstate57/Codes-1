def notNull(b):
    print(b.degree())
    pass


def euclid(A, B):
    a = A
    b = B
    c1 = 1
    c2 = 0
    d1 = 0
    d2 = 1

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




