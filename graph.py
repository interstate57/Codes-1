from bch import BCH

import matplotlib.pyplot as plt


def points(n):
    x_values = range(1, n // 2)
    y_values = []
    for t in x_values:
        bch = BCH(n, t)
        y_values.append(bch.k / bch.n)
    return x_values, y_values


if __name__ == "__main__":

    x1, y1 = points(7)
    x2, y2 = points(15)
    x3, y3 = points(31)
    x4, y4 = points(63)
    x5, y5 = points(127)

    plt.plot(x1, y1, 'r^:', x2, y2, 'bD:', x3, y3, 'go:', x4, y4, 'y^:', x5, y5, 'rD:')
    plt.show()
