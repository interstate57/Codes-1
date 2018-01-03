from bch import BCH
from GF import to_rev_poly, set_degree

from time import time, sleep
import numpy as np


n = 15
t = 3


def run_load_test(method):
    all_messages = [list(set_degree(to_rev_poly(dec), n)) for dec in range(2**n-1)]
    bch = BCH(n, t)
    start = time()
    rejects = 0
    for message in all_messages:
        result = bch._decode_one(message, 'euclid')
        if np.isnan(result[0]):
            rejects += 1
    print ("Rejects part:", rejects / (2**n - 1))
    print (time() - start)




run_load_test('euclid')
sleep(5)
run_load_test('pgz')

