# coding: utf-8

from sympy.core import Symbol

from matplotlib import pyplot as plt

from gelato.core import glt_symbol_m, glt_symbol_s, glt_symbol_a, glt_symbol_b
from gelato.core import plot_stiffness_symbols

# ...
def test_glt_symbol_1():
    print('============ test_glt_symbol_1 ==============')

    t = Symbol('t')

    n = 1 ; p = 1
    mt = glt_symbol_m(n, p, t)
    print(mt)

    plot_stiffness_symbols(n=64, degrees=[2,3,4,5])
# ...

# .....................................................
if __name__ == '__main__':
    test_glt_symbol_1()
