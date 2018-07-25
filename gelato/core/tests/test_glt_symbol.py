# coding: utf-8

from sympy.core import Symbol
from sympy import cos, sin

from matplotlib import pyplot as plt

from gelato.core import glt_symbol_m, glt_symbol_s, glt_symbol_a, glt_symbol_b

# ...
def test_glt_symbol_1():
    print('============ test_glt_symbol_1 ==============')

    t = Symbol('t')
    n = Symbol('n')

    one = 1
    two = 2
    three = 3

    # ... linear splines
    assert( glt_symbol_m(one, t, n) == (cos(t)/3 + 2/3)/n )
    assert( glt_symbol_s(one, t, n) == n*(-2*cos(t) + 2) )
    assert( glt_symbol_a(one, t, n) == -sin(t) )
    assert( glt_symbol_b(one, t, n) == 0 )
    # ...

    # ... quadratic splines
    assert( glt_symbol_m(two, t, n) == (13*cos(t)/30 + cos(2*t)/60 + 11/20)/n )
    assert( glt_symbol_s(two, t, n) == n*(-2*cos(t)/3 - cos(2*t)/3 + 1) )
    assert( glt_symbol_a(two, t, n) == -5*sin(t)/6 - sin(2*t)/12 )
    assert( glt_symbol_b(two, t, n) == n**3*(-8*cos(t) + 2*cos(2*t) + 6) )
    # ...

    # ... cubic splines
    assert( glt_symbol_m(three, t, n) == (397*cos(t)/840 + cos(2*t)/21 + cos(3*t)/2520 + 4793650774/9999999959)/n )
    assert( glt_symbol_s(three, t, n) == n*(-cos(t)/4 - 2*cos(2*t)/5 - cos(3*t)/60 + 2/3))
    assert( glt_symbol_a(three, t, n) == -49*sin(t)/72 - 7*sin(2*t)/45 - sin(3*t)/360 )
    assert( glt_symbol_b(three, t, n) == n**3*(-3*cos(t) + cos(3*t)/3 + 8/3) )
    # ...

#    p = 3
#    print(glt_symbol_m(p, n, t))
#    print(glt_symbol_s(p, n, t))
#    print(glt_symbol_a(p, n, t))
#    print(glt_symbol_b(p, n, t))
# ...

# .....................................................
if __name__ == '__main__':
    test_glt_symbol_1()
