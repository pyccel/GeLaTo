# coding: utf-8

from sympy.core import Symbol
from sympy import cos, sin, Rational as frac

from gelato import Mass, Stiffness, Advection, Bilaplacian

#==============================================================================
def test_glt_symbol_1():

    t = Symbol('t')

    one = 1
    two = 2
    three = 3

    # ... linear splines
    assert( Mass(one, t) == cos(t)/3 + frac(2, 3) )
    assert( Stiffness(one, t) == -2*cos(t) + 2 )
    assert( Advection(one, t) == -sin(t) )
#    assert( Bilaplacian(one, t) == 0 )
    # ...

    # ... quadratic splines
    assert( Mass(two, t) == 13*cos(t)/30 + cos(2*t)/60 + frac(11, 20) )
    assert( Stiffness(two, t) == -2*cos(t)/3 - cos(2*t)/3 + 1 )
    assert( Advection(two, t) == -5*sin(t)/6 - sin(2*t)/12 )
#    assert( Bilaplacian(two, t) == -8*cos(t) + 2*cos(2*t) + 6 )
    # ...

    # ... cubic splines
    assert( Mass(three, t) == 397*cos(t)/840 + cos(2*t)/21 + cos(3*t)/2520 + frac(151, 315) )
    assert( Stiffness(three, t) == -cos(t)/4 - 2*cos(2*t)/5 - cos(3*t)/60 + frac(2, 3))
    assert( Advection(three, t) == -49*sin(t)/72 - 7*sin(2*t)/45 - sin(3*t)/360 )
#    assert( Bilaplacian(three, t) == -3*cos(t) + cos(3*t)/3 + 8/3 )
    # ...

#    p = 3
#    print(Mass(p, t))
#    print(Stiffness(p, t))
#    print(Advection(p, t))
##    print(Bilaplacian(p, t))

#==============================================================================
#def test_glt_symbol_2():
#
#    from sympy import limit, sin, Symbol, oo
#    from sympy.abc import x
#
#    t = Symbol('t')
#    p = Symbol('p')
#
#    l = limit(Mass(p, t), p, oo)
#    print(l)

#==============================================================================
# CLEAN UP SYMPY NAMESPACE
#==============================================================================

def teardown_module():
    from sympy import cache
    cache.clear_cache()

def teardown_function():
    from sympy import cache
    cache.clear_cache()
