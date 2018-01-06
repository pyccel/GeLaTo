# coding: utf-8

import numpy as np

from sympy.core.containers import Tuple
from sympy import symbols
from sympy import Symbol
from sympy import Lambda
from sympy import IndexedBase

from gelato.expression import gelatize
from gelato.expression import normalize_weak_from
from gelato.expression import initialize_weak_form
from gelato.calculus   import (Dot, Cross, Grad, Curl, Rot, Div)


DIM = 3

# ...
def test_3d_1():
    x,y,z = symbols('x y z')

    u = Symbol('u')
    v = Symbol('v')

    a = Lambda((x,y,z,v,u), Dot(Grad(u), Grad(v)))
    print '> input       := {0}'.format(a)

    expr = gelatize(a, dim=DIM)
    print '> gelatized   := {0}'.format(expr)

    expr, info = initialize_weak_form(expr, dim=DIM)
    print '> temp form   := {0}'.format(expr)

    expr = normalize_weak_from(expr)
    print '> normal form := {0}'.format(expr)

    print('')
# ...

# ...
def test_3d_2():
    x,y,z = symbols('x y z')

    u = IndexedBase('u')
    v = IndexedBase('v')


    a = Lambda((x,y,z,v,u), Div(u) * Div(v) + 0.2 * Dot(u, v))
    print '> input       := {0}'.format(a)

    expr = gelatize(a, dim=DIM)
    print '> gelatized   := {0}'.format(expr)

    expr, info = initialize_weak_form(expr, dim=DIM)
    print '> temp form   :='
    # for a nice printing, we print the dictionary entries one by one
    for key, value in expr.items():
        print '\t\t', key, '\t', value

    expr = normalize_weak_from(expr)
    print '> normal form := {0}'.format(expr)

    print('')
# ...

# ...
def test_3d_3():
    x,y,z = symbols('x y z')

    u = IndexedBase('u')
    v = IndexedBase('v')

    a = Lambda((x,y,z,v,u), Dot(Curl(u), Curl(v)) + 0.2 * Dot(u, v))
    print '> input       := {0}'.format(a)

    expr = gelatize(a, dim=DIM)
    print '> gelatized   := {0}'.format(expr)

    expr, info = initialize_weak_form(expr, dim=DIM)
    print '> temp form   :='
    # for a nice printing, we print the dictionary entries one by one
    for key, value in expr.items():
        print '\t\t', key, '\t', value

    expr = normalize_weak_from(expr)
    print '> normal form := {0}'.format(expr)

    print('')
# ...

# ...
def test_3d_4a():
    x,y,z = symbols('x y z')

    u = IndexedBase('u')
    v = IndexedBase('v')

    b = Tuple(1.0, 0., 0.)

    a = Lambda((x,y,z,v,u), Dot(Curl(Cross(b,u)), Curl(Cross(b,v))) + 0.2 * Dot(u, v))
    print '> input       := {0}'.format(a)

    expr = gelatize(a, dim=DIM)
    print '> gelatized   := {0}'.format(expr)

    expr, info = initialize_weak_form(expr, dim=DIM)
    print '> temp form   :='
    # for a nice printing, we print the dictionary entries one by one
    for key, value in expr.items():
        print '\t\t', key, '\t', value

    expr = normalize_weak_from(expr)
    print '> normal form := {0}'.format(expr)

    print('')
# ...

# ...
def test_3d_4b():
    """Alfven operator."""
    x,y,z = symbols('x y z')

    u = IndexedBase('u')
    v = IndexedBase('v')

    b = IndexedBase('b')

    c0,c1,c2 = symbols('c0 c1 c2')

    a = Lambda((x,y,z,v,u), (  c0 * Dot(u, v)
                             - c1 * Div(u) * Div(v)
                             + c2 *Dot(Curl(Cross(b,u)), Curl(Cross(b,v)))))
    print '> input       := {0}'.format(a)

    expr = gelatize(a, dim=DIM)
    print '> gelatized   := {0}'.format(expr)

    # TODO: fix, not working
#    expr, info = initialize_weak_form(expr, dim=DIM)
#    print '> temp form   :='
#    # for a nice printing, we print the dictionary entries one by one
#    for key, value in expr.items():
#        print '\t\t', key, '\t', value
#
#    expr = normalize_weak_from(expr)
#    print '> normal form := {0}'.format(expr)

    print('')
# ...

# .....................................................
if __name__ == '__main__':
    test_3d_1()
    test_3d_2()
    test_3d_3()
    test_3d_4a()
    test_3d_4b()
