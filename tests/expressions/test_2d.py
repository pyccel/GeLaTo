# coding: utf-8

import numpy as np

from sympy.core.containers import Tuple
from sympy import symbols
from sympy import Symbol
from sympy import Lambda
from sympy import IndexedBase

from gelato.expression import glt_symbol
from gelato.expression import gelatize
from gelato.expression import normalize_weak_from
from gelato.expression import initialize_weak_form
from gelato.calculus   import (Dot, Cross, Grad, Curl, Rot, Div)


DIM = 2

# ...
def test_2d_1():
    x,y = symbols('x y')

    u = Symbol('u')
    v = Symbol('v')

    a = Lambda((x,y,v,u), Dot(Grad(u), Grad(v)) + u*v)
    print '> input       := {0}'.format(a)

    expr = gelatize(a, dim=DIM)
    print '> gelatized   := {0}'.format(expr)

    expr = normalize_weak_from(expr)
    print '> normal form := {0}'.format(expr)

    # ... create a glt symbol from a string without evaluation
    #     a discretization is defined as a dictionary
    discretization = {"n_elements": [16, 16], "degrees": [3, 3]}

    expr = glt_symbol(expr, dim=DIM, discretization=discretization, evaluate=False)
    print '> glt symbol  := {0}'.format(expr)
    # ...

    print('')
# ...

# ...
def test_2d_2():
    x,y = symbols('x y')

    u = IndexedBase('u')
    v = IndexedBase('v')

    a = Lambda((x,y,v,u), Rot(u) * Rot(v) + Div(u) * Div(v) + 0.2 * Dot(u, v))
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

    # ... create a glt symbol from a string without evaluation
    #     a discretization is defined as a dictionary
    discretization = {"n_elements": [16, 16], "degrees": [3, 3]}

    expr = glt_symbol(expr, dim=DIM, discretization=discretization, evaluate=False)
    print '> glt symbol  := {0}'.format(expr)
    # ...

    print('')
# ...

# ...
def test_2d_3():
    x,y = symbols('x y')

    u = Symbol('u')
    v = Symbol('v')

    a = Lambda((x,y,v,u), Cross(Curl(u), Curl(v)) + 0.2 * u * v)
    print '> input       := {0}'.format(a)

    expr = gelatize(a, dim=DIM)
    print '> gelatized   := {0}'.format(expr)

    expr, info = initialize_weak_form(expr, dim=DIM)
    print '> temp form   := {0}'.format(expr)

    expr = normalize_weak_from(expr)
    print '> normal form := {0}'.format(expr)

    # ... create a glt symbol from a string without evaluation
    #     a discretization is defined as a dictionary
    discretization = {"n_elements": [16, 16], "degrees": [3, 3]}

    expr = glt_symbol(expr, dim=DIM, discretization=discretization, evaluate=False)
    print '> glt symbol  := {0}'.format(expr)
    # ...

    print('')
# ...

# .....................................................
if __name__ == '__main__':
    test_2d_1()
    test_2d_2()
    test_2d_3()
