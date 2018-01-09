# coding: utf-8

import numpy as np

from sympy.core.containers import Tuple
from sympy import symbols
from sympy import Symbol
from sympy import Lambda
from sympy import IndexedBase

from gelato.expression import construct_weak_form
from gelato.calculus   import (Dot, Cross, Grad, Curl, Rot, Div)
from gelato.calculus   import Constant


DIM = 2

# ...
def test_2d_1():
    x,y = symbols('x y')

    u = Symbol('u')
    v = Symbol('v')

    a = Lambda((x,y,v,u), Dot(Grad(u), Grad(v)) + u*v)
    print '> input       := {0}'.format(a)

    # ...
    expr = construct_weak_form(a, dim=DIM)
    print '> weak form := {0}'.format(expr)
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

    # ...
    expr = construct_weak_form(a, dim=DIM, is_block=True, verbose=True)
    print '> weak form := {0}'.format(expr)
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

    # ...
    expr = construct_weak_form(a, dim=DIM, is_block=True)
    print '> weak form := {0}'.format(expr)
    # ...

    print('')
# ...

# ...
def test_2d_4():
    x,y = symbols('x y')

    u = Symbol('u')
    v = Symbol('v')

    bx = Constant('bx')
    by = Constant('by')
    b = Tuple(bx, by)

    a = Lambda((x,y,v,u), 0.2 * u * v + Dot(b, Grad(v)) * u)
    print '> input       := {0}'.format(a)

    # ...
    expr = construct_weak_form(a, dim=DIM, is_block=False)
    print '> weak form := {0}'.format(expr)
    # ...

    print('')
# ...

# .....................................................
if __name__ == '__main__':
    test_2d_1()
    test_2d_2()
    test_2d_3()
    test_2d_4()
