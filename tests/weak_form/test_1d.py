# coding: utf-8

import numpy as np

from sympy.core.containers import Tuple
from sympy import symbols
from sympy import Symbol
from sympy import Lambda

from gelato.expression import construct_weak_form
from gelato.calculus   import (Dot, Cross, Grad, Curl, Rot, Div, dx)
from gelato.calculus   import Constant


DIM = 1

# ...
def test_1d_1():
    x,y = symbols('x y')

    u = Symbol('u')
    v = Symbol('v')

    a = Lambda((x,y,v,u), Dot(Grad(u), Grad(v)) + u*v)
    print('> input     := {0}'.format(a))

    # ...
    expr = construct_weak_form(a, dim=DIM)
    print('> weak form := {0}'.format(expr))
    # ...

    print('')
# ...

# ...
def test_1d_2():
    x,y = symbols('x y')

    u = Symbol('u')
    v = Symbol('v')

#    b = Function('b')
    b = Constant('b')

    a = Lambda((x,y,v,u), Dot(Grad(b*u), Grad(v)) + u*v)
    print('> input       := {0}'.format(a))

    # ...
    expr = construct_weak_form(a, dim=DIM)
    print('> weak form := {0}'.format(expr))
    # ...

    print('')
# ...

# ...
def test_1d_3():
    x = symbols('x')

    u0, u1 = symbols('u0 u1')
    v0, v1 = symbols('v0 v1')

    a = Lambda((x,v0,v1,u0,u1), dx(u0)*dx(v0) + dx(u1)*v0 + u0*dx(v1) + u1*v1)
    print('> input       := {0}'.format(a))

    # ...
    expr = construct_weak_form(a, dim=DIM, is_block=True)
    print('> weak form := {0}'.format(expr))
    # ...

    print('')
# ...

# .....................................................
if __name__ == '__main__':

    test_1d_1()
    test_1d_2()
    test_1d_3()
