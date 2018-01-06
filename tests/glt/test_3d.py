# coding: utf-8

import numpy as np

from sympy.core.containers import Tuple
from sympy import symbols
from sympy import Symbol
from sympy import Lambda
from sympy import IndexedBase

from gelato.glt import glt_symbol
from gelato.calculus   import (Dot, Cross, Grad, Curl, Rot, Div)


DIM = 3

# ...
def test_3d_1():
    x,y,z = symbols('x y z')

    u = Symbol('u')
    v = Symbol('v')

    a = Lambda((x,y,z,v,u), Dot(Grad(u), Grad(v)))
    print '> input       := {0}'.format(a)

    # ... create a glt symbol from a string without evaluation
    #     a discretization is defined as a dictionary
    discretization = {"n_elements": [16, 16, 16], "degrees": [3, 3, 3]}

    expr = glt_symbol(a, dim=DIM, discretization=discretization, evaluate=False)
    print '> glt symbol  := {0}'.format(expr)
    # ...

    print('')
# ...

# ...
def test_3d_2():
    x,y,z = symbols('x y z')

    u = IndexedBase('u')
    v = IndexedBase('v')


    a = Lambda((x,y,z,v,u), Div(u) * Div(v) + 0.2 * Dot(u, v))
    print '> input       := {0}'.format(a)

    # ... create a glt symbol from a string without evaluation
    #     a discretization is defined as a dictionary
    discretization = {"n_elements": [16, 16, 16], "degrees": [3, 3, 3]}

    expr = glt_symbol(a, dim=DIM, discretization=discretization, evaluate=False, is_block=True)
    print '> glt symbol  := {0}'.format(expr)
    # ...

    print('')
# ...

# ...
def test_3d_3():
    x,y,z = symbols('x y z')

    u = IndexedBase('u')
    v = IndexedBase('v')

    a = Lambda((x,y,z,v,u), Dot(Curl(u), Curl(v)) + 0.2 * Dot(u, v))
    print '> input       := {0}'.format(a)

    # ... create a glt symbol from a string without evaluation
    #     a discretization is defined as a dictionary
    discretization = {"n_elements": [16, 16, 16], "degrees": [3, 3, 3]}

    expr = glt_symbol(a, dim=DIM, discretization=discretization, evaluate=False, is_block=True)
    print '> glt symbol  := {0}'.format(expr)
    # ...

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

    # ... create a glt symbol from a string without evaluation
    #     a discretization is defined as a dictionary
    discretization = {"n_elements": [16, 16, 16], "degrees": [3, 3, 3]}

    expr = glt_symbol(a, dim=DIM, discretization=discretization, evaluate=False, is_block=True)
    print '> glt symbol  := {0}'.format(expr)
    # ...

    print('')
# ...

## ... TODO fix it
#def test_3d_4b():
#    """Alfven operator."""
#    x,y,z = symbols('x y z')
#
#    u = IndexedBase('u')
#    v = IndexedBase('v')
#
#    b = IndexedBase('b')
#
#    c0,c1,c2 = symbols('c0 c1 c2')
#
#    a = Lambda((x,y,z,v,u), (  c0 * Dot(u, v)
#                             - c1 * Div(u) * Div(v)
#                             + c2 *Dot(Curl(Cross(b,u)), Curl(Cross(b,v)))))
#    print '> input       := {0}'.format(a)
#
#    # ... create a glt symbol from a string without evaluation
#    #     a discretization is defined as a dictionary
#    discretization = {"n_elements": [16, 16, 16], "degrees": [3, 3, 3]}
#
#    expr = glt_symbol(a, dim=DIM, discretization=discretization, evaluate=False, is_block=True)
#    print '> glt symbol  := {0}'.format(expr)
#    # ...
#
#    print('')
## ...

# .....................................................
if __name__ == '__main__':
    test_3d_1()
    test_3d_2()
    test_3d_3()
    test_3d_4a()
#    test_3d_4b()
