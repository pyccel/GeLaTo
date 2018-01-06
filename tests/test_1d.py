# coding: utf-8

import numpy as np

from sympy.core.containers import Tuple
from sympy import symbols
from sympy import Symbol
from sympy import Lambda

from gelato.expression import glt_symbol
from gelato.expression import gelatize
from gelato.expression import normalize_weak_from
from gelato.expression import initialize_weak_form
from gelato.calculus   import (Dot, Cross, Grad, Curl, Rot, Div)
from gelato.calculus   import LinearOperator


DIM = 1

# ...
def test_1d_1():
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
    discretization = {"n_elements": [16], "degrees": [3]}

    expr = glt_symbol(expr, dim=DIM, discretization=discretization, evaluate=False)
    print '> glt symbol  := {0}'.format(expr)
    # ...

    print('')
# ...

# .....................................................
if __name__ == '__main__':

    test_1d_1()
