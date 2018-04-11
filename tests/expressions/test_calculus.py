# coding: utf-8

import numpy as np

from sympy.core.containers import Tuple
from sympy import symbols
from sympy import Symbol
from sympy import Lambda

from gelato.calculus import (dx, dy, dz)
from gelato.calculus import LinearOperator


# ...
def test_0():
    x,y, a = symbols('x y a')

    # ...
    expr = x+y
    print('> expr := {0}'.format(expr))

    expr = LinearOperator(expr)
    print('> gelatized := {0}'.format(expr))
    print('')
    # ...

    # ...
    expr = 2*x+y
    print('> expr := {0}'.format(expr))

    expr = LinearOperator(expr)
    print('> gelatized := {0}'.format(expr))
    print('')
    # ...

    # ...
    expr = a*x+y
    print('> expr := {0}'.format(expr))

    expr = LinearOperator(expr)
    print('> gelatized := {0}'.format(expr))
    # ...

    # ...
    expr = 2*a*x+y
    print('> expr := {0}'.format(expr))

    expr = LinearOperator(expr)
    print('> gelatized := {0}'.format(expr))
    # ...
# ...

# ...
def test_1():
    u, v, a = symbols('u v a')

    # ...
    expr = u+v
    print('> expr := {0}'.format(expr))

    expr = dx(expr)
    print('> gelatized := {0}'.format(expr))
    print('')
    # ...

    # ...
    expr = 2*u*v
    print('> expr := {0}'.format(expr))

    expr = dx(expr)
    print('> gelatized := {0}'.format(expr))
    print('')
    # ...

    # ... dx should not operate on u^2,
    #     since we consider only linearized weak formulations
    expr = u*u
    print('> expr := {0}'.format(expr))

    expr = dx(expr)
    print('> gelatized := {0}'.format(expr))
    print('')
    # ...
# ...

# .....................................................
if __name__ == '__main__':
    test_0()
    test_1()
