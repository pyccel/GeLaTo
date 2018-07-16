# coding: utf-8

from sympy import symbols

from gelato.core import (dx, dy, dz)
from gelato.core import grad, dot, inner


# ...
def test_partial_derivatives():
    print('============ test_partial_derivatives ==============')

    x, y, z = symbols('x y z')
    u = symbols('u')

    # ...
    expr = dx(x**2)
    print(expr)
    # ...

    print('')
# ...

# .....................................................
if __name__ == '__main__':
    test_partial_derivatives()
