# coding: utf-8

from sympy import symbols
from sympy import srepr

from gelato.core import (dx, dy, dz)
from gelato.core import grad, dot, inner
from gelato.core import Field


# ...
def test_partial_derivatives():
    print('============ test_partial_derivatives ==============')

    x, y, z = symbols('x y z')
    u = symbols('u')
    F = Field('F')

    # ...
    assert(dx(x**2) == 2*x)
    assert(dy(x**2) == 0)
    assert(dz(x**2) == 0)

    assert(dx(x*F) == F + x*dx(F))
    # ...

    expr = dx(x*F)
#    print(expr)
#    print('> ', srepr(expr))

    print('')
# ...

# .....................................................
if __name__ == '__main__':
    test_partial_derivatives()
