# coding: utf-8

from sympy import symbols
from sympy import Tuple
from sympy import Matrix
from sympy import srepr

from gelato.core import (dx, dy, dz)
from gelato.core import grad, dot, inner
from gelato.core import Field


# ...
def test_partial_derivatives():
    print('============ test_partial_derivatives ==============')

    x, y, z = symbols('x y z')
    alpha, beta = symbols('alpha beta')

    F = Field('F')

    u = Field('u')
    v = Field('v')
    w = Field('w')
    uvw = Tuple(u,v,w)

    # ...
    assert(dx(x**2) == 2*x)
    assert(dy(x**2) == 0)
    assert(dz(x**2) == 0)

    assert(dx(x*F) == F + x*dx(F))
    assert(dx(uvw) == Matrix([[dx(u), dx(v), dx(w)]]))
    assert(dx(uvw) + dy(uvw) == Matrix([[dx(u) + dy(u),
                                         dx(v) + dy(v),
                                         dx(w) + dy(w)]]))

    expected = Matrix([[alpha*dx(u) + beta*dy(u),
                        alpha*dx(v) + beta*dy(v),
                        alpha*dx(w) + beta*dy(w)]])
    assert(alpha * dx(uvw) + beta * dy(uvw) == expected)
    # ...

#    expr = alpha * dx(uvw) + beta * dy(uvw)
#    print(expr)

#    print('> ', srepr(expr))
#    print('')
# ...

# .....................................................
if __name__ == '__main__':
    test_partial_derivatives()
