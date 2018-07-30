# coding: utf-8

from sympy import Symbol
from sympy.core.containers import Tuple
from sympy import symbols
from sympy import pi, cos, sin
from sympy import srepr
from sympy import I

from symfe.core import dx, dy, dz
from symfe.core import Constant
from symfe.core import Field
from symfe.core import grad, dot, inner, cross, rot, curl, div
from symfe.core import H1Space
from symfe.core import TestFunction
from symfe.core import VectorTestFunction
from symfe.core import BilinearForm

from gelato.core import gelatize
from gelato.core import (Mass,
                         Stiffness,
                         Advection,
                         Bilaplacian)

# ...
def test_gelatize_1d_1():
    print('============ test_gelatize_1d_1 =============')

    V = H1Space('V', ldim=1)

    v = TestFunction(V, name='v')
    u = TestFunction(V, name='u')

    nx = symbols('nx', integer=True)
    px = symbols('px', integer=True)
    tx = symbols('tx')

    c1 = Constant('c1')
    c2 = Constant('c2')
    c3 = Constant('c3')
    c4 = Constant('c4')

    # ...
    expected = Mass(px,tx)/nx
    assert(gelatize(BilinearForm((v,u), u*v)) == expected)
    # ...

    # ...
    expected = nx*Stiffness(px,tx)
    assert(gelatize(BilinearForm((v,u), dot(grad(v), grad(u)))) == expected)
    # ...

    # ...
    expected = I*Advection(px,tx)
    assert(gelatize(BilinearForm((v,u), dx(u) * v)) == expected)
    # ...

    # ...
    expected = c1*Mass(px,tx)/nx + c2*I*Advection(px,tx) - c3*I*Advection(px,tx) + c4*nx*Stiffness(px,tx)
    assert(gelatize(BilinearForm((v,u), c1*v*u + c2*dx(u)*v + c3*dx(v)*u + c4*dx(v)*dx(u))) == expected)
    # ...

#    expr = c1*v*u + c2*dx(u)*v + c3*dx(v)*u + c4*dx(v)*dx(u)
#
#    expr = BilinearForm((v,u), expr)
#    print('> input     >>> {0}'.format(expr))
#    print('> gelatized >>> {0}'.format(gelatize(expr)))
# ...


# .....................................................
if __name__ == '__main__':
    test_gelatize_1d_1()
