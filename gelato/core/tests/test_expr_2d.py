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
from gelato.core import Glt
from gelato.core import (Mass,
                         Stiffness,
                         Advection,
                         Bilaplacian)


# ...
def test_gelatize_2d_1():
    print('============ test_gelatize_2d_1 =============')

    V = H1Space('V', ldim=2)

    v = TestFunction(V, name='v')
    u = TestFunction(V, name='u')

    nx, ny = symbols('nx ny', integer=True)
    px, py = symbols('px py', integer=True)
    tx, ty = symbols('tx ty')


    c = Constant('c')

    bx = Constant('bx')
    by = Constant('by')
    b = Tuple(bx, by)

    # ...
    expected = Mass(px,tx)*Mass(py,ty)/(nx*ny)
    assert(gelatize(BilinearForm((v,u), u*v), evaluate=False) == expected)
    # ...

    # ...
    expected = Mass(py,ty)*nx*Stiffness(px,tx)/ny
    assert(gelatize(BilinearForm((v,u), dx(u)*dx(v)), evaluate=False) == expected)
    # ...

    # ...
    expected = I*Advection(py,ty)*Mass(px,tx)/nx
    assert(gelatize(BilinearForm((v,u), dy(u) * v), evaluate=False) == expected)
    # ...

    # ...
    expected = I*Advection(px,tx)*Mass(py,ty)/ny
    assert(gelatize(BilinearForm((v,u), dx(u) * v), evaluate=False) == expected)
    # ...

    # ...
    expected = Mass(px,tx)*ny*Stiffness(py,ty)/nx + Mass(py,ty)*nx*Stiffness(px,tx)/ny
    assert(gelatize(BilinearForm((v,u), dot(grad(v), grad(u))), evaluate=False) == expected)
    # ...

    # ...
    expected = (nx*Mass(py,ty)*Stiffness(px,tx)/ny +
                I*Advection(px,tx)*Mass(py,ty)/ny +
                ny*Mass(px,tx)*Stiffness(py,ty)/nx +
                I*Advection(py,ty)*Mass(px,tx)/nx)
    assert(gelatize(BilinearForm((v,u), dot(grad(v), grad(u)) + dx(u)*v + dy(u)*v), evaluate=False) == expected)
    # ...

    # ...
    expected = -bx*I*Advection(px,tx)*Mass(py,ty)/ny - by*I*Advection(py,ty)*Mass(px,tx)/nx
    assert(gelatize(BilinearForm((v,u), dot(b, grad(v)) * u), evaluate=False) == expected)
    # ...

    # ...
    expected = bx**2*nx*Mass(py,ty)*Stiffness(px,tx)/ny + 2*bx*by*Advection(px,tx)*Advection(py,ty) + by**2*ny*Mass(px,tx)*Stiffness(py,ty)/nx
    assert(gelatize(BilinearForm((v,u), dot(b, grad(v)) * dot(b, grad(u))), evaluate=False) == expected)
    # ...

    degrees = None
#    degrees = [2, 1]

#    evaluate = True
    evaluate = False

#    expr = dot(b, grad(v)) * dot(b, grad(u))
#
#    expr = BilinearForm((v,u), expr)
#    print('> input     >>> {0}'.format(expr))
#    print('> gelatized >>> {0}'.format(gelatize(expr, degrees, evaluate=evaluate)))
# ...

# ...
def test_gelatize_2d_2():
    print('============ test_gelatize_2d_2 =============')

    V = H1Space('V', ldim=2)

    v = TestFunction(V, name='v')
    u = TestFunction(V, name='u')

    mass = BilinearForm((v,u), u*v)
    laplace = BilinearForm((v,u), dot(grad(v), grad(u)))

    # ...
    degrees = [1,1]
#    degrees = None
    symbol_m = Glt(mass, degrees=degrees, n_elements=[4,4])
    print(symbol_m)
    # ...

    # ...
    degrees = [2,2]
    degrees = None
    symbol_l = Glt(laplace, degrees=degrees)
    print(symbol_l)
    # ...

#    from sympy import simplify
#    print(symbol_l/symbol_m)

# ...

# .....................................................
if __name__ == '__main__':
    test_gelatize_2d_1()
    test_gelatize_2d_2()
