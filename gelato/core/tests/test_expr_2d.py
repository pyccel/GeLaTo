# coding: utf-8

from sympy import Symbol
from sympy.core.containers import Tuple
from sympy import symbols
from sympy import pi, cos, sin
from sympy import srepr
from sympy import I

from sympde.core import dx, dy, dz
from sympde.core import Domain
from sympde.core import Constant
from sympde.core import Field
from sympde.core import grad, dot, inner, cross, rot, curl, div
from sympde.core import FunctionSpace
from sympde.core import VectorFunctionSpace
from sympde.core import TestFunction
from sympde.core import VectorTestFunction
from sympde.core import BilinearForm

from gelato.core import gelatize
from gelato.core import (Mass,
                         Stiffness,
                         Advection,
                         Bilaplacian)

DIM = 2
domain = Domain('Omega', dim=DIM)

# ...
def test_gelatize_2d_1():
    print('============ test_gelatize_2d_1 =============')

    V = FunctionSpace('V', domain)

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
    assert(gelatize(BilinearForm((v,u), u*v)) == expected)
    # ...

    # ...
    expected = Mass(py,ty)*nx*Stiffness(px,tx)/ny
    assert(gelatize(BilinearForm((v,u), dx(u)*dx(v))) == expected)
    # ...

    # ...
    expected = I*Advection(py,ty)*Mass(px,tx)/nx
    assert(gelatize(BilinearForm((v,u), dy(u) * v)) == expected)
    # ...

    # ...
    expected = I*Advection(px,tx)*Mass(py,ty)/ny
    assert(gelatize(BilinearForm((v,u), dx(u) * v)) == expected)
    # ...

    # ...
    expected = Mass(px,tx)*ny*Stiffness(py,ty)/nx + Mass(py,ty)*nx*Stiffness(px,tx)/ny
    assert(gelatize(BilinearForm((v,u), dot(grad(v), grad(u)))) == expected)
    # ...

    # ...
    expected = (nx*Mass(py,ty)*Stiffness(px,tx)/ny +
                I*Advection(px,tx)*Mass(py,ty)/ny +
                ny*Mass(px,tx)*Stiffness(py,ty)/nx +
                I*Advection(py,ty)*Mass(px,tx)/nx)
    assert(gelatize(BilinearForm((v,u), dot(grad(v), grad(u)) + dx(u)*v + dy(u)*v)) == expected)
    # ...

    # ...
    expected = -bx*I*Advection(px,tx)*Mass(py,ty)/ny - by*I*Advection(py,ty)*Mass(px,tx)/nx
    assert(gelatize(BilinearForm((v,u), dot(b, grad(v)) * u)) == expected)
    # ...

    # ...
    expected = bx**2*nx*Mass(py,ty)*Stiffness(px,tx)/ny + 2*bx*by*Advection(px,tx)*Advection(py,ty) + by**2*ny*Mass(px,tx)*Stiffness(py,ty)/nx
    assert(gelatize(BilinearForm((v,u), dot(b, grad(v)) * dot(b, grad(u)))) == expected)
    # ...

    degrees = None
#    degrees = [2, 1]

#    evaluate = True
    evaluate = False

#    expr = dot(b, grad(v)) * dot(b, grad(u))
#    expr = dot(grad(v), grad(u)) + dx(u)*v + dy(u)*v
#
#    expr = BilinearForm((v,u), expr)
#    print('> input     >>> {0}'.format(expr))
#    print('> gelatized >>> {0}'.format(gelatize(expr, degrees=degrees)))
# ...

# ...
def test_gelatize_2d_3():
    print('============ test_gelatize_2d_3 =============')

    V = VectorFunctionSpace('V', domain)

    v = VectorTestFunction(V, name='v')
    u = VectorTestFunction(V, name='u')

    c = Constant('c')

#    a = BilinearForm((v,u), div(v) * div(u) + rot(v) * rot(u))

    degrees = None

    expr = div(v) * div(u) + c * rot(v) * rot(u)
    expr = BilinearForm((v,u), expr)
    print('> input     >>> {0}'.format(expr))
    print('> gelatized >>> {0}'.format(gelatize(expr, degrees=degrees)))
# ...

# ...
def test_gelatize_2d_4():
    print('============ test_gelatize_2d_4 =============')

    V = FunctionSpace('V', domain)

    v = TestFunction(V, name='v')
    u = TestFunction(V, name='u')

    expr = BilinearForm((v,u), dot(grad(v), grad(u)))

    print('> input     >>> {0}'.format(expr))
    print('> gelatized >>> {0}'.format(gelatize(expr)))
# ...

# ...
def test_gelatize_2d_5():
    print('============ test_gelatize_2d_5 =============')

    V = FunctionSpace('V', domain)

    v = TestFunction(V, name='v')
    u = TestFunction(V, name='u')

    c = Constant('c')

    expr = BilinearForm((v,u), dot(grad(v), grad(u)) + c*v*u)

    print('> input     >>> {0}'.format(expr))
    print('> gelatized >>> {0}'.format(gelatize(expr)))
# ...

# .....................................................
if __name__ == '__main__':
    test_gelatize_2d_1()
#    test_gelatize_2d_3() # TODO debug
    test_gelatize_2d_4()
    test_gelatize_2d_5()
