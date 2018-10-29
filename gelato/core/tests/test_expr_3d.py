# coding: utf-8

from sympy import Symbol
from sympy.core.containers import Tuple
from sympy import symbols
from sympy import pi, cos, sin
from sympy import srepr
from sympy import I

from sympde.core import dx, dy, dz
from sympde.core import Constant
from sympde.core import Field
from sympde.core import grad, dot, inner, cross, rot, curl, div
from sympde.core import H1Space
from sympde.core import TestFunction
from sympde.core import VectorTestFunction
from sympde.core import BilinearForm

from gelato.core import gelatize
from gelato.core import (Mass,
                         Stiffness,
                         Advection,
                         Bilaplacian)

# ...
def test_gelatize_3d_1():
    print('============ test_gelatize_3d_1 =============')

    V = H1Space('V', ldim=3)

    v = TestFunction(V, name='v')
    u = TestFunction(V, name='u')

    nx, ny, nz = symbols('nx ny nz', integer=True)
    px, py, pz = symbols('px py pz', integer=True)
    tx, ty, tz = symbols('tx ty tz')

    c = Constant('c')

    bx = Constant('bx')
    by = Constant('by')
    bz = Constant('bz')
    b = Tuple(bx, by, bz)

    # ...
    expected = Mass(px,tx)*Mass(py,ty)*Mass(pz,tz)/(nx*ny*nz)
    assert(gelatize(BilinearForm((v,u), u*v)) == expected)
    # ...

    # ...
    expected = nx*Mass(py,ty)*Mass(pz,tz)*Stiffness(px,tx)/(ny*nz)
    assert(gelatize(BilinearForm((v,u), dx(u)*dx(v))) == expected)
    # ...

    # ...
    expected = I*Advection(py,ty)*Mass(px,tx)*Mass(pz,tz)/(nx*nz)
    assert(gelatize(BilinearForm((v,u), dy(u) * v)) == expected)
    # ...

    # ...
    expected = I*Advection(px,tx)*Mass(py,ty)*Mass(pz,tz)/(ny*nz)
    assert(gelatize(BilinearForm((v,u), dx(u) * v)) == expected)
    # ...

    # ...
    expected = ( nx*Mass(py,ty)*Mass(pz,tz)*Stiffness(px,tx)/(ny*nz) +
                ny*Mass(px,tx)*Mass(pz,tz)*Stiffness(py,ty)/(nx*nz) +
                nz*Mass(px,tx)*Mass(py,ty)*Stiffness(pz,tz)/(nx*ny))
    assert(gelatize(BilinearForm((v,u), dot(grad(v), grad(u)))) == expected)
    # ...

    # ...
    expected = (nx*Mass(py,ty)*Mass(pz,tz)*Stiffness(px,tx)/(ny*nz) +
                I*Advection(px,tx)*Mass(py,ty)*Mass(pz,tz)/(ny*nz) +
                ny*Mass(px,tx)*Mass(pz,tz)*Stiffness(py,ty)/(nx*nz) +
                I*Advection(py,ty)*Mass(px,tx)*Mass(pz,tz)/(nx*nz) +
                nz*Mass(px,tx)*Mass(py,ty)*Stiffness(pz,tz)/(nx*ny))
    assert(gelatize(BilinearForm((v,u), dot(grad(v), grad(u)) + dx(u)*v + dy(u)*v)) == expected)
    # ...

    # ...
    expected = (-bx*I*Advection(px,tx)*Mass(py,ty)*Mass(pz,tz)/(ny*nz) -
                by*I*Advection(py,ty)*Mass(px,tx)*Mass(pz,tz)/(nx*nz) -
                bz*I*Advection(pz,tz)*Mass(px,tx)*Mass(py,ty)/(nx*ny))
    assert(gelatize(BilinearForm((v,u), dot(b, grad(v)) * u)) == expected)
    # ...

    # ...
    expected = (bx**2*nx*Mass(py,ty)*Mass(pz,tz)*Stiffness(px,tx)/(ny*nz) +
                2*bx*by*Advection(px,tx)*Advection(py,ty)*Mass(pz,tz)/nz +
                2*bx*bz*Advection(px,tx)*Advection(pz,tz)*Mass(py,ty)/ny +
                by**2*ny*Mass(px,tx)*Mass(pz,tz)*Stiffness(py,ty)/(nx*nz) +
                2*by*bz*Advection(py,ty)*Advection(pz,tz)*Mass(px,tx)/nx +
                bz**2*nz*Mass(px,tx)*Mass(py,ty)*Stiffness(pz,tz)/(nx*ny))
    assert(gelatize(BilinearForm((v,u), dot(b, grad(v)) * dot(b, grad(u)))) == expected)
    # ...

    degrees = None
#    degrees = [2, 1, 1]

#    evaluate = True
    evaluate = False

#    expr = u*v
#
#    expr = BilinearForm((v,u), expr)
#    print('> input     >>> {0}'.format(expr))
#    print('> gelatized >>> {0}'.format(gelatize(expr, degrees, evaluate=evaluate)))
# ...


# .....................................................
if __name__ == '__main__':

    test_gelatize_3d_1()
