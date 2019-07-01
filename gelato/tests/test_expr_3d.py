# coding: utf-8

from sympy import Symbol
from sympy.core.containers import Tuple
from sympy import symbols
from sympy import pi, cos, sin
from sympy import srepr
from sympy import I

from sympde.core import Constant
from sympde.calculus import grad, dot, inner, cross, rot, curl, div
from sympde.calculus import laplace, hessian, bracket, convect
from sympde.topology import (dx, dy, dz)
from sympde.topology import FunctionSpace, VectorFunctionSpace
from sympde.topology import Field, TestFunction
from sympde.topology import Domain
from sympde.topology import Trace, trace_0, trace_1
from sympde.topology import Mapping
from sympde.topology import Square
from sympde.expr.expr import LinearForm, BilinearForm

from gelato import gelatize
from gelato import (Mass,
                    Stiffness,
                    Advection,
                    Bilaplacian)

DIM = 3
domain = Domain('Omega', dim=DIM)

#==============================================================================
def test_gelatize_3d_1():

    V = FunctionSpace('V', domain)

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
    assert(gelatize(BilinearForm((u,v), u*v)) == expected)
    # ...

    # ...
    expected = nx*Mass(py,ty)*Mass(pz,tz)*Stiffness(px,tx)/(ny*nz)
    assert(gelatize(BilinearForm((u,v), dx(u)*dx(v))) == expected)
    # ...

    # ...
    expected = I*Advection(py,ty)*Mass(px,tx)*Mass(pz,tz)/(nx*nz)
    assert(gelatize(BilinearForm((u,v), dy(u) * v)) == expected)
    # ...

    # ...
    expected = I*Advection(px,tx)*Mass(py,ty)*Mass(pz,tz)/(ny*nz)
    assert(gelatize(BilinearForm((u,v), dx(u) * v)) == expected)
    # ...

    # ...
    expected = ( nx*Mass(py,ty)*Mass(pz,tz)*Stiffness(px,tx)/(ny*nz) +
                ny*Mass(px,tx)*Mass(pz,tz)*Stiffness(py,ty)/(nx*nz) +
                nz*Mass(px,tx)*Mass(py,ty)*Stiffness(pz,tz)/(nx*ny))
    assert(gelatize(BilinearForm((u,v), dot(grad(v), grad(u)))) == expected)
    # ...

    # ...
    expected = (nx*Mass(py,ty)*Mass(pz,tz)*Stiffness(px,tx)/(ny*nz) +
                I*Advection(px,tx)*Mass(py,ty)*Mass(pz,tz)/(ny*nz) +
                ny*Mass(px,tx)*Mass(pz,tz)*Stiffness(py,ty)/(nx*nz) +
                I*Advection(py,ty)*Mass(px,tx)*Mass(pz,tz)/(nx*nz) +
                nz*Mass(px,tx)*Mass(py,ty)*Stiffness(pz,tz)/(nx*ny))
    assert(gelatize(BilinearForm((u,v), dot(grad(v), grad(u)) + dx(u)*v + dy(u)*v)) == expected)
    # ...

    # ...
    expected = (-bx*I*Advection(px,tx)*Mass(py,ty)*Mass(pz,tz)/(ny*nz) -
                by*I*Advection(py,ty)*Mass(px,tx)*Mass(pz,tz)/(nx*nz) -
                bz*I*Advection(pz,tz)*Mass(px,tx)*Mass(py,ty)/(nx*ny))
    assert(gelatize(BilinearForm((u,v), dot(b, grad(v)) * u)) == expected)
    # ...

#    # ... TODO
#    expected = (bx**2*nx*Mass(py,ty)*Mass(pz,tz)*Stiffness(px,tx)/(ny*nz) +
#                2*bx*by*Advection(px,tx)*Advection(py,ty)*Mass(pz,tz)/nz +
#                2*bx*bz*Advection(px,tx)*Advection(pz,tz)*Mass(py,ty)/ny +
#                by**2*ny*Mass(px,tx)*Mass(pz,tz)*Stiffness(py,ty)/(nx*nz) +
#                2*by*bz*Advection(py,ty)*Advection(pz,tz)*Mass(px,tx)/nx +
#                bz**2*nz*Mass(px,tx)*Mass(py,ty)*Stiffness(pz,tz)/(nx*ny))
#    assert(gelatize(BilinearForm((u,v), dot(b, grad(v)) * dot(b, grad(u)))) == expected)
#    # ...

    degrees = None
#    degrees = [2, 1, 1]

#    evaluate = True
    evaluate = False

#    expr = u*v
#
#    expr = BilinearForm((u,v), expr)
#    print('> input     >>> {0}'.format(expr))
#    print('> gelatized >>> {0}'.format(gelatize(expr, degrees, evaluate=evaluate)))

#==============================================================================
def test_gelatize_3d_5_mapping():

    M = Mapping('M', DIM)

    V = FunctionSpace('V', domain)

    v = TestFunction(V, name='v')
    u = TestFunction(V, name='u')

    c = Constant('c')

    expr = BilinearForm((u,v), dot(grad(v), grad(u)) + c*v*u)

    print('> input     >>> {0}'.format(expr))
    print('> gelatized >>> {0}'.format(gelatize(expr, mapping=M, human=True)))

#==============================================================================
# TODO it takes some time => optimize LogicalExpr
#      using nodes to describe subexpresions
def test_gelatize_3d_3_mapping():

    M = Mapping('M', DIM)

    V = VectorFunctionSpace('V', domain)

    v = TestFunction(V, name='v')
    u = TestFunction(V, name='u')

    c = Constant('c')

    expr = c * div(v) * div(u) + dot(curl(v), curl(u))
    expr = BilinearForm((u,v), expr)
    print('> input     >>> {0}'.format(expr))
    print('> gelatized >>> {0}'.format(gelatize(expr, mapping=M, human=True)))

#==============================================================================
# CLEAN UP SYMPY NAMESPACE
#==============================================================================

def teardown_module():
    from sympy import cache
    cache.clear_cache()

def teardown_function():
    from sympy import cache
    cache.clear_cache()
