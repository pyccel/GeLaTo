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
from sympde.topology import ScalarFunctionSpace, VectorFunctionSpace
from sympde.topology import TestFunction
#from sympde.topology import element_of # TODO not working yet
from sympde.topology import Domain
from sympde.topology import Mapping
from sympde.expr.expr import BilinearForm

from gelato import gelatize, GltExpr
from gelato import (Mass,
                    Stiffness,
                    Advection,
                    Bilaplacian)

DIM = 2

#==============================================================================
def test_gelatize_2d_1():
    domain = Domain('Omega', dim=DIM)

    V = ScalarFunctionSpace('V', domain)

    v = TestFunction(V, name='v')
    u = TestFunction(V, name='u')

    nx, ny = symbols('nx ny', integer=True)
    px, py = symbols('px py', integer=True)
    tx, ty = symbols('tx ty')

    expected = Mass(px,tx)*Mass(py,ty)/(nx*ny)
    assert(gelatize(BilinearForm((u,v), u*v)) == expected)

#==============================================================================
def test_gelatize_2d_2():
    domain = Domain('Omega', dim=DIM)

    V = ScalarFunctionSpace('V', domain)

    v = TestFunction(V, name='v')
    u = TestFunction(V, name='u')

    nx, ny = symbols('nx ny', integer=True)
    px, py = symbols('px py', integer=True)
    tx, ty = symbols('tx ty')

    expected = Mass(py,ty)*nx*Stiffness(px,tx)/ny
    assert(gelatize(BilinearForm((u,v), dx(u)*dx(v))) == expected)

#==============================================================================
def test_gelatize_2d_3():
    domain = Domain('Omega', dim=DIM)

    V = ScalarFunctionSpace('V', domain)

    v = TestFunction(V, name='v')
    u = TestFunction(V, name='u')

    nx, ny = symbols('nx ny', integer=True)
    px, py = symbols('px py', integer=True)
    tx, ty = symbols('tx ty')

    expected = I*Advection(py,ty)*Mass(px,tx)/nx
    assert(gelatize(BilinearForm((u,v), dy(u) * v)) == expected)

#==============================================================================
def test_gelatize_2d_4():
    domain = Domain('Omega', dim=DIM)

    V = ScalarFunctionSpace('V', domain)

    v = TestFunction(V, name='v')
    u = TestFunction(V, name='u')

    nx, ny = symbols('nx ny', integer=True)
    px, py = symbols('px py', integer=True)
    tx, ty = symbols('tx ty')

    expected = I*Advection(px,tx)*Mass(py,ty)/ny
    assert(gelatize(BilinearForm((u,v), dx(u) * v)) == expected)

#==============================================================================
def test_gelatize_2d_5():
    domain = Domain('Omega', dim=DIM)

    V = ScalarFunctionSpace('V', domain)

    v = TestFunction(V, name='v')
    u = TestFunction(V, name='u')

    nx, ny = symbols('nx ny', integer=True)
    px, py = symbols('px py', integer=True)
    tx, ty = symbols('tx ty')

    expected = Mass(px,tx)*ny*Stiffness(py,ty)/nx + Mass(py,ty)*nx*Stiffness(px,tx)/ny
    assert(gelatize(BilinearForm((u,v), dot(grad(v), grad(u)))) == expected)

#==============================================================================
def test_gelatize_2d_6():
    domain = Domain('Omega', dim=DIM)

    V = ScalarFunctionSpace('V', domain)

    v = TestFunction(V, name='v')
    u = TestFunction(V, name='u')

    nx, ny = symbols('nx ny', integer=True)
    px, py = symbols('px py', integer=True)
    tx, ty = symbols('tx ty')

    expected = (nx*Mass(py,ty)*Stiffness(px,tx)/ny +
                I*Advection(px,tx)*Mass(py,ty)/ny +
                ny*Mass(px,tx)*Stiffness(py,ty)/nx +
                I*Advection(py,ty)*Mass(px,tx)/nx)
    assert(gelatize(BilinearForm((u,v), dot(grad(v), grad(u)) + dx(u)*v + dy(u)*v)) == expected)

#==============================================================================
def test_gelatize_2d_7():
    domain = Domain('Omega', dim=DIM)

    V = ScalarFunctionSpace('V', domain)

    v = TestFunction(V, name='v')
    u = TestFunction(V, name='u')

    nx, ny = symbols('nx ny', integer=True)
    px, py = symbols('px py', integer=True)
    tx, ty = symbols('tx ty')

    bx = Constant('bx')
    by = Constant('by')
    b = Tuple(bx, by)

    expected = -bx*I*Advection(px,tx)*Mass(py,ty)/ny - by*I*Advection(py,ty)*Mass(px,tx)/nx
    assert(gelatize(BilinearForm((u,v), dot(b, grad(v)) * u)) == expected)

#==============================================================================
## TODO
#def test_gelatize_2d_8():
#    domain = Domain('Omega', dim=DIM)
#
#    V = ScalarFunctionSpace('V', domain)
#
#    v = TestFunction(V, name='v')
#    u = TestFunction(V, name='u')
#
#    nx, ny = symbols('nx ny', integer=True)
#    px, py = symbols('px py', integer=True)
#    tx, ty = symbols('tx ty')
#
#    bx = Constant('bx')
#    by = Constant('by')
#    b = Tuple(bx, by)
#
#    # ...
#    expected = bx**2*nx*Mass(py,ty)*Stiffness(px,tx)/ny + 2*bx*by*Advection(px,tx)*Advection(py,ty) + by**2*ny*Mass(px,tx)*Stiffness(py,ty)/n
#    assert(gelatize(BilinearForm((u,v), dot(b, grad(v)) * dot(b, grad(u)))) == expected)
#    # ...
#
#    degrees = None
##    degrees = [2, 1]
#
##    evaluate = True
#    evaluate = False
#
#    expr = dot(b, grad(v)) * dot(b, grad(u))
#    expr = dot(grad(v), grad(u)) + dx(u)*v + dy(u)*v
#
#    expr = BilinearForm((u,v), expr)
#    print('> input     >>> {0}'.format(expr))
#    print('> gelatized >>> {0}'.format(gelatize(expr, degrees=degrees)))

#==============================================================================
# CLEAN UP SYMPY NAMESPACE
#==============================================================================

def teardown_module():
    from sympy import cache
    cache.clear_cache()

def teardown_function():
    from sympy import cache
    cache.clear_cache()
