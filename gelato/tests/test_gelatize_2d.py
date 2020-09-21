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
from sympde.topology import dx1, dx2, dx3
from sympde.topology import ScalarFunctionSpace, VectorFunctionSpace
from sympde.topology import Domain
from sympde.topology import Mapping
from sympde.topology import elements_of
from sympde.expr import BilinearForm
from sympde.expr import integral

from gelato import gelatize
from gelato import (Mass,
                    Stiffness,
                    Advection,
                    Bilaplacian)

DIM = 2

#==============================================================================
def test_gelatize_2d_1():
    domain = Domain('Omega', dim=DIM)

    V = ScalarFunctionSpace('V', domain)

    u,v = elements_of(V, names='u,v')

    nx, ny = symbols('nx ny', integer=True)
    px, py = symbols('px py', integer=True)
    tx, ty = symbols('tx ty')

    expected = Mass(px,tx)*Mass(py,ty)/(nx*ny)

    expr = u*v
    expr = BilinearForm((u,v), integral(domain, expr))
    assert(gelatize(expr) == expected)


#==============================================================================
def test_gelatize_2d_2():
    domain = Domain('Omega', dim=DIM)

    V = ScalarFunctionSpace('V', domain)

    u,v = elements_of(V, names='u,v')

    nx, ny = symbols('nx ny', integer=True)
    px, py = symbols('px py', integer=True)
    tx, ty = symbols('tx ty')

    expected = Mass(py,ty)*nx*Stiffness(px,tx)/ny

    expr = dx1(u)*dx1(v)
    expr = BilinearForm((u,v), integral(domain, expr))
    assert(gelatize(expr) == expected)

#==============================================================================
def test_gelatize_2d_3():
    domain = Domain('Omega', dim=DIM)

    V = ScalarFunctionSpace('V', domain)

    u,v = elements_of(V, names='u,v')

    nx, ny = symbols('nx ny', integer=True)
    px, py = symbols('px py', integer=True)
    tx, ty = symbols('tx ty')

    expected = I*Advection(py,ty)*Mass(px,tx)/nx

    expr = dx2(u) * v
    expr = BilinearForm((u,v), integral(domain, expr))
    assert(gelatize(expr) == expected)

#==============================================================================
def test_gelatize_2d_4():
    domain = Domain('Omega', dim=DIM)

    V = ScalarFunctionSpace('V', domain)

    u,v = elements_of(V, names='u,v')

    nx, ny = symbols('nx ny', integer=True)
    px, py = symbols('px py', integer=True)
    tx, ty = symbols('tx ty')

    expected = I*Advection(px,tx)*Mass(py,ty)/ny

    expr = dx1(u) * v
    expr = BilinearForm((u,v), integral(domain, expr))
    assert(gelatize(expr) == expected)

#==============================================================================
def test_gelatize_2d_5():
    domain = Domain('Omega', dim=DIM)

    V = ScalarFunctionSpace('V', domain)

    u,v = elements_of(V, names='u,v')

    nx, ny = symbols('nx ny', integer=True)
    px, py = symbols('px py', integer=True)
    tx, ty = symbols('tx ty')

    expected = Mass(px,tx)*ny*Stiffness(py,ty)/nx + Mass(py,ty)*nx*Stiffness(px,tx)/ny

    expr = dot(grad(v), grad(u))
    expr = BilinearForm((u,v), integral(domain, expr))
    assert(gelatize(expr) == expected)

#==============================================================================
def test_gelatize_2d_6():
    domain = Domain('Omega', dim=DIM)

    V = ScalarFunctionSpace('V', domain)

    u,v = elements_of(V, names='u,v')

    nx, ny = symbols('nx ny', integer=True)
    px, py = symbols('px py', integer=True)
    tx, ty = symbols('tx ty')

    expected = (nx*Mass(py,ty)*Stiffness(px,tx)/ny +
                I*Advection(px,tx)*Mass(py,ty)/ny +
                ny*Mass(px,tx)*Stiffness(py,ty)/nx +
                I*Advection(py,ty)*Mass(px,tx)/nx)

    expr = dot(grad(v), grad(u)) + dx1(u)*v + dx2(u)*v
    expr = BilinearForm((u,v), integral(domain, expr))
    assert(gelatize(expr) == expected)

#==============================================================================
def test_gelatize_2d_7():
    domain = Domain('Omega', dim=DIM)

    V = ScalarFunctionSpace('V', domain)

    u,v = elements_of(V, names='u,v')

    nx, ny = symbols('nx ny', integer=True)
    px, py = symbols('px py', integer=True)
    tx, ty = symbols('tx ty')

    bx = Constant('bx')
    by = Constant('by')
    b = Tuple(bx, by)

    expected = -bx*I*Advection(px,tx)*Mass(py,ty)/ny - by*I*Advection(py,ty)*Mass(px,tx)/nx

    expr = dot(b, grad(v)) * u
    expr = BilinearForm((u,v), integral(domain, expr))
    assert(gelatize(expr) == expected)

#==============================================================================
## TODO
#def test_gelatize_2d_8():
#    domain = Domain('Omega', dim=DIM)
#
#    V = ScalarFunctionSpace('V', domain)
#
#    u,v = elements_of(V, names='u,v')
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
#    expr = BilinearForm((u,v), integral(domain, expr))
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
