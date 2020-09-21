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

DIM = 3

#==============================================================================
def test_gelatize_3d_1():
    domain = Domain('Omega', dim=DIM)

    V = ScalarFunctionSpace('V', domain)

    u,v = elements_of(V, names='u,v')

    nx, ny, nz = symbols('nx ny nz', integer=True)
    px, py, pz = symbols('px py pz', integer=True)
    tx, ty, tz = symbols('tx ty tz')

    expected = Mass(px,tx)*Mass(py,ty)*Mass(pz,tz)/(nx*ny*nz)

    expr = u*v
    expr = BilinearForm((u,v), integral(domain, expr))
    assert(gelatize(expr) == expected)

#==============================================================================
def test_gelatize_3d_2():
    domain = Domain('Omega', dim=DIM)

    V = ScalarFunctionSpace('V', domain)

    u,v = elements_of(V, names='u,v')

    nx, ny, nz = symbols('nx ny nz', integer=True)
    px, py, pz = symbols('px py pz', integer=True)
    tx, ty, tz = symbols('tx ty tz')

    expected = nx*Mass(py,ty)*Mass(pz,tz)*Stiffness(px,tx)/(ny*nz)

    expr = dx1(u)*dx1(v)
    expr = BilinearForm((u,v), integral(domain, expr))
    assert(gelatize(expr) == expected)

#==============================================================================
def test_gelatize_3d_3():
    domain = Domain('Omega', dim=DIM)

    V = ScalarFunctionSpace('V', domain)

    u,v = elements_of(V, names='u,v')

    nx, ny, nz = symbols('nx ny nz', integer=True)
    px, py, pz = symbols('px py pz', integer=True)
    tx, ty, tz = symbols('tx ty tz')

    expected = I*Advection(py,ty)*Mass(px,tx)*Mass(pz,tz)/(nx*nz)

    expr = dx2(u) * v
    expr = BilinearForm((u,v), integral(domain, expr))
    assert(gelatize(expr) == expected)

#==============================================================================
def test_gelatize_3d_4():
    domain = Domain('Omega', dim=DIM)

    V = ScalarFunctionSpace('V', domain)

    u,v = elements_of(V, names='u,v')

    nx, ny, nz = symbols('nx ny nz', integer=True)
    px, py, pz = symbols('px py pz', integer=True)
    tx, ty, tz = symbols('tx ty tz')

    expected = I*Advection(px,tx)*Mass(py,ty)*Mass(pz,tz)/(ny*nz)

    expr = dx1(u) * v
    expr = BilinearForm((u,v), integral(domain, expr))
    assert(gelatize(expr) == expected)

#==============================================================================
def test_gelatize_3d_5():
    domain = Domain('Omega', dim=DIM)

    V = ScalarFunctionSpace('V', domain)

    u,v = elements_of(V, names='u,v')

    nx, ny, nz = symbols('nx ny nz', integer=True)
    px, py, pz = symbols('px py pz', integer=True)
    tx, ty, tz = symbols('tx ty tz')

    expected = ( nx*Mass(py,ty)*Mass(pz,tz)*Stiffness(px,tx)/(ny*nz) +
                ny*Mass(px,tx)*Mass(pz,tz)*Stiffness(py,ty)/(nx*nz) +
                nz*Mass(px,tx)*Mass(py,ty)*Stiffness(pz,tz)/(nx*ny))

    expr = dot(grad(v), grad(u))
    expr = BilinearForm((u,v), integral(domain, expr))
    assert(gelatize(expr) == expected)

#==============================================================================
def test_gelatize_3d_6():
    domain = Domain('Omega', dim=DIM)

    V = ScalarFunctionSpace('V', domain)

    u,v = elements_of(V, names='u,v')

    nx, ny, nz = symbols('nx ny nz', integer=True)
    px, py, pz = symbols('px py pz', integer=True)
    tx, ty, tz = symbols('tx ty tz')

    expected = (nx*Mass(py,ty)*Mass(pz,tz)*Stiffness(px,tx)/(ny*nz) +
                I*Advection(px,tx)*Mass(py,ty)*Mass(pz,tz)/(ny*nz) +
                ny*Mass(px,tx)*Mass(pz,tz)*Stiffness(py,ty)/(nx*nz) +
                I*Advection(py,ty)*Mass(px,tx)*Mass(pz,tz)/(nx*nz) +
                nz*Mass(px,tx)*Mass(py,ty)*Stiffness(pz,tz)/(nx*ny))

    expr = dot(grad(v), grad(u)) + dx1(u)*v + dx2(u)*v
    expr = BilinearForm((u,v), integral(domain, expr))
    assert(gelatize(expr) == expected)

#==============================================================================
def test_gelatize_3d_7():
    domain = Domain('Omega', dim=DIM)

    V = ScalarFunctionSpace('V', domain)

    u,v = elements_of(V, names='u,v')

    nx, ny, nz = symbols('nx ny nz', integer=True)
    px, py, pz = symbols('px py pz', integer=True)
    tx, ty, tz = symbols('tx ty tz')

    bx = Constant('bx')
    by = Constant('by')
    bz = Constant('bz')
    b = Tuple(bx, by, bz)

    expected = (-bx*I*Advection(px,tx)*Mass(py,ty)*Mass(pz,tz)/(ny*nz) -
                by*I*Advection(py,ty)*Mass(px,tx)*Mass(pz,tz)/(nx*nz) -
                bz*I*Advection(pz,tz)*Mass(px,tx)*Mass(py,ty)/(nx*ny))

    expr = dot(b, grad(v)) * u
    expr = BilinearForm((u,v), integral(domain, expr))
    assert(gelatize(expr) == expected)


#==============================================================================
## TODO
#def test_gelatize_3d_8():
#    domain = Domain('Omega', dim=DIM)
#
#    V = ScalarFunctionSpace('V', domain)
#
#    u,v = elements_of(V, names='u,v')
#
#    nx, ny, nz = symbols('nx ny nz', integer=True)
#    px, py, pz = symbols('px py pz', integer=True)
#    tx, ty, tz = symbols('tx ty tz')
#
#    c = Constant('c')
#
#    bx = Constant('bx')
#    by = Constant('by')
#    bz = Constant('bz')
#    b = Tuple(bx, by, bz)
#
#    # ...
#    expected = (bx**2*nx*Mass(py,ty)*Mass(pz,tz)*Stiffness(px,tx)/(ny*nz) +
#                2*bx*by*Advection(px,tx)*Advection(py,ty)*Mass(pz,tz)/nz +
#                2*bx*bz*Advection(px,tx)*Advection(pz,tz)*Mass(py,ty)/ny +
#                by**2*ny*Mass(px,tx)*Mass(pz,tz)*Stiffness(py,ty)/(nx*nz) +
#                2*by*bz*Advection(py,ty)*Advection(pz,tz)*Mass(px,tx)/nx +
#                bz**2*nz*Mass(px,tx)*Mass(py,ty)*Stiffness(pz,tz)/(nx*ny))
#    expr = dot(b, grad(v)) * dot(b, grad(u))
#    expr = BilinearForm((u,v), integral(domain, expr))
#    assert(gelatize(expr) == expected)
#    # ...
#
##    degrees = None
#    degrees = [2, 1, 1]
#
#    evaluate = True
##    evaluate = False
#
#    expr = u*v
#
#    expr = BilinearForm((u,v), expr)
#    print('> input     >>> {0}'.format(expr))
#    print('> gelatized >>> {0}'.format(gelatize(expr, degrees, evaluate=evaluate)))

#==============================================================================
# CLEAN UP SYMPY NAMESPACE
#==============================================================================

def teardown_module():
    from sympy import cache
    cache.clear_cache()

def teardown_function():
    from sympy import cache
    cache.clear_cache()
