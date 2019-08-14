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
from sympde.topology import ScalarFunctionSpace
from sympde.topology import TestFunction
#from sympde.topology import element_of # TODO not working yet
from sympde.topology import Domain
from sympde.expr.expr import BilinearForm

from gelato import gelatize
from gelato import (Mass,
                    Stiffness,
                    Advection,
                    Bilaplacian)

DIM = 1

#==============================================================================
def test_gelatize_1d_1():
    domain = Domain('Omega', dim=DIM)

    V = ScalarFunctionSpace('V', domain)

    v = TestFunction(V, name='v')
    u = TestFunction(V, name='u')

    nx = symbols('nx', integer=True)
    px = symbols('px', integer=True)
    tx = symbols('tx')

    expected = Mass(px,tx)/nx
    assert(gelatize(BilinearForm((u,v), u*v)) == expected)

#==============================================================================
def test_gelatize_1d_2():
    domain = Domain('Omega', dim=DIM)

    V = ScalarFunctionSpace('V', domain)

    v = TestFunction(V, name='v')
    u = TestFunction(V, name='u')

    nx = symbols('nx', integer=True)
    px = symbols('px', integer=True)
    tx = symbols('tx')

    expected = nx*Stiffness(px,tx)
    assert(gelatize(BilinearForm((u,v), dot(grad(v), grad(u)))) == expected)

#==============================================================================
def test_gelatize_1d_3():
    domain = Domain('Omega', dim=DIM)

    V = ScalarFunctionSpace('V', domain)

    v = TestFunction(V, name='v')
    u = TestFunction(V, name='u')

    nx = symbols('nx', integer=True)
    px = symbols('px', integer=True)
    tx = symbols('tx')

    expected = I*Advection(px,tx)
    assert(gelatize(BilinearForm((u,v), dx(u) * v)) == expected)

#==============================================================================
def test_gelatize_1d_4():
    domain = Domain('Omega', dim=DIM)

    V = ScalarFunctionSpace('V', domain)

    v = TestFunction(V, name='v')
    u = TestFunction(V, name='u')

    nx = symbols('nx', integer=True)
    px = symbols('px', integer=True)
    tx = symbols('tx')

    c1 = Constant('c1')
    c2 = Constant('c2')
    c3 = Constant('c3')
    c4 = Constant('c4')

    expected = c1*Mass(px,tx)/nx + c2*I*Advection(px,tx) - c3*I*Advection(px,tx) + c4*nx*Stiffness(px,tx)
    assert(gelatize(BilinearForm((u,v), c1*v*u + c2*dx(u)*v + c3*dx(v)*u + c4*dx(v)*dx(u))) == expected)

#==============================================================================
# CLEAN UP SYMPY NAMESPACE
#==============================================================================

def teardown_module():
    from sympy import cache
    cache.clear_cache()

def teardown_function():
    from sympy import cache
    cache.clear_cache()
