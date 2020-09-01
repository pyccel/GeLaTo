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
from sympde.topology import Domain
from sympde.topology import Mapping
from sympde.topology import elements_of
from sympde.expr import BilinearForm
from sympde.expr import integral

from gelato import gelatize

DIM = 2

#==============================================================================
def test_bilinear_2d_1():
    domain = Domain('Omega', dim=DIM)

    V = VectorFunctionSpace('V', domain)

    u,v = elements_of(V, names='u,v')

    c = Constant('c')

    expr = c * div(v) * div(u) + curl(v) * curl(u)
    expr = BilinearForm((u,v), integral(domain, expr))

    print('> input     >>> {0}'.format(expr))
    print('> gelatized >>> {0}'.format(gelatize(expr)))

#==============================================================================
def test_bilinear_2d_2():
    domain = Domain('Omega', dim=DIM)

    V = ScalarFunctionSpace('V', domain)

    u,v = elements_of(V, names='u,v')

    expr = dot(grad(v), grad(u))
    expr = BilinearForm((u,v), integral(domain, expr))

    print('> input     >>> {0}'.format(expr))
    print('> gelatized >>> {0}'.format(gelatize(expr)))

#==============================================================================
def test_bilinear_2d_3():
    domain = Domain('Omega', dim=DIM)

    V = ScalarFunctionSpace('V', domain)

    u,v = elements_of(V, names='u,v')

    c = Constant('c')

    expr = dot(grad(v), grad(u)) + c*v*u
    expr = BilinearForm((u,v), integral(domain, expr))

    print('> input     >>> {0}'.format(expr))
    print('> gelatized >>> {0}'.format(gelatize(expr)))

#==============================================================================
def test_bilinear_2d_mapping_1():
    domain = Domain('Omega', dim=DIM)

    M = Mapping('M', DIM)

    mapped_domain = M(domain)

    V = ScalarFunctionSpace('V', mapped_domain)

    u,v = elements_of(V, names='u,v')

    c = Constant('c')

    expr = dot(grad(v), grad(u)) + c*v*u
    expr = BilinearForm((u,v), integral(mapped_domain, expr))

    print('> input     >>> {0}'.format(expr))
    print('> gelatized >>> {0}'.format(gelatize(expr, mapping=M, human=True)))

#==============================================================================
def test_bilinear_2d_mapping_2():
    domain = Domain('Omega', dim=DIM)

    M = Mapping('M', DIM)

    mapped_domain = M(domain)

    V = VectorFunctionSpace('V', mapped_domain)

    u,v = elements_of(V, names='u,v')

    c = Constant('c')

    expr = c * div(v) * div(u) + curl(v) * curl(u)
    expr = BilinearForm((u,v), integral(mapped_domain, expr))

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
