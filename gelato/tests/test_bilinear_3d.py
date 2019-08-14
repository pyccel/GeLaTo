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

from gelato import gelatize
from gelato import (Mass,
                    Stiffness,
                    Advection,
                    Bilaplacian)

DIM = 3

#==============================================================================
def test_bilinear_3d_mapping_1():
    domain = Domain('Omega', dim=DIM)

    M = Mapping('M', DIM)

    V = ScalarFunctionSpace('V', domain)

    v = TestFunction(V, name='v')
    u = TestFunction(V, name='u')

    c = Constant('c')

    expr = BilinearForm((u,v), dot(grad(v), grad(u)) + c*v*u)

    print('> input     >>> {0}'.format(expr))
    print('> gelatized >>> {0}'.format(gelatize(expr, mapping=M, human=True)))

#==============================================================================
# TODO it takes some time => optimize LogicalExpr
#      using nodes to describe subexpresions
def test_bilinear_3d_mapping_2():
    domain = Domain('Omega', dim=DIM)

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
