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
def test_glt_expr_2d_1():
    domain = Domain('Omega', dim=DIM)

    V = ScalarFunctionSpace('V', domain)

    v = TestFunction(V, name='v')
    u = TestFunction(V, name='u')

    c = Constant('c')

    a = BilinearForm((u,v), dot(grad(v), grad(u)) + c*v*u)
    glt = GltExpr(a)
    print(glt)
    print(glt(degrees=[2,2]))
    print(glt(tx=0.1, ty=0.2, degrees=[2,2]))


#==============================================================================
# CLEAN UP SYMPY NAMESPACE
#==============================================================================

def teardown_module():
    from sympy import cache
    cache.clear_cache()

def teardown_function():
    from sympy import cache
    cache.clear_cache()
