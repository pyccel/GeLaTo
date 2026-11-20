# coding: utf-8

from sympy import Symbol
from sympy.core.containers import Tuple
from sympy import symbols
from sympy import pi, cos, sin
from sympy import srepr
from sympy import I

from sympde.api import Constant
from sympde.api import grad, dot, inner, cross, rot, curl, div
from sympde.api import laplace, hessian, bracket, convect
from sympde.api import dx, dy, dz
from sympde.api import ScalarFunctionSpace, VectorFunctionSpace
from sympde.api import Domain
from sympde.api import Mapping
from sympde.api import elements_of
from sympde.api import BilinearForm
from sympde.api import integral

from gelato import GltExpr

DIM = 2

#==============================================================================
def test_glt_expr_2d_1():
    domain = Domain('Omega', dim=DIM)

    V = ScalarFunctionSpace('V', domain)

    u,v = elements_of(V, names='u,v')

    c = Constant('c')

    expr = dot(grad(v), grad(u)) + c*v*u
    a = BilinearForm((u,v), integral(domain, expr))

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
