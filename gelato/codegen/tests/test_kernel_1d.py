# coding: utf-8

import os

from sympy import Symbol
from sympy.core.containers import Tuple
from sympy import symbols
from sympy import IndexedBase
from sympy import Matrix
from sympy import Function
from sympy import pi, cos, sin
from sympy import S

from sympde.core import dx, dy, dz
from sympde.core import Constant
from sympde.core import Field
from sympde.core import grad, dot, inner, cross, rot, curl, div
from sympde.core import FunctionSpace
from sympde.core import TestFunction
from sympde.core import VectorTestFunction
from sympde.core import BilinearForm, LinearForm, Integral
from sympde.core import Mapping
from sympde.core import Domain
from sympde.core import Boundary, trace_0, trace_1
from sympde.core import evaluate

from gelato.codegen.ast import Kernel

from spl.api.codegen.printing import pycode
from spl.fem.splines import SplineSpace

from numpy import linspace, zeros, allclose

# ...
def create_discrete_space():
    # ... discrete spaces
    # Input data: degree, number of elements
    p  = 1
    ne = 2**1

    # Create uniform grid
    grid = linspace( 0., 1., num=ne+1 )

    # Create finite element space and precompute quadrature data
    V = SplineSpace( p, grid=grid )
    V.init_fem()
    # ...

    return V
# ...


#DEBUG = False
DEBUG = True
DIM = 1

domain = Domain('\Omega', dim=DIM)

def test_kernel_1d_scalar_1(mapping=False):
    print('============ test_kernel_1d_scalar_1 =============')

    if mapping: mapping = Mapping('M', rdim=DIM, domain=domain)

    U = FunctionSpace('U', domain)
    V = FunctionSpace('V', domain)

    v = TestFunction(V, name='v')
    u = TestFunction(U, name='u')

    expr = dot(grad(v), grad(u))
    a = BilinearForm((v,u), expr, mapping=mapping)

    # ... discrete spaces
    Vh = create_discrete_space()
    # ...

    kernel = Kernel(a, Vh, name='kernel')
    code = pycode(kernel.func)
    if DEBUG: print(code)

def test_kernel_1d_scalar_2(mapping=False):
    print('============ test_kernel_1d_scalar_2 =============')

    if mapping: mapping = Mapping('M', rdim=DIM, domain=domain)

    U = FunctionSpace('U', domain)
    V = FunctionSpace('V', domain)

    v = TestFunction(V, name='v')
    u = TestFunction(U, name='u')

    c = Constant('c', real=True)

    expr = dot(grad(v), grad(u)) + c*v*u
    a = BilinearForm((v,u), expr, mapping=mapping)

    # ... discrete spaces
    Vh = create_discrete_space()
    # ...


    kernel = Kernel(a, Vh, name='kernel')
    code = pycode(kernel.func)
    if DEBUG: print(code)


#................................
if __name__ == '__main__':

    # .................................
    # without mapping
    test_kernel_1d_scalar_1(mapping=False)
    test_kernel_1d_scalar_2(mapping=False)
    # .................................
