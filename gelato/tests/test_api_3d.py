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

from gelato.api import DiscreteSymbol
from gelato.printing.pycode import pycode

from spl.fem.splines import SplineSpace
from spl.fem.tensor  import TensorFemSpace

from numpy import linspace, zeros, allclose

# ...
def create_discrete_space():
    # ... discrete spaces
    # Input data: degree, number of elements
    p1  = 1 ; p2  = 1 ; p3  = 1
    ne1 = 2**1 ; ne2 = 2**1 ; ne3 = 2**1

    # Create uniform grid
    grid_1 = linspace( 0., 1., num=ne1+1 )
    grid_2 = linspace( 0., 1., num=ne2+1 )
    grid_3 = linspace( 0., 1., num=ne3+1 )

    # Create 1D finite element spaces and precompute quadrature data
    V1 = SplineSpace( p1, grid=grid_1 ); V1.init_fem()
    V2 = SplineSpace( p2, grid=grid_2 ); V2.init_fem()
    V3 = SplineSpace( p3, grid=grid_3 ); V3.init_fem()

    # Create 3D tensor product finite element space
    V = TensorFemSpace( V1, V2, V3 )
    # ...

    return V
# ...


#DEBUG = False
DEBUG = True
DIM = 3

domain = Domain('\Omega', dim=DIM)

def test_interface_3d_scalar_1(mapping=False):
    print('============ test_interface_3d_scalar_1 =============')

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

    symbol = DiscreteSymbol(a, Vh)

    t1 = linspace(0, 1, 100)
    t2 = linspace(0, 1, 100)
    t3 = linspace(0, 1, 100)
    M = symbol.evaluate(t1, t2, t3)
    print(M.shape)

def test_interface_3d_scalar_2(mapping=False):
    print('============ test_interface_3d_scalar_2 =============')

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

    symbol = DiscreteSymbol(a, Vh)

    t1 = linspace(0, 1, 100)
    t2 = linspace(0, 1, 100)
    t3 = linspace(0, 1, 100)
    M = symbol.evaluate(t1, t2, t3, c=0.5)
    print(M.shape)


#................................
if __name__ == '__main__':

    # .................................
    # without mapping
    test_interface_3d_scalar_1(mapping=False)
    test_interface_3d_scalar_2(mapping=False)
    # .................................
