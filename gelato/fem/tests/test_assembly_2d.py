# coding: utf-8

# TODO split the asserts between algebraic and weak formulations ones

from sympy.core.containers import Tuple
from sympy import symbols
from sympy import Symbol
from sympy import Function
from sympy import pi, cos

from gelato.core import dx, dy, dz
from gelato.core import Constant
from gelato.core import Field
from gelato.core import grad, dot, inner, cross, rot, curl, div
from gelato.core import H1Space
from gelato.core import TestFunction
from gelato.core import VectorTestFunction
from gelato.core import BilinearForm, LinearForm
from gelato.core import atomize, normalize
from gelato.core import gelatize

from gelato.fem  import discretize

from numpy import linspace

from spl.fem.splines import SplineSpace
from spl.fem.tensor  import TensorFemSpace

# ...
def test_assembly_bilinear_2d_scalar_1():
    print('============ test_assembly_bilinear_2d_scalar_1 =============')

    U = H1Space('U', ldim=2)
    V = H1Space('V', ldim=2)

    v = TestFunction(V, name='v')
    u = TestFunction(U, name='u')

    expr = dot(grad(v), grad(u))

    a = BilinearForm((v,u), expr)
    print('> input      >>> {0}'.format(a))

    # ... discretization
    # Input data: degree, number of elements
    p1  = 3  ; p2  = 3
    ne1 = 16 ; ne2 = 16

    # Create uniform grid
    grid_1 = linspace( 0., 1., num=ne1+1 )
    grid_2 = linspace( 0., 1., num=ne2+1 )

    # Create 1D finite element spaces and precompute quadrature data
    V1 = SplineSpace( p1, grid=grid_1 ); V1.init_fem()
    V2 = SplineSpace( p2, grid=grid_2 ); V2.init_fem()

    # Create 2D tensor product finite element space
    V = TensorFemSpace( V1, V2 )
    # ...

    # ...
    discretize( a, [V, V] )
    # ...

    # ...
    M = a.assemble()
    # ...
# ...

# ...
def test_assembly_linear_2d_scalar_1():
    print('============ test_assembly_linear_2d_scalar_1 =============')

    V = H1Space('V', ldim=2)

    v = TestFunction(V, name='v')

    x,y = V.coordinates

    #expr = cos(2*pi*x)*v
    expr = x*(1.-x)*y*(1.-y)*v

    a = LinearForm(v, expr)
    print('> input      >>> {0}'.format(a))

    # ... discretization
    # Input data: degree, number of elements
    p1  = 3  ; p2  = 3
    ne1 = 16 ; ne2 = 16

    # Create uniform grid
    grid_1 = linspace( 0., 1., num=ne1+1 )
    grid_2 = linspace( 0., 1., num=ne2+1 )

    # Create 1D finite element spaces and precompute quadrature data
    V1 = SplineSpace( p1, grid=grid_1 ); V1.init_fem()
    V2 = SplineSpace( p2, grid=grid_2 ); V2.init_fem()

    # Create 2D tensor product finite element space
    V = TensorFemSpace( V1, V2 )
    # ...

    # ...
    discretize( a, V )
    # ...

    # ...
    b = a.assemble()
    # ...
# ...


# .....................................................
if __name__ == '__main__':
    test_assembly_bilinear_2d_scalar_1()

    test_assembly_linear_2d_scalar_1()
