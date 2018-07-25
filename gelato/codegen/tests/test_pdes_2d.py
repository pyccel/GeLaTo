# coding: utf-8

from numpy import linspace

from sympy import Symbol
from sympy.core.containers import Tuple
from sympy import symbols
from sympy import pi, cos, sin
from sympy import srepr

from symfe.core import dx, dy, dz
from symfe.core import Constant
from symfe.core import Field
from symfe.core import grad, dot, inner, cross, rot, curl, div
from symfe.core import H1Space
from symfe.core import TestFunction
from symfe.core import VectorTestFunction
from symfe.core import BilinearForm

from gelato.codegen import discretize

from spl.fem.splines import SplineSpace
from spl.fem.tensor  import TensorFemSpace

# ...
def test_pdes_2d_1():
    print('============ test_pdes_2d_1 =============')

    # ... abstract model
    V = H1Space('V', ldim=2)
    V_0 = H1Space('V_0', ldim=1, coordinates=['x'])
    V_1 = H1Space('V_1', ldim=1, coordinates=['y'])

    v = TestFunction(V, name='v')
    u = TestFunction(V, name='u')

    a = BilinearForm((v,u), u*v)
    # ...

    # ... discretization
    # Input data: degree, number of elements
    p1  = 3 ; p2  = 3
    ne1 = 2**4 ; ne2 = 2**4

    # Create uniform grid
    grid_1 = linspace( 0., 1., num=ne1+1 )
    grid_2 = linspace( 0., 1., num=ne2+1 )

    # Create 1D finite element spaces and precompute quadrature data
    V1 = SplineSpace( p1, grid=grid_1 ); V1.init_fem()
    V2 = SplineSpace( p2, grid=grid_2 ); V2.init_fem()

    # Create 2D tensor product finite element space
    V = TensorFemSpace( V1, V2 )
    # ...

    discretize( a, [p1, p2] )

# ...

# .....................................................
if __name__ == '__main__':
    test_pdes_2d_1()
