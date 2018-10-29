# coding: utf-8

from numpy import linspace, zeros, pi

from sympy import Symbol
from sympy.core.containers import Tuple
from sympy import symbols
from sympy import srepr

from sympde.core import dx, dy, dz
from sympde.core import Constant
from sympde.core import Field
from sympde.core import grad, dot, inner, cross, rot, curl, div
from sympde.core import H1Space
from sympde.core import TestFunction
from sympde.core import VectorTestFunction
from sympde.core import BilinearForm

from gelato.codegen import compile_symbol
from gelato.codegen import discretize_symbol

from spl.fem.splines import SplineSpace
from spl.fem.tensor  import TensorFemSpace

# ...
def test_pdes_3d_1():
    print('============ test_pdes_3d_1 =============')

    # ... abstract model
    V = H1Space('V', ldim=3)

    v = TestFunction(V, name='v')
    u = TestFunction(V, name='u')

    a = BilinearForm((v,u), dot(grad(v), grad(u)))
    # ...

    # ... discretization
    # Input data: degree, number of elements
    p1  = 1 ; p2  = 1 ; p3  = 1
    ne1 = 4 ; ne2 = 4 ; ne3 = 4

    # Create uniform grid
    grid_1 = linspace( 0., 1., num=ne1+1 )
    grid_2 = linspace( 0., 1., num=ne2+1 )
    grid_3 = linspace( 0., 1., num=ne3+1 )

    # Create 1D finite element spaces and precompute quadrature data
    V1 = SplineSpace( p1, grid=grid_1 ); V1.init_fem()
    V2 = SplineSpace( p2, grid=grid_2 ); V2.init_fem()
    V3 = SplineSpace( p3, grid=grid_3 ); V3.init_fem()

    # Create 2D tensor product finite element space
    V = TensorFemSpace( V1, V2, V3 )
    # ...

    # ...
    discretize_symbol( a, [V,V] )
#    print(mass.symbol.__doc__)
    # ...

    # ...
    n1 = 21 ; n2 = 21 ; n3 = 21

    t1 = linspace(-pi, pi, n1)
    t2 = linspace(-pi, pi, n2)
    t3 = linspace(-pi, pi, n3)
    x1 = linspace(0.,1., n1)
    x2 = linspace(0.,1., n2)
    x3 = linspace(0.,1., n3)

    xs = [x1, x2, x3]
    ts = [t1, t2, t3]

    e = a.symbol(*xs, *ts)
    print('> ', e.min(), e.max())
    # ...
# ...

# ...
def test_pdes_3d_2():
    print('============ test_pdes_3d_2 =============')

    # ... abstract model
    V = H1Space('V', ldim=3)

    v = TestFunction(V, name='v')
    u = TestFunction(V, name='u')

    c = Constant('c', real=True, label='mass stabilization')

    a = BilinearForm((v,u), dot(grad(v), grad(u)) + c*v*u)
    # ...

    # ... discretization
    # Input data: degree, number of elements
    p1  = 1 ; p2  = 1 ; p3  = 1
    ne1 = 4 ; ne2 = 4 ; ne3 = 4

    # Create uniform grid
    grid_1 = linspace( 0., 1., num=ne1+1 )
    grid_2 = linspace( 0., 1., num=ne2+1 )
    grid_3 = linspace( 0., 1., num=ne3+1 )

    # Create 1D finite element spaces and precompute quadrature data
    V1 = SplineSpace( p1, grid=grid_1 ); V1.init_fem()
    V2 = SplineSpace( p2, grid=grid_2 ); V2.init_fem()
    V3 = SplineSpace( p3, grid=grid_3 ); V3.init_fem()

    # Create 2D tensor product finite element space
    V = TensorFemSpace( V1, V2, V3 )
    # ...

    # ...
    discretize_symbol( a, [V,V] )
#    print(mass.symbol.__doc__)
    # ...

    # ...
    n1 = 21 ; n2 = 21 ; n3 = 21

    t1 = linspace(-pi, pi, n1)
    t2 = linspace(-pi, pi, n2)
    t3 = linspace(-pi, pi, n3)
    x1 = linspace(0.,1., n1)
    x2 = linspace(0.,1., n2)
    x3 = linspace(0.,1., n3)

    xs = [x1, x2, x3]
    ts = [t1, t2, t3]

    e = a.symbol(*xs, *ts, 0.25)
    print('> c = 0.25 :: ', e.min(), e.max())

    e = a.symbol(*xs, *ts, 0.6)
    print('> c = 0.6 :: ', e.min(), e.max())
    # ...
# ...

# .....................................................
if __name__ == '__main__':
    test_pdes_3d_1()
    test_pdes_3d_2()
