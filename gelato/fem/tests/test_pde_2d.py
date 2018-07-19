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
from spl.fem.basic   import FemField

# ...
def test_pde_2d_scalar_1():
    print('============ test_pde_2d_scalar_1 =============')

    # ... abstract model
    U = H1Space('U', ldim=2)
    V = H1Space('V', ldim=2)

    x,y = V.coordinates

    v = TestFunction(V, name='v')
    u = TestFunction(U, name='u')

    a = BilinearForm((v,u), dot(grad(v), grad(u)))
    b = LinearForm(v, x*(1.-x)*y*(1.-y)*v)

    print('> input bilinear-form  >>> {0}'.format(a))
    print('> input linear-form    >>> {0}'.format(b))
    # ...

    # ... discretization
    # Input data: degree, number of elements
    p1  = 3  ; p2  = 3
    ne1 = 4 ; ne2 = 4

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
    discretize( b, V )
    # ...

    # ...
    M   = a.assemble()
    rhs = b.assemble()
    # ...
# ...


# ...
def test_pde_2d_scalar_2():
    print('============ test_pde_2d_scalar_2 =============')

    # ... abstract model
    U = H1Space('U', ldim=2)
    V = H1Space('V', ldim=2)

    x,y = V.coordinates

    v = TestFunction(V, name='v')
    u = TestFunction(U, name='u')

    c = Constant('c', real=True, label='mass stabilization')

    a = BilinearForm((v,u), dot(grad(v), grad(u)) + c*v*u)
    b = LinearForm(v, x*(1.-x)*y*(1.-y)*v)

    print('> input bilinear-form  >>> {0}'.format(a))
    print('> input linear-form    >>> {0}'.format(b))
    # ...

    # ... discretization
    # Input data: degree, number of elements
    p1  = 3  ; p2  = 3
    ne1 = 4 ; ne2 = 4

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
    discretize( b, V )
    # ...

    # ...
    M   = a.assemble(0.1)
    rhs = b.assemble()
    # ...
# ...

# ...
def test_pde_2d_scalar_3():
    print('============ test_pde_2d_scalar_3 =============')

    # ... abstract model
    U = H1Space('U', ldim=2)
    V = H1Space('V', ldim=2)

    x,y = V.coordinates

    v = TestFunction(V, name='v')
    u = TestFunction(U, name='u')

    F = Field('F', space=V)

    a = BilinearForm((v,u), dot(grad(v), grad(u)) + F*v*u)
    b = LinearForm(v, x*(1.-x)*y*(1.-y)*v)

    print('> input bilinear-form  >>> {0}'.format(a))
    print('> input linear-form    >>> {0}'.format(b))
    # ...

    # ... discretization
    # Input data: degree, number of elements
    p1  = 3  ; p2  = 3
    ne1 = 4 ; ne2 = 4

    # Create uniform grid
    grid_1 = linspace( 0., 1., num=ne1+1 )
    grid_2 = linspace( 0., 1., num=ne2+1 )

    # Create 1D finite element spaces and precompute quadrature data
    V1 = SplineSpace( p1, grid=grid_1 ); V1.init_fem()
    V2 = SplineSpace( p2, grid=grid_2 ); V2.init_fem()

    # Create 2D tensor product finite element space
    V = TensorFemSpace( V1, V2 )

    # Define a field
    phi = FemField( V, 'phi' )
    phi._coeffs[:,:] = 1.
    # ...

    # ...
    discretize( a, [V, V] )
    discretize( b, V )
    # ...

    # ...
    M   = a.assemble(phi)
    rhs = b.assemble()
    # ...
# ...

# ...
def test_pde_2d_scalar_4():
    print('============ test_pde_2d_scalar_4 =============')

    # ... abstract model
    U = H1Space('U', ldim=2)
    V = H1Space('V', ldim=2)

    x,y = V.coordinates

    v = TestFunction(V, name='v')
    u = TestFunction(U, name='u')

    F = Field('F', space=V)

    a = BilinearForm((v,u), dot(grad(F*v), grad(u)))
    b = LinearForm(v, x*(1.-x)*y*(1.-y)*v)

    print('> input bilinear-form  >>> {0}'.format(a))
    print('> input linear-form    >>> {0}'.format(b))
    # ...

    # ... discretization
    # Input data: degree, number of elements
    p1  = 3  ; p2  = 3
    ne1 = 4 ; ne2 = 4

    # Create uniform grid
    grid_1 = linspace( 0., 1., num=ne1+1 )
    grid_2 = linspace( 0., 1., num=ne2+1 )

    # Create 1D finite element spaces and precompute quadrature data
    V1 = SplineSpace( p1, grid=grid_1 ); V1.init_fem()
    V2 = SplineSpace( p2, grid=grid_2 ); V2.init_fem()

    # Create 2D tensor product finite element space
    V = TensorFemSpace( V1, V2 )

    # Define a field
    phi = FemField( V, 'phi' )
    phi._coeffs[:,:] = 1.
    # ...

    # ...
    discretize( a, [V, V] )
    discretize( b, V )
    # ...

    # ...
    M   = a.assemble(phi)
    rhs = b.assemble()
    # ...
# ...

# ...
def test_pde_2d_block_1():
    print('============ test_pde_2d_block_1 =============')

    # ... abstract model
    V = H1Space('V', ldim=2, is_block=True, shape=2)
    U = H1Space('U', ldim=2, is_block=True, shape=2)

    x,y = V.coordinates

    v = VectorTestFunction(V, name='v')
    u = VectorTestFunction(U, name='u')

    a = BilinearForm((v,u), div(v) * div(u) + rot(v) * rot(u))
    b = LinearForm(v, x*(1.-y)*v[0] + (1.-x)*y*v[1])

    print('> input bilinear-form  >>> {0}'.format(a))
    print('> input linear-form    >>> {0}'.format(b))
    # ...

    # ... discretization
    # Input data: degree, number of elements
    p1  = 3  ; p2  = 3
    ne1 = 4 ; ne2 = 4

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
    discretize( b, V )
    # ...

    # ...
    M   = a.assemble()
    rhs = b.assemble()
    # ...
# ...

# .....................................................
if __name__ == '__main__':
    # ... scalar case
    test_pde_2d_scalar_1()
    test_pde_2d_scalar_2()
    test_pde_2d_scalar_3()
    test_pde_2d_scalar_4()
    # ...

    # ... block case
    test_pde_2d_block_1()
    # ...
