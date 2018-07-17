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
from spl.fem.basic   import FemField

# ...
def test_pde_1d_scalar_1():
    print('============ test_pde_1d_scalar_1 =============')

    # ... abstract model
    U = H1Space('U', ldim=1)
    V = H1Space('V', ldim=1)

    x = V.coordinates

    v = TestFunction(V, name='v')
    u = TestFunction(U, name='u')

    a = BilinearForm((v,u), dot(grad(v), grad(u)))
    b = LinearForm(v, x*(1.-x)*v)

    print('> input bilinear-form  >>> {0}'.format(a))
    print('> input linear-form    >>> {0}'.format(b))
    # ...

    # ... discretization
    # Input data: degree, number of elements
    p  = 3
    ne = 2**4

    # Create uniform grid
    grid = linspace( 0., 1., num=ne+1 )

    # Create finite element space and precompute quadrature data
    V = SplineSpace( p, grid=grid )
    V.init_fem()
    # ...

    # ...
    discretize( a, [V, V] )
    discretize( b, V )
    # ...

#    # ...
#    print(a.assemble.__doc__)
#    print(b.assemble.__doc__)
#    # ...

    # ...
    M   = a.assemble()
    rhs = b.assemble()
    # ...
# ...

# ...
def test_pde_1d_scalar_2():
    print('============ test_pde_1d_scalar_2 =============')

    # ... abstract model
    U = H1Space('U', ldim=1)
    V = H1Space('V', ldim=1)

    x = V.coordinates

    v = TestFunction(V, name='v')
    u = TestFunction(U, name='u')

    c = Constant('c', real=True, label='mass stabilization')

    a = BilinearForm((v,u), dot(grad(v), grad(u)) + c*v*u)
    b = LinearForm(v, x*(1.-x)*v)

    print('> input bilinear-form  >>> {0}'.format(a))
    print('> input linear-form    >>> {0}'.format(b))
    # ...

    # ... discretization
    # Input data: degree, number of elements
    p  = 3
    ne = 2**4

    # Create uniform grid
    grid = linspace( 0., 1., num=ne+1 )

    # Create finite element space and precompute quadrature data
    V = SplineSpace( p, grid=grid )
    V.init_fem()
    # ...

    # ...
    discretize( a, [V, V] )
    discretize( b, V )
    # ...

#    # ...
#    print(a.assemble.__doc__)
#    print(b.assemble.__doc__)
#    # ...

    # ...
    M   = a.assemble(0.1)
    rhs = b.assemble()
    # ...
# ...


# ...
def test_pde_1d_scalar_3():
    print('============ test_pde_1d_scalar_3 =============')

    # ... abstract model
    U = H1Space('U', ldim=1)
    V = H1Space('V', ldim=1)

    x = V.coordinates

    v = TestFunction(V, name='v')
    u = TestFunction(U, name='u')

    F = Field('F', space=V)

    a = BilinearForm((v,u), dot(grad(v), grad(u)) + F*v*u)
    b = LinearForm(v, x*(1.-x)*v)

    print('> input bilinear-form  >>> {0}'.format(a))
    print('> input linear-form    >>> {0}'.format(b))
    # ...

    # ... discretization
    # Input data: degree, number of elements
    p  = 3
    ne = 2**4

    # Create uniform grid
    grid = linspace( 0., 1., num=ne+1 )

    # Create finite element space and precompute quadrature data
    V = SplineSpace( p, grid=grid )
    V.init_fem()

    # Define a field
    phi = FemField( V, 'phi' )
    phi._coeffs[:] = 1.
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

# .....................................................
if __name__ == '__main__':
    test_pde_1d_scalar_1()
    test_pde_1d_scalar_2()
    test_pde_1d_scalar_3()
