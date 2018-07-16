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

# ...
def test_assembly_bilinear_1d_scalar_1():
    print('============ test_assembly_bilinear_1d_scalar_1 =============')

    U = H1Space('U', ldim=1)
    V = H1Space('V', ldim=1)

    v = TestFunction(V, name='v')
    u = TestFunction(U, name='u')

    expr = inner(grad(v), grad(u))

    a = BilinearForm((v,u), expr)
    print('> input      >>> {0}'.format(a))

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
    discretize(a,
               spaces=[V, V],
               backend='python', verbose=True)
    # ...

    # ...
    M = a.assemble()
    print(M.tocoo())
    # ...
# ...

# ...
def test_assembly_linear_1d_scalar_1():
    print('============ test_assembly_linear_1d_scalar_1 =============')

    V = H1Space('V', ldim=1)

    v = TestFunction(V, name='v')

    x = V.coordinates

    #expr = cos(2*pi*x)*v
    expr = x*(1.-x)*v

    a = LinearForm(v, expr)
    print('> input      >>> {0}'.format(a))

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
    discretize(a,
               spaces=V,
               backend='python', verbose=True)
    # ...

    # ...
    b = a.assemble()
    print(b.toarray())
    # ...
# ...


# .....................................................
if __name__ == '__main__':
    test_assembly_bilinear_1d_scalar_1()

    test_assembly_linear_1d_scalar_1()
