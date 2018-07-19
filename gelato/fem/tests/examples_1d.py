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

from numpy import linspace, array

from spl.fem.splines import SplineSpace
from spl.fem.basic   import FemField

import matplotlib.pyplot as plt

# ...
def test_poisson_1d_scalar_1():
    print('============ test_poisson_1d_scalar_1 =============')

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

    # ...
    M   = a.assemble()
    rhs = b.assemble()
    # ...

    # Apply homogeneous dirichlet boundary conditions
    from utils import apply_homogeneous_dirichlet_bc
    apply_homogeneous_dirichlet_bc(V, M, rhs)

    # Solve linear system
    from spl.linalg.solvers import cg
    x, info = cg( M, rhs, tol=1e-9, maxiter=1000, verbose=True )

    # ............................................
    #              PLOT
    # ............................................
    # Create potential field
    phi = FemField( V, 'phi' )
    phi.coeffs[:] = x[:]
    phi.coeffs.update_ghost_regions()

    # Plot solution on refined grid
    y      = linspace( grid[0], grid[-1], 101 )
    phi_y  = array( [phi(yj) for yj in y] )
    fig,ax = plt.subplots( 1, 1 )
    ax.plot( y, phi_y )
    ax.set_xlabel( 'x' )
    ax.set_ylabel( 'y' )
    ax.grid()

    # Show figure and keep it open if necessary
    fig.tight_layout()
    plt.show()
    # ............................................
# ...


# .....................................................
if __name__ == '__main__':
    # ... scalar case
    test_poisson_1d_scalar_1()
    # ...
