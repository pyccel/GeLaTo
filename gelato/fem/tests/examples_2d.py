# coding: utf-8

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
from gelato.core import BilinearForm, LinearForm, FunctionForm
from gelato.core import atomize, normalize
from gelato.core import gelatize

from gelato.fem  import discretize

from numpy import linspace, array, sqrt

from spl.fem.splines import SplineSpace
from spl.fem.tensor  import TensorFemSpace
from spl.fem.basic   import FemField

import matplotlib.pyplot as plt

# ...
def test_poisson_2d_scalar_1():
    print('============ test_poisson_2d_scalar_1 =============')

    # ... abstract model
    U = H1Space('U', ldim=2)
    V = H1Space('V', ldim=2)

    x,y = V.coordinates

    v = TestFunction(V, name='v')
    u = TestFunction(U, name='u')

    F = Field('F', space=V)

    a = BilinearForm((v,u), dot(grad(v), grad(u)))
    b = LinearForm(v, 4.*v)
    norm = FunctionForm((F-x*(1.-x)*y*(1.-y))**2)

    print('> input bilinear-form  >>> {0}'.format(a))
    print('> input linear-form    >>> {0}'.format(b))
    print('> input function-form  >>> {0}'.format(norm))
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

    # ...
    discretize( a, [V, V] )
    discretize( b, V )
    discretize( norm, V )
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
    phi.coeffs[:,:] = x[:,:]
    phi.coeffs.update_ghost_regions()

    # compute L2 error norm
    err = norm.assemble(phi)
    print('> L2 error =  ', sqrt(err))

    # Plot solution on refined grid
    xx = linspace( grid_1[0], grid_1[-1], num=101 )
    yy = linspace( grid_2[0], grid_2[-1], num=101 )
    zz = array( [[phi( xi,yi ) for yi in yy] for xi in xx] )
    fig, ax = plt.subplots( 1, 1 )
    im = ax.contourf( xx, yy, zz.transpose(), 40, cmap='jet' )
    fig.colorbar( im )
    ax.set_xlabel( r'$x$', rotation='horizontal' )
    ax.set_ylabel( r'$y$', rotation='horizontal' )
    ax.set_title ( r'$\phi(x,y)$' )
    ax.grid()

    # Show figure and keep it open if necessary
    fig.tight_layout()
    fig.show()
    plt.show()
    # ............................................
# ...

# .....................................................
if __name__ == '__main__':
    # ... scalar case
    test_poisson_2d_scalar_1()
    # ...
