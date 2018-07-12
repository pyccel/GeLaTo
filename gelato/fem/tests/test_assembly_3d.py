# coding: utf-8

# TODO - split the asserts between algebraic and weak formulations ones
#      - add assert for grad in vector case

from sympy import Symbol
from sympy.core.containers import Tuple
from sympy import symbols
from sympy import IndexedBase
from sympy import Matrix

from gelato.calculus import dx, dy, dz
from gelato.calculus import Constant
from gelato.calculus import Field
from gelato.calculus import grad, dot, inner, cross, rot, curl, div

from gelato.fem.core import FemSpace
from gelato.fem.core import TestFunction
from gelato.fem.core import VectorTestFunction
from gelato.fem.expr import BilinearForm, LinearForm
from gelato.fem.expr import atomize, normalize, matricize
from gelato.fem.expr import gelatize

from gelato.fem.utils    import compile_kernel

from numpy import linspace

from spl.fem.splines import SplineSpace


# ...
def test_kernel_3d_scalar_1():
    print('============ test_kernel_3d_scalar_1 =============')

    U = FemSpace('U', ldim=3)
    V = FemSpace('V', ldim=3)

    v = TestFunction(V, name='v')
    u = TestFunction(U, name='u')

    expr = inner(grad(v), grad(u))

    a = BilinearForm((v,u), expr)
    print('> input      >>> {0}'.format(a))

    a_expr = gelatize(a, basis={V: 'Nj', U: 'Ni'})
    print('> gelatized  >>> {0}'.format(a_expr))
    print('')

    # ...  create a finite element space
    p  = 3
    ne = 64

    print('> Grid   :: {ne}'.format(ne=ne))
    print('> Degree :: {p}'.format(p=p))

    grid = linspace(0., 1., ne+1)

    V = SplineSpace(p, grid=grid)
    # ...

    # ...
    kernel_py  = compile_kernel('kernel_3d_scalar_1', a, [V, V],
                                backend='python', verbose=True)
    # ...

# ...


# .....................................................
if __name__ == '__main__':
    test_kernel_3d_scalar_1()

