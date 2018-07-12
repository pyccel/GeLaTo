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

from gelato.fem.core import H1Space
from gelato.fem.core import TestFunction
from gelato.fem.core import VectorTestFunction
from gelato.fem.expr import BilinearForm, LinearForm
from gelato.fem.expr import atomize, normalize, matricize
from gelato.fem.expr import gelatize

from gelato.fem.utils    import compile_kernel

from numpy import linspace

from spl.fem.splines import SplineSpace
from spl.fem.tensor  import TensorFemSpace


# ...
def test_kernel_2d_scalar_1():
    print('============ test_kernel_2d_scalar_1 =============')

    U = H1Space('U', ldim=2)
    V = H1Space('V', ldim=2)

    v = TestFunction(V, name='v')
    u = TestFunction(U, name='u')

    expr = inner(grad(v), grad(u))

    a = BilinearForm((v,u), expr)
    print('> input      >>> {0}'.format(a))

    # ...  create a finite element space
    p1  = 2 ; p2  = 2
    ne1 = 8 ; ne2 = 8

    print('> Grid   :: [{ne1},{ne2}]'.format(ne1=ne1, ne2=ne2))
    print('> Degree :: [{p1},{p2}]'.format(p1=p1, p2=p2))

    grid_1 = linspace(0., 1., ne1+1)
    grid_2 = linspace(0., 1., ne2+1)

    V1 = SplineSpace(p1, grid=grid_1)
    V2 = SplineSpace(p2, grid=grid_2)

    V = TensorFemSpace(V1, V2)
    # ...

    # ...
    kernel_py  = compile_kernel('kernel_2d_scalar_1', a, [V, V],
                                backend='python', verbose=True)
    # ...

# ...


# .....................................................
if __name__ == '__main__':
    test_kernel_2d_scalar_1()

