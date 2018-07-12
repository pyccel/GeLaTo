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

from gelato.fem.core import BasicSobolevSpace
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
def test_kernel_3d_scalar_1():
    print('============ test_kernel_3d_scalar_1 =============')

    U = BasicSobolevSpace('U', ldim=3)
    V = BasicSobolevSpace('V', ldim=3)

    v = TestFunction(V, name='v')
    u = TestFunction(U, name='u')

    expr = inner(grad(v), grad(u))

    a = BilinearForm((v,u), expr)
    print('> input      >>> {0}'.format(a))

    # ...  create a finite element space
    p1  = 2 ; p2  = 2 ; p3  = 2
    ne1 = 2 ; ne2 = 2 ; ne3 = 2

    print('> Grid   :: [{},{},{}]'.format(ne1, ne2, ne3))
    print('> Degree :: [{},{},{}]'.format(p1, p2, p3))

    grid_1 = linspace(0., 1., ne1+1)
    grid_2 = linspace(0., 1., ne2+1)
    grid_3 = linspace(0., 1., ne3+1)

    V1 = SplineSpace(p1, grid=grid_1)
    V2 = SplineSpace(p2, grid=grid_2)
    V3 = SplineSpace(p3, grid=grid_3)

    V = TensorFemSpace(V1, V2, V3)
    # ...

    # ...
    kernel_py  = compile_kernel('kernel_3d_scalar_1', a, [V, V],
                                backend='python', verbose=True)
    # ...

# ...


# .....................................................
if __name__ == '__main__':
    test_kernel_3d_scalar_1()

