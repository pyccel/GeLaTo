# coding: utf-8

# TODO - split the asserts between algebraic and weak formulations ones
#      - add assert for grad in vector case

from sympy import Symbol
from sympy.core.containers import Tuple
from sympy import symbols
from sympy import IndexedBase
from sympy import Matrix
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
from gelato.core import atomize, normalize, matricize
from gelato.core import gelatize

from gelato.fem.utils    import compile_kernel

from numpy import linspace

# ...
def test_kernel_bilinear_3d_scalar_1():
    print('============ test_kernel_bilinear_3d_scalar_1 =============')

    U = H1Space('U', ldim=3)
    V = H1Space('V', ldim=3)

    v = TestFunction(V, name='v')
    u = TestFunction(U, name='u')

    expr = inner(grad(v), grad(u))

    a = BilinearForm((v,u), expr)
    print('> input      >>> {0}'.format(a))

    # ...
    kernel_py  = compile_kernel('kernel_bilinear_3d_scalar_1', a,
                                backend='python', verbose=True)
    # ...
# ...

# ...
def test_kernel_bilinear_3d_block_1():
    print('============ test_kernel_bilinear_3d_block_1 =============')

    V = H1Space('V', ldim=3, is_block=True, shape=3)
    U = H1Space('U', ldim=3, is_block=True, shape=3)

    v = VectorTestFunction(V, name='v')
    u = VectorTestFunction(U, name='u')

    expr = div(u) * div(v) + 0.2 * dot(u, v)

    a = BilinearForm((v,u), expr)
    print('> input      >>> {0}'.format(a))

    # ...
    kernel_py  = compile_kernel('kernel_bilinear_3d_block_1', a,
                                backend='python', verbose=True)
    # ...
# ...

# ...
def test_kernel_linear_3d_scalar_1():
    print('============ test_kernel_linear_3d_scalar_1 =============')

    V = H1Space('V', ldim=3)

    v = TestFunction(V, name='v')

    x,y,z = V.coordinates

    expr = cos(2*pi*x)*cos(4*pi*y)*cos(5*pi*z)*v

    a = LinearForm(v, expr)
    print('> input      >>> {0}'.format(a))

    # ...
    kernel_py  = compile_kernel('kernel_linear_3d_scalar_1', a,
                                backend='python', verbose=True)
    # ...
# ...

# ...
def test_kernel_linear_3d_block_1():
    print('============ test_kernel_linear_3d_block_1 =============')

    V = H1Space('V', ldim=3, is_block=True, shape=3)

    v = VectorTestFunction(V, name='v')

    x,y,z = V.coordinates

    b = Tuple(2, 3, 5)
    expr = v[0] - v[1] + v[2] #dot(v, b)

    a = LinearForm(v, expr)
    print('> input      >>> {0}'.format(a))

    # ...
    kernel_py  = compile_kernel('kernel_linear_3d_block_1', a,
                                backend='python', verbose=True)
    # ...
# ...

# .....................................................
if __name__ == '__main__':
    test_kernel_bilinear_3d_scalar_1()
    test_kernel_bilinear_3d_block_1()

    test_kernel_linear_3d_scalar_1()
    test_kernel_linear_3d_block_1()
