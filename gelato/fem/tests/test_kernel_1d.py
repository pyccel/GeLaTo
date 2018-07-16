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

from gelato.fem    import compile_kernel

from numpy import linspace

# ...
def test_kernel_bilinear_1d_scalar_1():
    print('============ test_kernel_bilinear_1d_scalar_1 =============')

    U = H1Space('U', ldim=1)
    V = H1Space('V', ldim=1)

    v = TestFunction(V, name='v')
    u = TestFunction(U, name='u')

    expr = inner(grad(v), grad(u))

    a = BilinearForm((v,u), expr)
    print('> input      >>> {0}'.format(a))

    # ...
    kernel_py  = compile_kernel('kernel_bilinear_1d_scalar_1', a,
                                backend='python', verbose=True)
    # ...
# ...

# ...
def test_kernel_linear_1d_scalar_1():
    print('============ test_kernel_linear_1d_scalar_1 =============')

    V = H1Space('V', ldim=1)

    v = TestFunction(V, name='v')

    x = V.coordinates

    expr = cos(2*pi*x)*v

    a = LinearForm(v, expr)
    print('> input      >>> {0}'.format(a))

    # ...
    kernel_py  = compile_kernel('kernel_linear_1d_scalar_1', a,
                                backend='python', verbose=True)
    # ...
# ...


# .....................................................
if __name__ == '__main__':
    test_kernel_bilinear_1d_scalar_1()

    test_kernel_linear_1d_scalar_1()
