# coding: utf-8

# TODO split the asserts between algebraic and weak formulations ones

from sympy.core.containers import Tuple
from sympy import symbols
from sympy import Symbol
from sympy import Function
from sympy import pi, cos

from gelato.calculus import dx, dy, dz
from gelato.calculus import Constant
from gelato.calculus import Field
from gelato.calculus import grad, dot, inner, cross, rot, curl, div

from gelato.fem.core import H1Space
from gelato.fem.core import TestFunction
from gelato.fem.core import VectorTestFunction
from gelato.fem.expr import BilinearForm, LinearForm
from gelato.fem.expr import atomize, normalize
from gelato.fem.expr import gelatize

from gelato.fem.utils    import compile_kernel

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

    x = Symbol('x')

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
#    test_kernel_bilinear_1d_scalar_1()

    test_kernel_linear_1d_scalar_1()
