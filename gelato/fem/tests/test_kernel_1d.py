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
from gelato.core import BilinearForm, LinearForm, FunctionForm
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
def test_kernel_bilinear_1d_scalar_2():
    print('============ test_kernel_bilinear_1d_scalar_2 =============')

    U = H1Space('U', ldim=1)
    V = H1Space('V', ldim=1)

    v = TestFunction(V, name='v')
    u = TestFunction(U, name='u')

    c = Constant('c', real=True, label='mass stabilization')

    expr = dot(grad(v), grad(u)) + c*v*u

    a = BilinearForm((v,u), expr)
    print('> input      >>> {0}'.format(a))

    # ...
    kernel_py  = compile_kernel('kernel_bilinear_1d_scalar_2', a,
                                backend='python', verbose=True)
    # ...
# ...

# ...
def test_kernel_bilinear_1d_scalar_3():
    print('============ test_kernel_bilinear_1d_scalar_3 =============')

    U = H1Space('U', ldim=1)
    V = H1Space('V', ldim=1)

    v = TestFunction(V, name='v')
    u = TestFunction(U, name='u')

    F = Field('F', space=V)

    expr = dot(grad(v), grad(u)) + F*v*u

    a = BilinearForm((v,u), expr)
    print('> input      >>> {0}'.format(a))

    # ...
    kernel_py  = compile_kernel('kernel_bilinear_1d_scalar_3', a,
                                backend='python', verbose=True)
    # ...
# ...

# ...
def test_kernel_bilinear_1d_scalar_4():
    print('============ test_kernel_bilinear_1d_scalar_4 =============')

    U = H1Space('U', ldim=1)
    V = H1Space('V', ldim=1)

    v = TestFunction(V, name='v')
    u = TestFunction(U, name='u')

    F = Field('F', space=V)

    expr = dot(grad(F*v), grad(u))

    a = BilinearForm((v,u), expr)
    print('> input      >>> {0}'.format(a))

    # ...
    kernel_py  = compile_kernel('kernel_bilinear_1d_scalar_4', a,
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

# ...
def test_kernel_linear_1d_scalar_2():
    print('============ test_kernel_linear_1d_scalar_2 =============')

    V = H1Space('V', ldim=1)

    v = TestFunction(V, name='v')

    c = Constant('c', real=True, label='mass stabilization')

    x = V.coordinates

    expr = c*cos(2*pi*x)*v

    a = LinearForm(v, expr)
    print('> input      >>> {0}'.format(a))

    # ...
    kernel_py  = compile_kernel('kernel_linear_1d_scalar_2', a,
                                backend='python', verbose=True)
    # ...
# ...

# ...
def test_kernel_linear_1d_scalar_3():
    print('============ test_kernel_linear_1d_scalar_3 =============')

    V = H1Space('V', ldim=1)

    v = TestFunction(V, name='v')

    F = Field('F', space=V)

    x = V.coordinates

    expr = F*v

    a = LinearForm(v, expr)
    print('> input      >>> {0}'.format(a))

    # ...
    kernel_py  = compile_kernel('kernel_linear_1d_scalar_3', a,
                                backend='python', verbose=True)
    # ...
# ...

# ...
def test_kernel_linear_1d_scalar_4():
    print('============ test_kernel_linear_1d_scalar_4 =============')

    V = H1Space('V', ldim=1)

    v = TestFunction(V, name='v')

    F = Field('F', space=V)

    x = V.coordinates

    expr = dx(F)*v

    a = LinearForm(v, expr)
    print('> input      >>> {0}'.format(a))

    # ...
    kernel_py  = compile_kernel('kernel_linear_1d_scalar_4', a,
                                backend='python', verbose=True)
    # ...
# ...

# ...
def test_kernel_function_1d_scalar_1():
    print('============ test_kernel_function_1d_scalar_1 =============')

    V = H1Space('V', ldim=1)

    F = Field('F', space=V)

    x = V.coordinates

    expr = F-cos(2*pi*x)

    a = FunctionForm(expr)
    print('> input      >>> {0}'.format(a))

    # ...
    kernel_py  = compile_kernel('kernel_function_1d_scalar_1', a,
                                backend='python', verbose=True)
    # ...
# ...

# ...
def test_kernel_function_1d_scalar_2():
    print('============ test_kernel_function_1d_scalar_2 =============')

    V = H1Space('V', ldim=1)

    F = Field('F', space=V)

    x = V.coordinates

    expr = dot(grad(F-cos(2*pi*x)), grad(F-cos(2*pi*x)))

    a = FunctionForm(expr)
    print('> input      >>> {0}'.format(a))

    # ...
    kernel_py  = compile_kernel('kernel_function_1d_scalar_2', a,
                                backend='python', verbose=True)
    # ...
# ...


# .....................................................
if __name__ == '__main__':
    # ... scalar case
    test_kernel_bilinear_1d_scalar_1()
    test_kernel_bilinear_1d_scalar_2()
    test_kernel_bilinear_1d_scalar_3()
    test_kernel_bilinear_1d_scalar_4()

    test_kernel_linear_1d_scalar_1()
    test_kernel_linear_1d_scalar_2()
    test_kernel_linear_1d_scalar_3()
    test_kernel_linear_1d_scalar_4()

    test_kernel_function_1d_scalar_1()
    test_kernel_function_1d_scalar_2()
    # ...
