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

from gelato.fem  import compile_assembly

from numpy import linspace

# ...
def test_assembly_bilinear_1d_scalar_1():
    print('============ test_assembly_bilinear_1d_scalar_1 =============')

    U = H1Space('U', ldim=1)
    V = H1Space('V', ldim=1)

    v = TestFunction(V, name='v')
    u = TestFunction(U, name='u')

    expr = dot(grad(v), grad(u))

    a = BilinearForm((v,u), expr)
    print('> input      >>> {0}'.format(a))

    # ...
    assembly_py  = compile_assembly('assembly_bilinear_1d_scalar_1', a,
                                    backend='python', verbose=True)
    # ...
# ...

# ...
def test_assembly_bilinear_1d_scalar_2():
    print('============ test_assembly_bilinear_1d_scalar_2 =============')

    U = H1Space('U', ldim=1)
    V = H1Space('V', ldim=1)

    v = TestFunction(V, name='v')
    u = TestFunction(U, name='u')

    c = Constant('c', real=True, label='mass stabilization')

    expr = dot(grad(v), grad(u)) + c*v*u

    a = BilinearForm((v,u), expr)
    print('> input      >>> {0}'.format(a))

    # ...
    assembly_py  = compile_assembly('assembly_bilinear_1d_scalar_2', a,
                                    backend='python', verbose=True)
    # ...
# ...

# ...
def test_assembly_bilinear_1d_scalar_3():
    print('============ test_assembly_bilinear_1d_scalar_3 =============')

    U = H1Space('U', ldim=1)
    V = H1Space('V', ldim=1)

    v = TestFunction(V, name='v')
    u = TestFunction(U, name='u')

    F = Field('F', space=V)

    expr = dot(grad(v), grad(u)) + F*v*u

    a = BilinearForm((v,u), expr)
    print('> input      >>> {0}'.format(a))

    # ...
    assembly_py  = compile_assembly('assembly_bilinear_1d_scalar_3', a,
                                    backend='python', verbose=True)
    # ...
# ...

# ...
def test_assembly_bilinear_1d_scalar_4():
    print('============ test_assembly_bilinear_1d_scalar_4 =============')

    U = H1Space('U', ldim=1)
    V = H1Space('V', ldim=1)

    v = TestFunction(V, name='v')
    u = TestFunction(U, name='u')

    F = Field('F', space=V)

    expr = dot(grad(F*v), grad(u))

    a = BilinearForm((v,u), expr)
    print('> input      >>> {0}'.format(a))

    # ...
    assembly_py  = compile_assembly('assembly_bilinear_1d_scalar_4', a,
                                    backend='python', verbose=True)
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

    # ...
    assembly_py  = compile_assembly('assembly_linear_1d_scalar_1', a,
                                    backend='python', verbose=True)
    # ...
# ...

# ...
def test_assembly_function_1d_scalar_1():
    print('============ test_assembly_function_1d_scalar_1 =============')

    V = H1Space('V', ldim=1)

    F = Field('F', space=V)

    x = V.coordinates

    expr = F-cos(2*pi*x)

    a = FunctionForm(expr)
    print('> input      >>> {0}'.format(a))

    # ...
    assembly_py  = compile_assembly('assembly_function_1d_scalar_1', a,
                                    backend='python', verbose=True)
    # ...
# ...


# .....................................................
if __name__ == '__main__':
    # ... scalar case
#    test_assembly_bilinear_1d_scalar_1()
#    test_assembly_bilinear_1d_scalar_2()
#    test_assembly_bilinear_1d_scalar_3()
#    test_assembly_bilinear_1d_scalar_4()
#
#    test_assembly_linear_1d_scalar_1()

    test_assembly_function_1d_scalar_1()
    # ...
