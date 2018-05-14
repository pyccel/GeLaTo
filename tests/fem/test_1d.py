# coding: utf-8

# TODO: - allow for giving a name for the trial/test basis
#       - Ni/Nj should be Ni_0/Nj_O

from numpy import linspace

from sympy.core.containers import Tuple
from sympy import symbols
from sympy import Symbol
from sympy import Lambda
from sympy import Function

from gelato.expression   import construct_weak_form
from gelato.calculus     import (Dot, Cross, Grad, Curl, Rot, Div, dx)
from gelato.calculus     import Constant
from gelato.calculus     import Field
from gelato.fem.assembly import assemble_matrix
from gelato.fem.utils    import compile_kernel

from spl.fem.splines import SplineSpace
from spl.fem.tensor  import TensorSpace
from spl.fem.vector  import VectorFemSpace
from spl.fem.splines import Spline

from utils import assert_identical_coo

# ...
def test_1d_scalar_1():
    print('============== test_1d_scalar_1 ================')

    # ... define the weak formulation
    x = Symbol('x')

    u = Symbol('u')
    v = Symbol('v')

    a = Lambda((x,v,u), Dot(Grad(u), Grad(v)) + u*v)
    # ...

    # ...  create a finite element space
    p  = 3
    ne = 64

    print('> Grid   :: {ne}'.format(ne=ne))
    print('> Degree :: {p}'.format(p=p))

    grid = linspace(0., 1., ne+1)

    V = SplineSpace(p, grid=grid)
    # ...

    # ...
    kernel_py  = compile_kernel('kernel_scalar_1', a, V, backend='python')
    kernel_f90 = compile_kernel('kernel_scalar_1', a, V, backend='fortran')

    M_py  = assemble_matrix(V, kernel_py)
    M_f90 = assemble_matrix(V, kernel_f90)
    # ...

    assert_identical_coo(M_py, M_f90)

# ...

# ...
def test_1d_scalar_2():
    print('============== test_1d_scalar_2 ================')

    # ... define the weak formulation
    x = Symbol('x')

    u = Symbol('u')
    v = Symbol('v')

    alpha = Constant('alpha')
    nu = Constant('nu')

    a = Lambda((x,v,u), alpha * Dot(Grad(u), Grad(v)) + nu*u*v)
    # ...

    # ...  create a finite element space
    p  = 3
    ne = 64

    print('> Grid   :: {ne}'.format(ne=ne))
    print('> Degree :: {p}'.format(p=p))

    grid = linspace(0., 1., ne+1)

    V = SplineSpace(p, grid=grid)
    # ...

    # ...
    kernel_py  = compile_kernel('kernel_scalar_2', a, V,
                                d_constants={'nu': 0.1},
                                d_args={'alpha': 'double'},
                                backend='python')
    kernel_f90 = compile_kernel('kernel_scalar_2', a, V,
                                d_constants={'nu': 0.1},
                                d_args={'alpha': 'double'},
                                backend='fortran')

    M_py  = assemble_matrix(V, kernel_py, args={'alpha': 2.0})
    M_f90 = assemble_matrix(V, kernel_f90, args={'alpha': 2.0})
    # ...

    assert_identical_coo(M_py, M_f90)

# ...

# ...
def test_1d_scalar_3():
    print('============== test_1d_scalar_3 ================')

    x = Symbol('x')

    u = Symbol('u')
    v = Symbol('v')

    b = Function('b')

    a = Lambda((x,v,u), Dot(Grad(u), Grad(v)) + b(x)*u*v)
    print('> input     := {0}'.format(a))

    # ...  create a finite element space
    p  = 3
    ne = 64

    print('> Grid   :: {ne}'.format(ne=ne))
    print('> Degree :: {p}'.format(p=p))

    grid = linspace(0., 1., ne+1)

    V = SplineSpace(p, grid=grid)
    # ...

    # ... user defined function
    def b(s):
        r = 1.+ s*(1.-s)
        return r
    # ...

    # ... create an interactive pyccel context
    from pyccel.epyccel import ContextPyccel

    context = ContextPyccel(name='context_3')
    context.insert_function(b, ['double'], kind='function', results=['double'])

    context.compile()
    # ...

    # ...
    kernel_py  = compile_kernel('kernel_scalar_3', a, V,
                                context=context,
                                verbose=True,
                                backend='python')

    kernel_f90 = compile_kernel('kernel_scalar_3', a, V,
                                context=context,
                                verbose=True,
                                backend='fortran')
    # ...

    # ...
    M_py  = assemble_matrix(V, kernel_py)
    M_f90 = assemble_matrix(V, kernel_f90)
    # ...

    assert_identical_coo(M_py, M_f90)

    print('')
# ...

# ...
def test_1d_scalar_4():
    print('============== test_1d_scalar_4 ================')

    # ... define the weak formulation
    x = Symbol('x')

    u = Symbol('u')
    v = Symbol('v')

    F = Field('F')

    a = Lambda((x,v,u), Dot(Grad(F*u), Grad(v)) + u*v)
    # ...

    # ...  create a finite element space
    p  = 3
    ne = 64

    print('> Grid   :: {ne}'.format(ne=ne))
    print('> Degree :: {p}'.format(p=p))

    grid = linspace(0., 1., ne+1)

    V = SplineSpace(p, grid=grid)
    # ...

    # ...
    kernel_py  = compile_kernel('kernel_scalar_4', a, V, backend='python')
    kernel_f90 = compile_kernel('kernel_scalar_4', a, V, backend='fortran')

    F = Spline(V)
    F.coeffs._data[:] = 1.

    M_py  = assemble_matrix(V, kernel_py, fields={'F': F})
    M_f90 = assemble_matrix(V, kernel_f90, fields={'F': F})
    # ...

    assert_identical_coo(M_py, M_f90)

# ...

# ...
def test_1d_scalar_5():
    print('============== test_1d_scalar_5 ================')

    # ... define the weak formulation
    x = Symbol('x')

    u = Symbol('u')
    v = Symbol('v')

    a = Lambda((x,v,u), dx(dx(u)) * dx(dx(v)))
    # ...

    # ...  create a finite element space
    p  = 3
    ne = 64

    print('> Grid   :: {ne}'.format(ne=ne))
    print('> Degree :: {p}'.format(p=p))

    grid = linspace(0., 1., ne+1)

    V = SplineSpace(p, grid=grid, nderiv=2)
    # ...

    # ...
    kernel_py  = compile_kernel('kernel_scalar_5', a, V, backend='python')
    kernel_f90 = compile_kernel('kernel_scalar_5', a, V, backend='fortran')

    M_py  = assemble_matrix(V, kernel_py)
    M_f90 = assemble_matrix(V, kernel_f90)
    # ...

    assert_identical_coo(M_py, M_f90)

# ...

# ...
def test_1d_block_1():
    print('============== test_1d_block_1 ================')

    # ... define the weak formulation
    x = Symbol('x')

    u0, u1 = symbols('u0 u1')
    v0, v1 = symbols('v0 v1')

    a = Lambda((x,v0,v1,u0,u1), dx(u0)*dx(v0) + dx(u1)*v0 + u0*dx(v1) + u1*v1)
    # ...

    # ...  create a finite element space
    p  = 3
    ne = 64

    print('> Grid   :: {ne}'.format(ne=ne))
    print('> Degree :: {p}'.format(p=p))

    grid = linspace(0., 1., ne+1)

    V = SplineSpace(p, grid=grid)
    # ...

    # ... Vector fem space
    V = VectorFemSpace(V, V)
    # ...

    # ...
    kernel_py  = compile_kernel('kernel_block_1', a, V, backend='python')
    kernel_f90 = compile_kernel('kernel_block_1', backend='fortran')

    M_py  = assemble_matrix(V, kernel_py)
    M_f90 = assemble_matrix(V, kernel_f90)
    # ...

    assert_identical_coo(M_py, M_f90)

# ...


# .....................................................
if __name__ == '__main__':

    test_1d_scalar_1()
    test_1d_scalar_2()
    test_1d_scalar_3()
    test_1d_scalar_4()
    test_1d_scalar_5()
#    test_1d_block_1()
