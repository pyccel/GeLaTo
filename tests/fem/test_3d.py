# coding: utf-8

# TODO: - allow for giving a name for the trial/test basis
#       - Ni/Nj should be Ni_0/Nj_O

from numpy import linspace

from sympy.core.containers import Tuple
from sympy import symbols
from sympy import Symbol
from sympy import Lambda
from sympy import Function
from sympy import IndexedBase

from gelato.expression   import construct_weak_form
from gelato.calculus     import (Dot, Cross, Grad, Curl, Rot, Div)
from gelato.calculus     import Constant
from gelato.calculus     import Field
from gelato.fem.assembly import assemble_matrix
from gelato.fem.utils    import compile_kernel
from gelato.fem.utils    import compile_eval_field

from spl.fem.splines import SplineSpace
from spl.fem.tensor  import TensorSpace
from spl.fem.vector  import VectorFemSpace

from utils import assert_identical_coo

# ...
def test_3d_1():
    # ... define the weak formulation
    x,y,z = symbols('x y z')

    u = Symbol('u')
    v = Symbol('v')

    a = Lambda((x,y,z,v,u), Dot(Grad(u), Grad(v)) + u*v)
    # ...

    # ...  create a finite element space
    p1  = 2 ; p2  = 2 ; p3  = 2
    ne1 = 2 ; ne2 = 2 ; ne3 = 2
    # ...

    print('> Grid   :: [{},{},{}]'.format(ne1, ne2, ne3))
    print('> Degree :: [{},{},{}]'.format(p1, p2, p3))

    grid_1 = linspace(0., 1., ne1+1)
    grid_2 = linspace(0., 1., ne2+1)
    grid_3 = linspace(0., 1., ne3+1)

    V1 = SplineSpace(p1, grid=grid_1)
    V2 = SplineSpace(p2, grid=grid_2)
    V3 = SplineSpace(p3, grid=grid_3)

    V = TensorSpace(V1, V2, V3)
    # ...

    # ...
    kernel_py  = compile_kernel('kernel_1', a, V, backend='python')
    kernel_f90 = compile_kernel('kernel_1', a, V, backend='fortran')

    M_py  = assemble_matrix(V, kernel_py)
    M_f90 = assemble_matrix(V, kernel_f90)
    # ...

    assert_identical_coo(M_py, M_f90)

# ...

# ...
def test_3d_2():
    # ... define the weak formulation
    x,y,z = symbols('x y z')

    u = Symbol('u')
    v = Symbol('v')

    alpha = Constant('alpha')
    nu = Constant('nu')

    a = Lambda((x,v,u), alpha * Dot(Grad(u), Grad(v)) + nu*u*v)
    # ...

    # ...  create a finite element space
    p1  = 2 ; p2  = 2 ; p3  = 2
    ne1 = 2 ; ne2 = 2 ; ne3 = 2
    # ...

    print('> Grid   :: [{},{},{}]'.format(ne1, ne2, ne3))
    print('> Degree :: [{},{},{}]'.format(p1, p2, p3))

    grid_1 = linspace(0., 1., ne1+1)
    grid_2 = linspace(0., 1., ne2+1)
    grid_3 = linspace(0., 1., ne3+1)

    V1 = SplineSpace(p1, grid=grid_1)
    V2 = SplineSpace(p2, grid=grid_2)
    V3 = SplineSpace(p3, grid=grid_3)

    V = TensorSpace(V1, V2, V3)
    # ...

    # ...
    kernel_py  = compile_kernel('kernel_2', a, V,
                                d_constants={'nu': 0.1},
                                d_args={'alpha': 'double'},
                                backend='python')
    kernel_f90 = compile_kernel('kernel_2', a, V,
                                d_constants={'nu': 0.1},
                                d_args={'alpha': 'double'},
                                backend='fortran')

    M_py  = assemble_matrix(V, kernel_py, args={'alpha': 2.0})
    M_f90 = assemble_matrix(V, kernel_f90, args={'alpha': 2.0})
    # ...

    assert_identical_coo(M_py, M_f90)

# ...

# ...
def test_3d_3():
    x,y,z = symbols('x y z')

    u = IndexedBase('u')
    v = IndexedBase('v')

    a = Lambda((x,y,z,v,u), Dot(Curl(u), Curl(v)) + 0.2 * Dot(u, v))

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

    W = TensorSpace(V1, V2, V3)
    # ...

    # ... vector space
    V = VectorFemSpace(W, W, W)
    # ...

    # ...
    kernel_py  = compile_kernel('kernel_3', a, V, backend='python')
    kernel_f90 = compile_kernel('kernel_3', a, V, backend='fortran')

    M_py  = assemble_matrix(V, kernel_py)
    M_f90 = assemble_matrix(V, kernel_f90)
    # ...

    assert_identical_coo(M_py, M_f90)

# ...

# ...
def test_3d_4():
    # ... define the weak formulation
    x,y,z = symbols('x y z')

    u = Symbol('u')
    v = Symbol('v')

    b = Function('b')

    a = Lambda((x,y,z,v,u), Dot(Grad(u), Grad(v)) + b(x,y,z)*u*v)
    # ...

    # ...  create a finite element space
    p1  = 2 ; p2  = 2 ; p3  = 2
    ne1 = 2 ; ne2 = 2 ; ne3 = 2
    # ...

    print('> Grid   :: [{},{},{}]'.format(ne1, ne2, ne3))
    print('> Degree :: [{},{},{}]'.format(p1, p2, p3))

    grid_1 = linspace(0., 1., ne1+1)
    grid_2 = linspace(0., 1., ne2+1)
    grid_3 = linspace(0., 1., ne3+1)

    V1 = SplineSpace(p1, grid=grid_1)
    V2 = SplineSpace(p2, grid=grid_2)
    V3 = SplineSpace(p3, grid=grid_3)

    V = TensorSpace(V1, V2, V3)
    # ...

    # ... user defined function
    def b(x,y,z):
        r = 1.+ x*(1.-x) + y*(1.-y) + z*(1.-z)
        return r
    # ...

    # ... create an interactive pyccel context
    from pyccel.epyccel import ContextPyccel

    context = ContextPyccel(name='context_4')
    context.insert_function(b, ['double', 'double', 'double'], kind='function', results=['double'])

    context.compile()
    # ...

    # ...
    kernel_py  = compile_kernel('kernel_4', a, V,
                                context=context,
                                verbose=True,
                                backend='python')

    kernel_f90 = compile_kernel('kernel_4', a, V,
                                context=context,
                                verbose=True,
                                backend='fortran')
    # ...

    # ...
    M_py  = assemble_matrix(V, kernel_py)
    M_f90 = assemble_matrix(V, kernel_f90)
    # ...

    assert_identical_coo(M_py, M_f90)

# ...


# ...
def test_3d_5():
    # ... define the weak formulation
    x,y,z = symbols('x y z')

    u = Symbol('u')
    v = Symbol('v')

    F = Field('F')

    a = Lambda((x,y,z,v,u), Dot(Grad(F*u), Grad(v)) + u*v)
    # ...

    # ...  create a finite element space
    p1  = 2 ; p2  = 2 ; p3  = 2
    ne1 = 2 ; ne2 = 2 ; ne3 = 2
    # ...

    print('> Grid   :: [{},{},{}]'.format(ne1, ne2, ne3))
    print('> Degree :: [{},{},{}]'.format(p1, p2, p3))

    grid_1 = linspace(0., 1., ne1+1)
    grid_2 = linspace(0., 1., ne2+1)
    grid_3 = linspace(0., 1., ne3+1)

    V1 = SplineSpace(p1, grid=grid_1)
    V2 = SplineSpace(p2, grid=grid_2)
    V3 = SplineSpace(p3, grid=grid_3)

    V = TensorSpace(V1, V2, V3)
    # ...

    eval_field_py  = compile_eval_field('kernel_5', a, V, backend='python')

#    # ...
#    kernel_py  = compile_kernel('kernel_1', a, V, backend='python')
#    kernel_f90 = compile_kernel('kernel_1', a, V, backend='fortran')
#
#    M_py  = assemble_matrix(V, kernel_py)
#    M_f90 = assemble_matrix(V, kernel_f90)
#    # ...
#
#    assert_identical_coo(M_py, M_f90)

# ...


# .....................................................
if __name__ == '__main__':

#    test_3d_1()
#    test_3d_2()
#    test_3d_3()
#    test_3d_4()
    test_3d_5()
