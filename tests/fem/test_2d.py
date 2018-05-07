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
from gelato.fem.assembly import assemble_matrix
from gelato.fem.utils    import compile_kernel

from spl.fem.splines import SplineSpace
from spl.fem.tensor  import TensorSpace
from spl.fem.vector  import VectorFemSpace

from utils import assert_identical_coo

# ...
def test_2d_1():
    # ... define the weak formulation
    x,y = symbols('x y')

    u = Symbol('u')
    v = Symbol('v')

    a = Lambda((x,y,v,u), Dot(Grad(u), Grad(v)) + u*v)
    # ...

    # ...  create a finite element space
    p1  = 2 ; p2  = 2
    ne1 = 8 ; ne2 = 8

    print('> Grid   :: [{ne1},{ne2}]'.format(ne1=ne1, ne2=ne2))
    print('> Degree :: [{p1},{p2}]'.format(p1=p1, p2=p2))

    grid_1 = linspace(0., 1., ne1+1)
    grid_2 = linspace(0., 1., ne2+1)

    V1 = SplineSpace(p1, grid=grid_1)
    V2 = SplineSpace(p2, grid=grid_2)

    V = TensorSpace(V1, V2)
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
def test_2d_2():
    # ... define the weak formulation
    x,y = symbols('x y')

    u = Symbol('u')
    v = Symbol('v')

    alpha = Symbol('alpha')
    nu = Symbol('nu')

    a = Lambda((x,v,u), alpha * Dot(Grad(u), Grad(v)) + nu*u*v)
    # ...

    # ...  create a finite element space
    p1  = 2 ; p2  = 2
    ne1 = 8 ; ne2 = 8

    print('> Grid   :: [{ne1},{ne2}]'.format(ne1=ne1, ne2=ne2))
    print('> Degree :: [{p1},{p2}]'.format(p1=p1, p2=p2))

    grid_1 = linspace(0., 1., ne1+1)
    grid_2 = linspace(0., 1., ne2+1)

    V1 = SplineSpace(p1, grid=grid_1)
    V2 = SplineSpace(p2, grid=grid_2)

    V = TensorSpace(V1, V2)
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
def test_2d_3():
    # ... define the weak formulation
    x,y = symbols('x y')

    u = IndexedBase('u')
    v = IndexedBase('v')

    a = Lambda((x,y,v,u), Rot(u) * Rot(v) + Div(u) * Div(v) + 0.2 * Dot(u, v))
    # ...

    # ...  create a finite element space
    p1  = 2 ; p2  = 2
    ne1 = 8 ; ne2 = 8

    print('> Grid   :: [{ne1},{ne2}]'.format(ne1=ne1, ne2=ne2))
    print('> Degree :: [{p1},{p2}]'.format(p1=p1, p2=p2))

    grid_1 = linspace(0., 1., ne1+1)
    grid_2 = linspace(0., 1., ne2+1)

    V1 = SplineSpace(p1, grid=grid_1)
    V2 = SplineSpace(p2, grid=grid_2)

    W = TensorSpace(V1, V2)
    # ...

    # ... vector space
    V = VectorFemSpace(W, W)
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
def test_2d_4():
    # ... define the weak formulation
    x,y = symbols('x y')

    u = Symbol('u')
    v = Symbol('v')

    b = Function('b')

    a = Lambda((x,y,v,u), Dot(Grad(u), Grad(v)) + b(x,y)*u*v)
    # ...

    # ...  create a finite element space
    p1  = 2 ; p2  = 2
    ne1 = 8 ; ne2 = 8

    print('> Grid   :: [{ne1},{ne2}]'.format(ne1=ne1, ne2=ne2))
    print('> Degree :: [{p1},{p2}]'.format(p1=p1, p2=p2))

    grid_1 = linspace(0., 1., ne1+1)
    grid_2 = linspace(0., 1., ne2+1)

    V1 = SplineSpace(p1, grid=grid_1)
    V2 = SplineSpace(p2, grid=grid_2)

    V = TensorSpace(V1, V2)
    # ...

    # ...
    header_b = '#$ header function b(double, double) results(double)'
    def b(x,y):
        r = 1.+ x*(1.-x) + y*(1.-y)
        return r

    kernel_py  = compile_kernel('kernel_4', a, V,
                                d_functions={'b': (b, header_b)},
                                verbose=True,
                                backend='python')

    kernel_f90 = compile_kernel('kernel_4', a, V,
                                d_functions={'b': (b, header_b)},
                                verbose=True,
                                backend='fortran')

    M_py  = assemble_matrix(V, kernel_py)
    M_f90 = assemble_matrix(V, kernel_f90)
    # ...

    assert_identical_coo(M_py, M_f90)

# ...


# .....................................................
if __name__ == '__main__':

    test_2d_1()
    test_2d_2()
    test_2d_3()
    test_2d_4()
