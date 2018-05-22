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
from gelato.calculus     import (Dot, Cross, Grad, Curl, Rot, Div, dx, dy)
from gelato.calculus     import Constant
from gelato.calculus     import Field
from gelato.fem.assembly import assemble_matrix
from gelato.fem.utils    import compile_kernel

from spl.fem.splines import SplineSpace
from spl.fem.tensor  import TensorFemSpace
from spl.fem.vector  import VectorFemSpace
from spl.fem.splines import Spline

from utils import assert_identical_coo

# ...
def test_2d_scalar_1():
    print('============== test_2d_scalar_1 ================')

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

    V = TensorFemSpace(V1, V2)
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
def test_2d_scalar_2():
    print('============== test_2d_scalar_2 ================')

    # ... define the weak formulation
    x,y = symbols('x y')

    u = Symbol('u')
    v = Symbol('v')

    alpha = Constant('alpha')
    nu = Constant('nu')

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

    V = TensorFemSpace(V1, V2)
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
def test_2d_scalar_3():
    print('============== test_2d_scalar_3 ================')

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

    V = TensorFemSpace(V1, V2)
    # ...

    # ... user defined function
    def b(x,y):
        r = 1.+ x*(1.-x) + y*(1.-y)
        return r
    # ...

    # ... create an interactive pyccel context
    from pyccel.epyccel import ContextPyccel

    context = ContextPyccel(name='context_3')
    context.insert_function(b, ['double', 'double'], kind='function', results=['double'])

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

# ...

# ...
def test_2d_scalar_4():
    print('============== test_2d_scalar_4 ================')

    # ... define the weak formulation
    x,y = symbols('x y')

    u = Symbol('u')
    v = Symbol('v')

    b0 = Function('b0')
    b1 = Function('b1')

    a = Lambda((x,y,v,u),
               (b0(x,y)*dx(v) + b1(x,y)*dy(v)) * (b0(x,y)*dx(u) + b1(x,y)*dy(u)))
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

    V = TensorFemSpace(V1, V2)
    # ...

    # ... user defined function
    def b0(x,y):
        from numpy import sin
        from scipy import pi

        r = 1.1659397624413860850012270020670 * (1.0 + 0.1 * sin(2*pi*y))
        return r

    def b1(x,y):
        from numpy import sin
        from scipy import pi

        r = 1.0 * (1.0 + 0.1 * sin(2*pi*y))
        return r
    # ...

    # ... create an interactive pyccel context
    from pyccel.epyccel import ContextPyccel

    context = ContextPyccel(name='context_4')
    context.insert_function(b0, ['double', 'double'], kind='function', results=['double'])
    context.insert_function(b1, ['double', 'double'], kind='function', results=['double'])

    context.compile()
    # ...

    # ...
    kernel_py  = compile_kernel('kernel_scalar_4', a, V,
                                context=context,
                                verbose=True,
                                backend='python')

    kernel_f90 = compile_kernel('kernel_scalar_4', a, V,
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
def test_2d_scalar_5():
    print('============== test_2d_scalar_5 ================')

    # ... define the weak formulation
    x,y = symbols('x y')

    u = Symbol('u')
    v = Symbol('v')

    F = Field('F')

    a = Lambda((x,y,v,u), Dot(Grad(F*u), Grad(v)) + u*v)
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

    V = TensorFemSpace(V1, V2)
    # ...

    F = Spline(V)
    F.coeffs._data[:,:] = 1.

    # ...
    kernel_py  = compile_kernel('kernel_scalar_5', a, V, backend='python')
    kernel_f90 = compile_kernel('kernel_scalar_5', a, V, backend='fortran')

    M_py  = assemble_matrix(V, kernel_py, fields={'F': F})
    M_f90 = assemble_matrix(V, kernel_f90, fields={'F': F})
    # ...

    assert_identical_coo(M_py, M_f90)

# ...

# ...
def test_2d_scalar_6():
    print('============== test_2d_scalar_6 ================')

    # ... define the weak formulation
    x,y = symbols('x y')

    u = Symbol('u')
    v = Symbol('v')

    a = Lambda((x,y,v,u), dx(dx(u))*dx(dx(v)) + dy(dy(u))*dy(dy(v)) + Dot(Grad(u), Grad(v)) + u*v)
    # ...

    # ...  create a finite element space
    p1  = 2 ; p2  = 2
    ne1 = 8 ; ne2 = 8

    print('> Grid   :: [{ne1},{ne2}]'.format(ne1=ne1, ne2=ne2))
    print('> Degree :: [{p1},{p2}]'.format(p1=p1, p2=p2))

    grid_1 = linspace(0., 1., ne1+1)
    grid_2 = linspace(0., 1., ne2+1)

    V1 = SplineSpace(p1, grid=grid_1, nderiv=2)
    V2 = SplineSpace(p2, grid=grid_2, nderiv=2)

    V = TensorFemSpace(V1, V2)
    # ...

    # ...
    kernel_py  = compile_kernel('kernel_scalar_6', a, V, backend='python')
    kernel_f90 = compile_kernel('kernel_scalar_6', a, V, backend='fortran')

    M_py  = assemble_matrix(V, kernel_py)
    M_f90 = assemble_matrix(V, kernel_f90)
    # ...

    assert_identical_coo(M_py, M_f90)

# ...

# ...
def test_2d_block_1():
    print('============== test_2d_block_1 ================')

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

    W = TensorFemSpace(V1, V2)
    # ...

    # ... vector space
    V = VectorFemSpace(W, W)
    # ...

    # ...
    kernel_py  = compile_kernel('kernel_block_1', a, V, backend='python')
    kernel_f90 = compile_kernel('kernel_block_1', a, V, backend='fortran')

    M_py  = assemble_matrix(V, kernel_py)
    M_f90 = assemble_matrix(V, kernel_f90)
    # ...

    assert_identical_coo(M_py, M_f90)

# ...

# ...
def test_2d_block_2():
    print('============== test_2d_block_2 ================')

    # ... define the weak formulation
    x,y = symbols('x y')

    u = IndexedBase('u')
    v = IndexedBase('v')

    F = Field('F')

    a = Lambda((x,y,v,u), Rot(u) * Rot(v) + Div(u) * Div(v) + 0.2 * Dot(u, v) + F*u[0]*v[0])
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

    W = TensorFemSpace(V1, V2)
    # ...

    # ... vector space
    V = VectorFemSpace(W, W)
    # ...

    F = Spline(W)
    F.coeffs._data[:,:] = 1.

    # ...
    kernel_py  = compile_kernel('kernel_block_2', a, V, backend='python')
    kernel_f90 = compile_kernel('kernel_block_2', a, V, backend='fortran')

    M_py  = assemble_matrix(V, kernel_py, fields={'F': F})
    M_f90 = assemble_matrix(V, kernel_f90, fields={'F': F})
    # ...

    assert_identical_coo(M_py, M_f90)

# ...

# ...
def test_2d_block_3():
    print('============== test_2d_block_3 ================')

    x, y = symbols('x y')

    u = Symbol('u')
    v = Symbol('v')

    epsilon = Constant('epsilon')

    Laplace = lambda v,u: Dot(Grad(v), Grad(u))
    Mass = lambda v,u: v*u

    u1,u2,p = symbols('u1 u2 p')
    v1,v2,q = symbols('v1 v2 q')

    a = Lambda((x,y,v1,v2,q,u1,u2,p),
                 Laplace(v1,u1) - dx(v1) * p
               + Laplace(v2,u2) - dy(v2) * p
               + q * (dx(u1) + dy(u2))
               + epsilon * Mass(q,p))

    print('> input       := {0}'.format(a))

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
    V = VectorFemSpace(V, V, V)
    # ...

    # ...
    kernel_py  = compile_kernel('kernel_block_3', a, V,
                                d_args={'epsilon': 'double'},
                                backend='python')
    kernel_f90 = compile_kernel('kernel_block_3', a, V,
                                d_args={'epsilon': 'double'},
                                backend='fortran')

    M_py  = assemble_matrix(V, kernel_py, args={'epsilon': 1.e-3})
    M_f90 = assemble_matrix(V, kernel_f90, args={'epsilon': 1.e-3})
    # ...

    assert_identical_coo(M_py, M_f90)

    print('')
# ...


# .....................................................
if __name__ == '__main__':

#    test_2d_scalar_1()
#    test_2d_scalar_2()
#    test_2d_scalar_3()
#    test_2d_scalar_4()
#    test_2d_scalar_5()
#    test_2d_scalar_6()
#    test_2d_block_1()
#    test_2d_block_2()
    test_2d_block_3()
