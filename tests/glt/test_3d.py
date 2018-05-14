# coding: utf-8

import numpy as np
from numpy import linspace, zeros, pi


from sympy.core.containers import Tuple
from sympy import symbols
from sympy import Symbol
from sympy import Lambda
from sympy import Function
from sympy import IndexedBase

from gelato.glt import glt_symbol
from gelato.calculus   import Constant
from gelato.calculus   import (Dot, Cross, Grad, Curl, Rot, Div, dx, dy, dz)
from gelato.fem.utils    import compile_symbol

from spl.fem.splines import SplineSpace
from spl.fem.tensor  import TensorSpace
from spl.fem.vector  import VectorFemSpace


# ...
def test_3d_scalar_1():
    print('============== test_3d_scalar_1 ================')

    x,y,z = symbols('x y z')

    u = Symbol('u')
    v = Symbol('v')

    a = Lambda((x,y,z,v,u), Dot(Grad(u), Grad(v)))
    print('> input       := {0}'.format(a))

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

    # ... create a glt symbol from a string without evaluation
    expr = glt_symbol(a, space=V)
    print('> glt symbol  := {0}'.format(expr))
    # ...

    # ...
    symbol_f90 = compile_symbol('symbol_scalar_1', a, V, backend='fortran')
    # ...

    # ... example of symbol evaluation
    t1 = linspace(-pi,pi, ne1+1)
    t2 = linspace(-pi,pi, ne2+1)
    t3 = linspace(-pi,pi, ne3+1)
    x1 = linspace(0.,1., ne1+1)
    x2 = linspace(0.,1., ne2+1)
    x3 = linspace(0.,1., ne3+1)
    e = zeros((ne1+1, ne2+1, ne3+1), order='F')
    symbol_f90(x1,x2,x3,t1,t2,t3, e)
    # ...

    print('')
# ...

# ...
def test_3d_scalar_2():
    print('============== test_3d_scalar_2 ================')

    x,y,z = symbols('x y z')

    u = Symbol('u')
    v = Symbol('v')

    c = Constant('c')

    b0 = Constant('b0')
    b1 = Constant('b1')
    b2 = Constant('b2')
    b = Tuple(b0, b1, b2)

    a = Lambda((x,y,z,v,u), c * u * v + Dot(b, Grad(v)) * u + Dot(b, Grad(u)) * v)
    print('> input       := {0}'.format(a))

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

    # ... create a glt symbol from a string without evaluation
    expr = glt_symbol(a, space=V)
    print('> glt symbol  := {0}'.format(expr))
    # ...

    # ...
    symbol_f90 = compile_symbol('symbol_scalar_2', a, V,
                                d_constants={'b0': 0.1, 'b1': 1., 'b2': 1., 'c': 0.2},
                                backend='fortran')
    # ...

    # ... example of symbol evaluation
    t1 = linspace(-pi,pi, ne1+1)
    t2 = linspace(-pi,pi, ne2+1)
    t3 = linspace(-pi,pi, ne3+1)
    x1 = linspace(0.,1., ne1+1)
    x2 = linspace(0.,1., ne2+1)
    x3 = linspace(0.,1., ne3+1)
    e = zeros((ne1+1, ne2+1, ne3+1), order='F')
    symbol_f90(x1,x2,x3,t1,t2,t3, e)
    # ...

    print('')
# ...

# ...
def test_3d_scalar_3():
    print('============== test_3d_scalar_3 ================')

    x,y,z = symbols('x y z')

    u = Symbol('u')
    v = Symbol('v')

    b = Function('b')

    a = Lambda((x,y,z,v,u), Dot(Grad(u), Grad(v)) + b(x,y,z)*u*v)
    print('> input       := {0}'.format(a))

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

    # ... create a glt symbol from a string without evaluation
    expr = glt_symbol(a, space=V)
    print('> glt symbol  := {0}'.format(expr))
    # ...

    # ... user defined function
    def b(x,y,z):
        r = 1.+ x*y*z
        return r
    # ...

    # ... create an interactive pyccel context
    from pyccel.epyccel import ContextPyccel

    context = ContextPyccel(name='context_scalar_3')
    context.insert_function(b, ['double', 'double', 'double'], kind='function', results=['double'])

    context.compile()
    # ...

    # ...
    symbol_f90 = compile_symbol('symbol_scalar_3', a, V,
                                context=context,
                                backend='fortran')
    # ...

    # ... example of symbol evaluation
    t1 = linspace(-pi,pi, ne1+1)
    t2 = linspace(-pi,pi, ne2+1)
    t3 = linspace(-pi,pi, ne3+1)
    x1 = linspace(0.,1., ne1+1)
    x2 = linspace(0.,1., ne2+1)
    x3 = linspace(0.,1., ne3+1)
    e = zeros((ne1+1, ne2+1, ne3+1), order='F')
    symbol_f90(x1,x2,x3,t1,t2,t3, e)
    # ...

    print('')
# ...

# ...
def test_3d_scalar_4():
    print('============== test_3d_scalar_4 ================')

    x,y,z = symbols('x y z')

    u = Symbol('u')
    v = Symbol('v')

    a = Lambda((x,y,z,v,u), dx(dx(u))*dx(dx(v)) + dy(dy(u))*dy(dy(v)) +
               dz(dz(u))*dz(dz(v)))
    print('> input       := {0}'.format(a))

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

    # ... create a glt symbol from a string without evaluation
    expr = glt_symbol(a, space=V)
    print('> glt symbol  := {0}'.format(expr))
    # ...

    # ...
    symbol_f90 = compile_symbol('symbol_scalar_4', a, V, backend='fortran')
    # ...

    # ... example of symbol evaluation
    t1 = linspace(-pi,pi, ne1+1)
    t2 = linspace(-pi,pi, ne2+1)
    t3 = linspace(-pi,pi, ne3+1)
    x1 = linspace(0.,1., ne1+1)
    x2 = linspace(0.,1., ne2+1)
    x3 = linspace(0.,1., ne3+1)
    e = zeros((ne1+1, ne2+1, ne3+1), order='F')
    symbol_f90(x1,x2,x3,t1,t2,t3, e)
    # ...

    print('')
# ...


# ...
def test_3d_block_1():
    print('============== test_3d_block_1 ================')

    x,y,z = symbols('x y z')

    u = IndexedBase('u')
    v = IndexedBase('v')


    a = Lambda((x,y,z,v,u), Div(u) * Div(v) + 0.2 * Dot(u, v))
    print('> input       := {0}'.format(a))

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

    Vx = TensorSpace(V1, V2, V3)
    Vy = TensorSpace(V1, V2, V3)
    Vz = TensorSpace(V1, V2, V3)

    V = VectorFemSpace(Vx, Vy, Vz)
    # ...

    # ... create a glt symbol from a string without evaluation
    expr = glt_symbol(a, space=V)
    print('> glt symbol  := {0}'.format(expr))
    # ...

    # ...
    symbol_f90 = compile_symbol('symbol_block_1', a, V, backend='fortran')
    # ...

    # ... example of symbol evaluation
    t1 = linspace(-pi,pi, ne1+1)
    t2 = linspace(-pi,pi, ne2+1)
    t3 = linspace(-pi,pi, ne3+1)
    x1 = linspace(0.,1., ne1+1)
    x2 = linspace(0.,1., ne2+1)
    x3 = linspace(0.,1., ne3+1)
    e = zeros((3, 3, ne1+1, ne2+1, ne3+1), order='F')
    symbol_f90(x1,x2,x3,t1,t2,t3, e)
    # ...

    print('')
# ...

# ...
def test_3d_block_2():
    print('============== test_3d_block_2 ================')

    x,y,z = symbols('x y z')

    u = IndexedBase('u')
    v = IndexedBase('v')

    a = Lambda((x,y,z,v,u), Dot(Curl(u), Curl(v)) + 0.2 * Dot(u, v))
    print('> input       := {0}'.format(a))

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

    Vx = TensorSpace(V1, V2, V3)
    Vy = TensorSpace(V1, V2, V3)
    Vz = TensorSpace(V1, V2, V3)

    V = VectorFemSpace(Vx, Vy, Vz)
    # ...

    # ... create a glt symbol from a string without evaluation
    expr = glt_symbol(a, space=V)
    print('> glt symbol  := {0}'.format(expr))
    # ...

    # ...
    symbol_f90 = compile_symbol('symbol_block_2', a, V, backend='fortran')
    # ...

    # ... example of symbol evaluation
    t1 = linspace(-pi,pi, ne1+1)
    t2 = linspace(-pi,pi, ne2+1)
    t3 = linspace(-pi,pi, ne3+1)
    x1 = linspace(0.,1., ne1+1)
    x2 = linspace(0.,1., ne2+1)
    x3 = linspace(0.,1., ne3+1)
    e = zeros((3, 3, ne1+1, ne2+1, ne3+1), order='F')
    symbol_f90(x1,x2,x3,t1,t2,t3, e)
    # ...

    print('')
# ...



# ...
def test_3d_block_3():
    print('============== test_3d_block_3 ================')

    x,y,z = symbols('x y z')

    u = IndexedBase('u')
    v = IndexedBase('v')

    b = Tuple(1.0, 0., 0.)

    a = Lambda((x,y,z,v,u), Dot(Curl(Cross(b,u)), Curl(Cross(b,v))) + 0.2 * Dot(u, v))
    print('> input       := {0}'.format(a))

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

    Vx = TensorSpace(V1, V2, V3)
    Vy = TensorSpace(V1, V2, V3)
    Vz = TensorSpace(V1, V2, V3)

    V = VectorFemSpace(Vx, Vy, Vz)
    # ...

    # ... create a glt symbol from a string without evaluation
    expr = glt_symbol(a, space=V)
    print('> glt symbol  := {0}'.format(expr))
    # ...

    # ...
    symbol_f90 = compile_symbol('symbol_block_3', a, V, backend='fortran')
    # ...

    # ... example of symbol evaluation
    t1 = linspace(-pi,pi, ne1+1)
    t2 = linspace(-pi,pi, ne2+1)
    t3 = linspace(-pi,pi, ne3+1)
    x1 = linspace(0.,1., ne1+1)
    x2 = linspace(0.,1., ne2+1)
    x3 = linspace(0.,1., ne3+1)
    e = zeros((3, 3, ne1+1, ne2+1, ne3+1), order='F')
    symbol_f90(x1,x2,x3,t1,t2,t3, e)
    # ...

    print('')
# ...

# ...
def test_3d_block_4():
    print('============== test_3d_block_4 ================')

    """Alfven operator."""
    x,y,z = symbols('x y z')

    u = IndexedBase('u')
    v = IndexedBase('v')

    bx = Constant('bx')
    by = Constant('by')
    bz = Constant('bz')
    b = Tuple(bx, by, bz)

    c0 = Constant('c0')
    c1 = Constant('c1')
    c2 = Constant('c2')

    a = Lambda((x,y,z,v,u), (  c0 * Dot(u, v)
                             - c1 * Div(u) * Div(v)
                             + c2 * Dot(Curl(Cross(b,u)), Curl(Cross(b,v)))))
    print('> input       := {0}'.format(a))

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

    Vx = TensorSpace(V1, V2, V3)
    Vy = TensorSpace(V1, V2, V3)
    Vz = TensorSpace(V1, V2, V3)

    V = VectorFemSpace(Vx, Vy, Vz)
    # ...

    # ... create a glt symbol from a string without evaluation
    expr = glt_symbol(a, space=V)
    print('> glt symbol  := {0}'.format(expr))
    # ...

    # ...
    symbol_f90 = compile_symbol('symbol_block_4', a, V,
                                d_constants={'bx': 0.1, 'by': 1., 'bz': 0.2,
                                             'c0': 0.1, 'c1': 1., 'c2': 1.},
                                backend='fortran')
    # ...

    # ... example of symbol evaluation
    t1 = linspace(-pi,pi, ne1+1)
    t2 = linspace(-pi,pi, ne2+1)
    t3 = linspace(-pi,pi, ne3+1)
    x1 = linspace(0.,1., ne1+1)
    x2 = linspace(0.,1., ne2+1)
    x3 = linspace(0.,1., ne3+1)
    e = zeros((3, 3, ne1+1, ne2+1, ne3+1), order='F')
    symbol_f90(x1,x2,x3,t1,t2,t3, e)
    # ...

    print('')
# ...

# .....................................................
if __name__ == '__main__':
    test_3d_scalar_1()
    test_3d_scalar_2()
    test_3d_scalar_3()
    test_3d_scalar_4()
    test_3d_block_1()
    test_3d_block_2()
    test_3d_block_3()
    test_3d_block_4()
