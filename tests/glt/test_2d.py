# coding: utf-8

# TODO block example test 2 and 3

import numpy as np
from numpy import linspace, zeros, pi


from sympy.core.containers import Tuple
from sympy import symbols
from sympy import Symbol
from sympy import Lambda
from sympy import Function
from sympy import IndexedBase

from gelato.glt import glt_symbol
from gelato.calculus   import (Dot, Cross, Grad, Curl, Rot, Div)
from gelato.calculus   import Constant
from gelato.fem.utils    import compile_symbol

from spl.fem.splines import SplineSpace
from spl.fem.tensor  import TensorSpace
from spl.fem.vector  import VectorFemSpace


# ...
def test_2d_1():
    x,y = symbols('x y')

    u = Symbol('u')
    v = Symbol('v')

    a = Lambda((x,y,v,u), Dot(Grad(u), Grad(v)) + u*v)
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

    V = TensorSpace(V1, V2)
    # ...

    # ... create a glt symbol from a string without evaluation
    expr = glt_symbol(a, space=V)
    print('> glt symbol  := {0}'.format(expr))
    # ...

    # ...
    symbol_f90 = compile_symbol('symbol_1', a, V, backend='fortran')
    # ...

    # ... example of symbol evaluation
    t1 = linspace(-pi,pi, ne1+1)
    t2 = linspace(-pi,pi, ne2+1)
    x1 = linspace(0.,1., ne1+1)
    x2 = linspace(0.,1., ne2+1)
    e = zeros((ne1+1, ne2+1), order='F')
    symbol_f90(x1,x2,t1,t2, e)
    # ...

    print('')
# ...

# ...
# TODO
def test_2d_2():
    x,y = symbols('x y')

    u = IndexedBase('u')
    v = IndexedBase('v')

    a = Lambda((x,y,v,u), Rot(u) * Rot(v) + Div(u) * Div(v) + 0.2 * Dot(u, v))
    print('> input       := {0}'.format(a))

    # ... create a glt symbol from a string without evaluation
    #     a discretization is defined as a dictionary
    discretization = {"n_elements": [16, 16], "degrees": [3, 3]}

    expr = glt_symbol(a, discretization=discretization, evaluate=False,
                      is_block=True)
    print('> glt symbol  := {0}'.format(expr))
    # ...

    print('')
# ...

# ...
# TODO
def test_2d_3():
    x,y = symbols('x y')

    u = Symbol('u')
    v = Symbol('v')

    a = Lambda((x,y,v,u), Cross(Curl(u), Curl(v)) + 0.2 * u * v)
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

    V = TensorSpace(V1, V2)
    # ...

    # ... create a glt symbol from a string without evaluation
    expr = glt_symbol(a, space=V)
    print('> glt symbol  := {0}'.format(expr))
    # ...

    # ...
    symbol_f90 = compile_symbol('symbol_3', a, V, backend='fortran')
    # ...

    # ... example of symbol evaluation
    t1 = linspace(-pi,pi, ne1+1)
    t2 = linspace(-pi,pi, ne2+1)
    x1 = linspace(0.,1., ne1+1)
    x2 = linspace(0.,1., ne2+1)
    e = zeros((ne1+1, ne2+1), order='F')
    symbol_f90(x1,x2,t1,t2, e)
    # ...

    print('')
# ...

# ...
def test_2d_4():
    x,y = symbols('x y')

    u = Symbol('u')
    v = Symbol('v')

    c = Constant('c')

    b0 = Constant('b0')
    b1 = Constant('b1')
    b = Tuple(b0, b1)

    a = Lambda((x,y,v,u), c * u * v + Dot(b, Grad(v)) * u + Dot(b, Grad(u)) * v)
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

    V = TensorSpace(V1, V2)
    # ...

    # ... create a glt symbol from a string without evaluation
    expr = glt_symbol(a, space=V)
    print('> glt symbol  := {0}'.format(expr))
    # ...

    # ...
    symbol_f90 = compile_symbol('symbol_4', a, V,
                                d_constants={'b0': 0.1, 'b1': 1., 'c': 0.2},
                                backend='fortran')
    # ...

    # ... example of symbol evaluation
    t1 = linspace(-pi,pi, ne1+1)
    t2 = linspace(-pi,pi, ne2+1)
    x1 = linspace(0.,1., ne1+1)
    x2 = linspace(0.,1., ne2+1)
    e = zeros((ne1+1, ne2+1), order='F')
    symbol_f90(x1,x2,t1,t2, e)
    # ...

    print('')
# ...

# ...
def test_2d_5():
    x,y = symbols('x y')

    u = Symbol('u')
    v = Symbol('v')

    b = Function('b')

    a = Lambda((x,y,v,u), Dot(Grad(u), Grad(v)) + b(x,y)*u*v)
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

    V = TensorSpace(V1, V2)
    # ...

    # ... create a glt symbol from a string without evaluation
    expr = glt_symbol(a, space=V)
    print('> glt symbol  := {0}'.format(expr))
    # ...

    # ... user defined function
    def b(x,y):
        r = 1.+ x*y
        return r
    # ...

    # ... create an interactive pyccel context
    from pyccel.epyccel import ContextPyccel

    context = ContextPyccel(name='context_5')
    context.insert_function(b, ['double', 'double'], kind='function', results=['double'])

    context.compile()
    # ...

    # ...
    symbol_f90 = compile_symbol('symbol_5', a, V,
                                context=context,
                                backend='fortran')
    # ...

    # ... example of symbol evaluation
    t1 = linspace(-pi,pi, ne1+1)
    t2 = linspace(-pi,pi, ne2+1)
    x1 = linspace(0.,1., ne1+1)
    x2 = linspace(0.,1., ne2+1)
    e = zeros((ne1+1, ne2+1), order='F')
    symbol_f90(x1,x2,t1,t2, e)
    # ...

    print('')
# ...


# .....................................................
if __name__ == '__main__':
    test_2d_1()
#    test_2d_2()
#    test_2d_3()
    test_2d_4()
    test_2d_5()
