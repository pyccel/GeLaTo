# coding: utf-8

# TODO - python version of symbols
#      - block example (test_3)

import numpy as np
from numpy import linspace, zeros, pi

from sympy.core.containers import Tuple
from sympy import symbols
from sympy import Symbol
from sympy import Lambda
from sympy import Function

from gelato.glt import glt_symbol
from gelato.calculus   import (Dot, Cross, Grad, Curl, Rot, Div, dx)
from gelato.calculus   import Constant
from gelato.fem.utils  import compile_symbol

from spl.fem.splines import SplineSpace
from spl.fem.tensor  import TensorSpace
from spl.fem.vector  import VectorFemSpace



# ...
def test_1d_1():
    x = Symbol('x')

    u = Symbol('u')
    v = Symbol('v')

    a = Lambda((x,v,u), Dot(Grad(u), Grad(v)) + u*v)
    print('> input       := {0}'.format(a))

    # ...  create a finite element space
    p  = 3
    ne = 64

    print('> Grid   :: {ne}'.format(ne=ne))
    print('> Degree :: {p}'.format(p=p))

    grid = linspace(0., 1., ne+1)

    V = SplineSpace(p, grid=grid)
    # ...

    # ... create a glt symbol from a string without evaluation
    expr = glt_symbol(a, space=V)
    print('> glt symbol  := {0}'.format(expr))
    # ...

    # ...
    symbol_f90 = compile_symbol('symbol_1', a, V, backend='fortran')
    # ...

    # ... example of symbol evaluation
    t1 = linspace(-pi,pi, ne+1)
    x1 = linspace(0.,1., ne+1)
    e = zeros(ne+1)
    symbol_f90(x1, t1, e)
    # ...

    print('')
# ...

# ...
def test_1d_2():
    x = Symbol('x')

    u = Symbol('u')
    v = Symbol('v')

    b = Constant('b')

    a = Lambda((x,v,u), Dot(Grad(b*u), Grad(v)) + u*v)
    print('> input       := {0}'.format(a))

    # ...  create a finite element space
    p  = 3
    ne = 64

    print('> Grid   :: {ne}'.format(ne=ne))
    print('> Degree :: {p}'.format(p=p))

    grid = linspace(0., 1., ne+1)

    V = SplineSpace(p, grid=grid)
    # ...

    # ... create a glt symbol from a string without evaluation
    expr = glt_symbol(a, space=V)
    print('> glt symbol  := {0}'.format(expr))
    # ...

    # ...
    symbol_f90 = compile_symbol('symbol_2', a, V,
                                d_constants={'b': 0.1},
                                backend='fortran')
    # ...

    # ... example of symbol evaluation
    t1 = linspace(-pi,pi, ne+1)
    x1 = linspace(0.,1., ne+1)
    e = zeros(ne+1)
    symbol_f90(x1, t1, e)
    # ...

    print('')
# ...

# ...
def test_1d_3():
    x = Symbol('x')

    u0, u1 = symbols('u0 u1')
    v0, v1 = symbols('v0 v1')

    a = Lambda((x,v0,v1,u0,u1), dx(u0)*dx(v0) + dx(u1)*v0 + u0*dx(v1) + u1*v1)
    print('> input       := {0}'.format(a))

    # ... create a glt symbol from a string without evaluation
    #     a discretization is defined as a dictionary
    discretization = {"n_elements": [16], "degrees": [3]}

    expr = glt_symbol(a, discretization=discretization, evaluate=False,
                      is_block=True)
    print('> glt symbol  := {0}'.format(expr))
    # ...

    print('')
# ...

# ...
def test_1d_4():
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

    # ... create a glt symbol from a string without evaluation
    expr = glt_symbol(a, space=V)
    print('> glt symbol  := {0}'.format(expr))
    # ...

    # ... user defined function
    def b(s):
        r = 1.+ s*(1.-s)
        return r
    # ...

    # ... create an interactive pyccel context
    from pyccel.epyccel import ContextPyccel

    context = ContextPyccel(name='context_4')
    context.insert_function(b, ['double'], kind='function', results=['double'])

    context.compile()
    # ...

    # ...
    symbol_f90 = compile_symbol('symbol_4', a, V,
                                context=context,
                                backend='fortran')
    # ...

    # ... example of symbol evaluation
    t1 = linspace(-pi,pi, ne+1)
    x1 = linspace(0.,1., ne+1)
    e = zeros(ne+1)
    symbol_f90(x1, t1, e)
    # ...

    print('')
# ...

# .....................................................
if __name__ == '__main__':

    test_1d_1()
    test_1d_2()
#    test_1d_3()
    test_1d_4()
