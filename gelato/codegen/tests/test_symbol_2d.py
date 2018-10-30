# coding: utf-8

from numpy import linspace, zeros, pi

from sympy import Symbol
from sympy.core.containers import Tuple
from sympy import symbols
from sympy import srepr

from sympde.core import dx, dy, dz
from sympde.core import Domain
from sympde.core import Constant
from sympde.core import Field
from sympde.core import grad, dot, inner, cross, rot, curl, div
from sympde.core import FunctionSpace
from sympde.core import TestFunction
from sympde.core import VectorTestFunction
from sympde.core import BilinearForm

from gelato.codegen import compile_symbol

DIM = 2
domain = Domain('Omega', dim=DIM)

# ...
def test_symbol_2d_1():
    print('============ test_symbol_2d_1 =============')

    # ... abstract model
    V = FunctionSpace('V', domain)

    v = TestFunction(V, name='v')
    u = TestFunction(V, name='u')

    mass = BilinearForm((v,u), u*v)
    laplace = BilinearForm((v,u), dot(grad(v), grad(u)))
    # ...

    # ...
    degrees = [1,1]
    n_elements = [4,4]

    symbol = compile_symbol('mass_symbol', mass, degrees,
                            n_elements=n_elements)
    # ...

    # ...
    n1 = 101 ; n2 = 101

    t1 = linspace(-pi, pi, n1)
    t2 = linspace(-pi, pi, n2)
    x1 = linspace(0.,1., n1)
    x2 = linspace(0.,1., n2)

    xs = [x1, x2]
    ts = [t1, t2]

    e = zeros((n1, n2))
    symbol(*xs, *ts, e)
    print('> ', e.min(), e.max())
    # ...
# ...

# ...
def test_symbol_2d_2():
    print('============ test_symbol_2d_2 =============')

    # ... abstract model
    V = FunctionSpace('V', domain)

    v = TestFunction(V, name='v')
    u = TestFunction(V, name='u')

    mass = BilinearForm((v,u), u*v)
    laplace = BilinearForm((v,u), dot(grad(v), grad(u)))
    # ...

    # ...
    degrees = [1,1]

    symbol = compile_symbol('mass_symbol', mass, degrees)
    # ...

    # ...
    n1 = 5 ; n2 = 5

    t1 = linspace(-pi, pi, n1)
    t2 = linspace(-pi, pi, n2)
    x1 = linspace(0.,1., n1)
    x2 = linspace(0.,1., n2)

    xs = [x1, x2]
    ts = [t1, t2]

    e = zeros((n1, n2))
    # ...

    # ...
    n_elements = [4, 4]
    symbol(*xs, *ts, e, *n_elements)
    print('> ', e.min(), e.max())
    # ...

    # ...
    n_elements = [4, 8]
    symbol(*xs, *ts, e, *n_elements)
    print('> ', e.min(), e.max())
    # ...

# ...

# .....................................................
if __name__ == '__main__':
    test_symbol_2d_1()
    test_symbol_2d_2()
