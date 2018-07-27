# coding: utf-8

from numpy import linspace, zeros, pi

from sympy import Symbol
from sympy.core.containers import Tuple
from sympy import symbols
from sympy import srepr

from symfe.core import dx, dy, dz
from symfe.core import Constant
from symfe.core import Field
from symfe.core import grad, dot, inner, cross, rot, curl, div
from symfe.core import H1Space
from symfe.core import TestFunction
from symfe.core import VectorTestFunction
from symfe.core import BilinearForm

from gelato.codegen import compile_symbol

# ...
def test_symbol_3d_1():
    print('============ test_symbol_3d_1 =============')

    # ... abstract model
    V = H1Space('V', ldim=3)

    v = TestFunction(V, name='v')
    u = TestFunction(V, name='u')

    mass = BilinearForm((v,u), u*v)
    laplace = BilinearForm((v,u), dot(grad(v), grad(u)))
    # ...

    # ...
    degrees = [1,1,1]
    n_elements = [4,4,4]

    symbol = compile_symbol('mass_symbol', mass, degrees,
                            n_elements=n_elements)
    # ...

    # ...
    n1 = 21 ; n2 = 21 ; n3 = 21

    t1 = linspace(-pi, pi, n1)
    t2 = linspace(-pi, pi, n2)
    t3 = linspace(-pi, pi, n3)
    x1 = linspace(0.,1., n1)
    x2 = linspace(0.,1., n2)
    x3 = linspace(0.,1., n3)

    xs = [x1, x2, x3]
    ts = [t1, t2, t3]

    e = zeros((n1, n2, n3))
    symbol(*xs, *ts, e)
    print('> ', e.min(), e.max())
    # ...
# ...

# ...
def test_symbol_3d_2():
    print('============ test_symbol_3d_2 =============')

    # ... abstract model
    V = H1Space('V', ldim=3)

    v = TestFunction(V, name='v')
    u = TestFunction(V, name='u')

    mass = BilinearForm((v,u), u*v)
    laplace = BilinearForm((v,u), dot(grad(v), grad(u)))
    # ...

    # ...
    degrees = [1,1,1]

    symbol = compile_symbol('mass_symbol', mass, degrees)
    # ...

    # ...
    n1 = 21 ; n2 = 21 ; n3 = 21

    t1 = linspace(-pi, pi, n1)
    t2 = linspace(-pi, pi, n2)
    t3 = linspace(-pi, pi, n3)
    x1 = linspace(0.,1., n1)
    x2 = linspace(0.,1., n2)
    x3 = linspace(0.,1., n3)

    xs = [x1, x2, x3]
    ts = [t1, t2, t3]

    e = zeros((n1, n2, n3))
    # ...

    # ...
    n_elements = [4, 4, 4]
    symbol(*xs, *ts, e, *n_elements)
    print('> ', e.min(), e.max())
    # ...

    # ...
    n_elements = [4, 8, 8]
    symbol(*xs, *ts, e, *n_elements)
    print('> ', e.min(), e.max())
    # ...

# ...

# .....................................................
if __name__ == '__main__':
    test_symbol_3d_1()
    test_symbol_3d_2()
