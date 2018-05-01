# coding: utf-8

# TODO: - allow for giving a name for the trial/test basis
#       - Ni/Nj should be Ni_0/Nj_O

from numpy import linspace

from sympy.core.containers import Tuple
from sympy import symbols
from sympy import Symbol
from sympy import Lambda

from gelato.expression   import construct_weak_form
from gelato.calculus     import (Dot, Cross, Grad, Curl, Rot, Div)
from gelato.calculus     import Constant
from gelato.fem.assembly import assemble_matrix
from gelato.fem.utils    import compile_kernel

from spl.fem.splines import SplineSpace
from spl.fem.tensor  import TensorSpace

# ...
def test_1d_1():
    # ... define the weak formulation
    x = Symbol('x')

    u = Symbol('u')
    v = Symbol('v')

    expr = Lambda((x,v,u), Dot(Grad(u), Grad(v)) + u*v)
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
    kernel = compile_kernel('kernel', expr, V)
    M = assemble_matrix(V, kernel).tocsr()
    # ...

    print(M.shape)
# ...

# ...
def test_1d_2():
    # ... define the weak formulation
    x = Symbol('x')

    u = Symbol('u')
    v = Symbol('v')

    alpha = Symbol('alpha')
    nu = Symbol('nu')

    expr = Lambda((x,v,u), alpha * Dot(Grad(u), Grad(v)) + nu*u*v)
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
    kernel = compile_kernel('kernel', expr, V,
                            d_constants={'nu': 0.1},
                            d_args={'alpha': 'double'})
    M = assemble_matrix(V, kernel, args={'alpha': 2.0}).tocsr()
    # ...

    print(M.shape)
# ...


# .....................................................
if __name__ == '__main__':

    test_1d_1()
    test_1d_2()
