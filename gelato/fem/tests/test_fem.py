# coding: utf-8

# TODO add asserts

from gelato.fem.core import SplineFemSpace
from gelato.fem.core import TensorFemSpace
from gelato.fem.core import VectorFemSpace

from gelato.fem.core import TestFunction
from gelato.fem.core import TrialFunction

from gelato.calculus import grad, dot, inner

# ...
def test_fem_space_spline():
    print('============ test_fem_space_spline ==============')

    V = SplineFemSpace('V')

    print('> space :: ', V)
    print('> degree :: ', V.degree)
    print('> n_elements :: ', V.n_elements)

    v = TestFunction(V, name='v')
    u = TrialFunction(V, name='u')

    print('> test function  :: ', v)
    print('> trial function :: ', u)
# ...

# ...
def test_fem_space_tensor():
    print('============ test_fem_space_tensor ==============')

    V1 = SplineFemSpace('V1')
    V2 = SplineFemSpace('V2')
    V = TensorFemSpace('V', V1, V2)

    print('> space :: ', V)
    print('> degree :: ', V.degree)
    print('> n_elements :: ', V.n_elements)

    v = TestFunction(V, name='v')
    u = TrialFunction(V, name='u')

    print('> test function  :: ', v)
    print('> trial function :: ', u)

    # ...
    expr = inner(grad(v), grad(u))
    print('> expr := {0}'.format(expr))
    # ...

# ...

# ...
def test_fem_space_vector():
    print('============ test_fem_space_vector ==============')

    V1 = SplineFemSpace('V1')
    V2 = SplineFemSpace('V2')
    Vx = TensorFemSpace('Vx', V1, V2)
    Vy = TensorFemSpace('Vy', V2, V1)

    V = VectorFemSpace('V', Vx, Vy)

    print('> space :: ', V)
    print('> degree :: ', V.degree)
    print('> n_elements :: ', V.n_elements)

    v = TestFunction(V, name='v')
    u = TrialFunction(V, name='u')

    print('> test function  :: ', v)
    print('> trial function :: ', u)
# ...

# .....................................................
if __name__ == '__main__':

    test_fem_space_spline()
    test_fem_space_tensor()
    test_fem_space_vector()
