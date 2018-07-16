# coding: utf-8

# TODO add asserts

from gelato.core import BasicSobolevSpace
from gelato.core import SplineFemSpace
from gelato.core import TensorBasicSobolevSpace
from gelato.core import VectorBasicSobolevSpace
from gelato.core import TestFunction
from gelato.core import VectorTestFunction
from gelato.core import grad, dot, inner, cross, rot, curl, div
from gelato.core import dx, dy, dz

# ...
def test_fem_space_spline():
    print('============ test_fem_space_spline ==============')

    V = SplineFemSpace('V')

    print('> space :: ', V)
    print('> logical dim :: ', V.ldim)
    print('> degree :: ', V.degree)
    print('> n_elements :: ', V.n_elements)

    v = TestFunction(V, name='v')
    u = TestFunction(V, name='u')

    print('> test function  :: ', v)
    print('> trial function :: ', u)
# ...

# ...
def test_fem_space_tensor():
    print('============ test_fem_space_tensor ==============')

    V1 = SplineFemSpace('V1')
    V2 = SplineFemSpace('V2')
    V = TensorBasicSobolevSpace('V', V1, V2)

    print('> space :: ', V)
    print('> logical dim :: ', V.ldim)
    print('> degree :: ', V.degree)
    print('> n_elements :: ', V.n_elements)

    v = TestFunction(V, name='v')
    u = TestFunction(V, name='u')

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
    Vx = TensorBasicSobolevSpace('Vx', V1, V2)
    Vy = TensorBasicSobolevSpace('Vy', V2, V1)

    V = VectorBasicSobolevSpace('V', Vx, Vy)

    print('> space :: ', V)
    print('> logical dim :: ', V.ldim)
    print('> degree :: ', V.degree)
    print('> n_elements :: ', V.n_elements)

    v = TestFunction(V, name='v')
    u = TestFunction(V, name='u')

    print('> test function  :: ', v)
    print('> trial function :: ', u)
# ...

# .....................................................
if __name__ == '__main__':

    test_fem_space_spline()
    test_fem_space_tensor()
    test_fem_space_vector()
