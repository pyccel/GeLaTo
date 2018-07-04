# coding: utf-8

# TODO add asserts

from sympy import Symbol
from sympy.core.containers import Tuple

from gelato.fem.core import FemSpace
from gelato.fem.core import SplineFemSpace
from gelato.fem.core import TensorFemSpace
from gelato.fem.core import VectorFemSpace
from gelato.fem.core import TestFunction
from gelato.fem.core import TrialFunction
from gelato.fem.core import VectorTestFunction
from gelato.fem.core import VectorTrialFunction
from gelato.fem.expr import gelatize
from gelato.fem.expr import normalize_weak_from

from gelato.calculus import grad, dot, inner, cross, rot, curl, div
from gelato.calculus import dx, dy, dz

# ...
def test_fem_space_spline():
    print('============ test_fem_space_spline ==============')

    V = SplineFemSpace('V')

    print('> space :: ', V)
    print('> logical dim :: ', V.ldim)
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
    print('> logical dim :: ', V.ldim)
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
    print('> logical dim :: ', V.ldim)
    print('> degree :: ', V.degree)
    print('> n_elements :: ', V.n_elements)

    v = TestFunction(V, name='v')
    u = TrialFunction(V, name='u')

    print('> test function  :: ', v)
    print('> trial function :: ', u)
# ...

# ...
def test_trialtest_2d_1():
    print('============ test_trialtest_2d_1 =============')

    V = FemSpace('V', ldim=2, is_vector=True, shape=2)

    v = VectorTestFunction(V, name='v')

    assert(gelatize(rot(v)) == -dx(v[1]) + dy(v[0]))
    assert(gelatize(div(v)) == dx(v[0]) + dy(v[1]))

    v0_x = Symbol('v_x[0]')
    v1_x = Symbol('v_x[1]')
    v0_y = Symbol('v_y[0]')
    v1_y = Symbol('v_y[1]')

    assert(normalize_weak_from(rot(v)) == -v1_x + v0_y)
    assert(normalize_weak_from(div(v)) == v0_x + v1_y)

#    expr = div(v)
#    print('> input         >>> {0}'.format(expr))
#    print('> gelatized     >>> {0}'.format(gelatize(expr)))
#    print('> normal form   >>> {0}'.format(normalize_weak_from(expr)))
# ...

# ...
def test_trialtest_3d_1():
    print('============ test_trialtest_3d_1 =============')

    V = FemSpace('V', ldim=3, is_vector=True, shape=3)

    v = VectorTestFunction(V, name='v')

    assert(gelatize(curl(v)) == Tuple( dy(v[2]) - dz(v[1]),
                                      -dx(v[2]) + dz(v[0]),
                                       dx(v[1]) - dy(v[0])))
    assert(gelatize(div(v)) == dx(v[0]) + dy(v[1]) + dz(v[2]))

    v0_x = Symbol('v_x[0]')
    v1_x = Symbol('v_x[1]')
    v2_x = Symbol('v_x[2]')
    v0_y = Symbol('v_y[0]')
    v1_y = Symbol('v_y[1]')
    v2_y = Symbol('v_y[2]')
    v0_z = Symbol('v_z[0]')
    v1_z = Symbol('v_z[1]')
    v2_z = Symbol('v_z[2]')

    assert(normalize_weak_from(curl(v)) == Tuple( v2_y - v1_z,
                                                 -v2_x + v0_z,
                                                  v1_x - v0_y))
    assert(normalize_weak_from(div(v)) == v0_x + v1_y + v2_z)

#    expr = curl(v)
#    print('> input         >>> {0}'.format(expr))
#    print('> gelatized     >>> {0}'.format(gelatize(expr)))
#    print('> normal form   >>> {0}'.format(normalize_weak_from(expr)))
# ...

# .....................................................
if __name__ == '__main__':

#    test_fem_space_spline()
#    test_fem_space_tensor()
#    test_fem_space_vector()
#    test_trialtest_2d_1()
    test_trialtest_3d_1()
