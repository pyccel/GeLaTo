# coding: utf-8

from sympy import Symbol
from sympy.core.containers import Tuple
from sympy import symbols

from gelato.calculus import dx, dy, dz
from gelato.calculus import Constant
from gelato.calculus import Field
from gelato.calculus import grad, dot, inner, cross, rot, curl, div

from gelato.fem.core import FemSpace
from gelato.fem.core import TestFunction
from gelato.fem.core import TrialFunction
from gelato.fem.core import VectorTestFunction
from gelato.fem.core import VectorTrialFunction
from gelato.fem.expr import BilinearForm
from gelato.fem.expr import gelatize, normalize_weak_from


# ...
def test_gelatize_3d_1():
    print('============ test_gelatize_3d_1 =============')

    V = FemSpace('V', ldim=3, is_vector=True, shape=3)

    v = VectorTestFunction(V, name='v')

    assert(gelatize(curl(v)) == Tuple( dy(v[2]) - dz(v[1]),
                                      -dx(v[2]) + dz(v[0]),
                                       dx(v[1]) - dy(v[0])))
    assert(gelatize(div(v)) == dx(v[0]) + dy(v[1]) + dz(v[2]))

#    expr = curl(v)
#    print('> input         >>> {0}'.format(expr))
#    print('> gelatized     >>> {0}'.format(gelatize(expr)))
# ...

# ...
def test_normalize_3d_1():
    print('============ test_normalize_3d_1 =============')

    V = FemSpace('V', ldim=3, is_vector=True, shape=3)

    v = VectorTestFunction(V, name='v')

    v0_x = Symbol('v_x[0]')
    v1_x = Symbol('v_x[1]')
    v2_x = Symbol('v_x[2]')
    v0_y = Symbol('v_y[0]')
    v1_y = Symbol('v_y[1]')
    v2_y = Symbol('v_y[2]')
    v0_z = Symbol('v_z[0]')
    v1_z = Symbol('v_z[1]')
    v2_z = Symbol('v_z[2]')

#    assert(normalize_weak_from(curl(v)) == Tuple( v2_y - v1_z,
#                                                 -v2_x + v0_z,
#                                                  v1_x - v0_y))
#    assert(normalize_weak_from(div(v)) == v0_x + v1_y + v2_z)

#    expr = curl(v)
#    print('> input         >>> {0}'.format(expr))
#    print('> normal form   >>> {0}'.format(normalize_weak_from(expr)))
# ...

# ...
def test_bilinear_form_3d_1():
    print('============ test_bilinear_form_3d_1 =============')

    W = FemSpace('W', ldim=3)
    V = FemSpace('V', ldim=3)

    w = TestFunction(W, name='w')
    v = TrialFunction(V, name='v')

    expr = inner(grad(w), grad(v))

    a = BilinearForm(expr, trial_space=V, test_space=W)
    print('> input         >>> {0}'.format(a))
    print('> gelatized     >>> {0}'.format(gelatize(a)))
    print('> normal form   >>> {0}'.format(normalize_weak_from(a)))

    a_expr = normalize_weak_from(a, basis={V: 'Nj', W: 'Ni'})
    print('> basis  form   >>> {0}'.format(a_expr))
    print('')
# ...

# ...
def test_bilinear_form_3d_2():
    print('============ test_bilinear_form_3d_2 =============')

    W = FemSpace('W', ldim=3, is_vector=True, shape=3)
    V = FemSpace('V', ldim=3, is_vector=True, shape=3)

    w = VectorTestFunction(W, name='w')
    v = VectorTrialFunction(V, name='v')

    expr = div(w) * div(v) + 0.2 * dot(w, v)

    a = BilinearForm(expr, trial_space=V, test_space=W)
    print('> input         >>> {0}'.format(a))
    print('> gelatized     >>> {0}'.format(gelatize(a)))
    print('> normal form   >>> {0}'.format(normalize_weak_from(a)))
    print('')
# ...

# ...
def test_bilinear_form_3d_3():
    print('============ test_bilinear_form_3d_3 =============')

    W = FemSpace('W', ldim=3, is_vector=True, shape=3)
    V = FemSpace('V', ldim=3, is_vector=True, shape=3)

    w = VectorTestFunction(W, name='w')
    v = VectorTrialFunction(V, name='v')

    expr = dot(curl(w), curl(v)) + 0.2 * dot(w, v)

    a = BilinearForm(expr, trial_space=V, test_space=W)
    print('> input         >>> {0}'.format(a))
    print('> gelatized     >>> {0}'.format(gelatize(a)))
    print('> normal form   >>> {0}'.format(normalize_weak_from(a)))
    print('')
# ...

# ...
def test_bilinear_form_3d_4():
    print('============ test_bilinear_form_3d_4 =============')

    W = FemSpace('W', ldim=3, is_vector=True, shape=3)
    V = FemSpace('V', ldim=3, is_vector=True, shape=3)

    w = VectorTestFunction(W, name='w')
    v = VectorTrialFunction(V, name='v')

    bx = Constant('bx')
    by = Constant('by')
    bz = Constant('bz')
    b = Tuple(bx, by, bz)

    expr = dot(curl(cross(b,w)), curl(cross(b,v))) + 0.2 * dot(w, v)

    a = BilinearForm(expr, trial_space=V, test_space=W)
    print('> input         >>> {0}'.format(a))
    print('> gelatized     >>> {0}'.format(gelatize(a)))
    print('> normal form   >>> {0}'.format(normalize_weak_from(a)))
    print('')
# ...

# ...
def test_bilinear_form_3d_5():
    print('============ test_bilinear_form_3d_5 =============')

    W = FemSpace('W', ldim=3, is_vector=True, shape=3)
    V = FemSpace('V', ldim=3, is_vector=True, shape=3)

    w = VectorTestFunction(W, name='w')
    v = VectorTrialFunction(V, name='v')

    bx = Constant('bx')
    by = Constant('by')
    bz = Constant('bz')
    b = Tuple(bx, by, bz)

    c0,c1,c2 = symbols('c0 c1 c2')

    expr = c0 * dot(w, v) - c1 * div(w) * div(v) + c2 *dot(curl(cross(b,w)), curl(cross(b,v)))

    a = BilinearForm(expr, trial_space=V, test_space=W)
    print('> input         >>> {0}'.format(a))
    print('> gelatized     >>> {0}'.format(gelatize(a)))
    print('> normal form   >>> {0}'.format(normalize_weak_from(a)))
    print('')
# ...

# ... TODO debug
def test_bilinear_form_3d_10():
    print('============ test_bilinear_form_3d_10 =============')

    W = FemSpace('W', ldim=3)
    V = FemSpace('V', ldim=3)

    w = TestFunction(W, name='w')
    v = TrialFunction(V, name='v')

    a = BilinearForm(inner(grad(w), grad(v)), trial_space=V, test_space=W)
    b = BilinearForm(w*v, trial_space=V, test_space=W)
    c = a + b
    print('> input         >>> {0}'.format(c))
    print('> gelatized     >>> {0}'.format(gelatize(c)))
    print('> normal form   >>> {0}'.format(normalize_weak_from(c)))
    print('')
# ...


# .....................................................
if __name__ == '__main__':
    test_gelatize_3d_1()
    test_normalize_3d_1()

#    test_bilinear_form_3d_1()
#    test_bilinear_form_3d_2()
#    test_bilinear_form_3d_3()
#    test_bilinear_form_3d_4()
#    test_bilinear_form_3d_5()
##    test_bilinear_form_3d_10()