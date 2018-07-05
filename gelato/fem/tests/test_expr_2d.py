# coding: utf-8

# TODO split the asserts between algebraic and weak formulations ones

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
from gelato.fem.expr import gelatize, normalize
from gelato.fem.expr import normalize_weak_from


# ...
def test_gelatize_2d_1():
    print('============ test_gelatize_2d_1 =============')

    V = FemSpace('V', ldim=2)

    v = TestFunction(V, name='v')
    w = TestFunction(V, name='w')
    c = Constant('c')
    F = Field('F')

    # ... expressions that can be normalized (valid for a weak formulation)
    assert(gelatize(grad(v)) == Tuple(dx(v), dy(v)))
    assert(gelatize(grad(c*v)) == Tuple(c*dx(v), c*dy(v)))
    assert(gelatize(grad(F*v)) == Tuple(F*dx(v) + v*dx(F), F*dy(v) + v*dy(F)))

    assert(gelatize(dot(grad(v), grad(w))) == dx(v)*dx(w) + dy(v)*dy(w))
    # ...

#    expr = grad(F*v)
#    print('> input         >>> {0}'.format(expr))
#    print('> gelatized     >>> {0}'.format(gelatize(expr)))
# ...

# ...
def test_normalize_2d_1():
    print('============ test_normalize_2d_1 =============')

    V = FemSpace('V', ldim=2)
    W = FemSpace('W', ldim=2)

    v = TestFunction(V, name='v')
    u = TrialFunction(W, name='u')
    c = Constant('c')
    F = Field('F')

    Ni, Ni_x, Ni_y = symbols('Ni Ni_x Ni_y')
    Nj, Nj_x, Nj_y = symbols('Nj Nj_x Nj_y')

    bx, by = symbols('bx by')
    b = Tuple(bx, by)

    # ...
    assert(normalize(grad(v), basis={V: 'Ni'}) == Tuple(Ni_x, Ni_y))
    assert(normalize(grad(c*v), basis={V: 'Ni'}) == Tuple(c*Ni_x, c*Ni_y))
    assert(normalize(dot(b, grad(v)), basis={V: 'Ni'}) == Ni_x*bx + Ni_y*by)
    assert(normalize(dot(b, grad(v)) + c*v, basis={V: 'Ni'}) == Ni_x*bx + Ni_y*by + c*Ni)
    assert(normalize(dot(grad(v), grad(u)), basis={V: 'Ni', W: 'Nj'}) == Ni_x*Nj_x + Ni_y*Nj_y)
    assert(normalize(dot(grad(v), grad(u)) + c*v*u, basis={V: 'Ni', W: 'Nj'}) == Ni_x*Nj_x + Ni_y*Nj_y + c*Ni*Nj)
    # ...

#    expr = dot(grad(v), b)
#    print('> input         >>> {0}'.format(expr))
#
#    print('> normal form   >>> {0}'.format(normalize(expr, basis={V: 'Ni'})))
#    print('> normal form   >>> {0}'.format(normalize(dot(grad(v), grad(u)),
#                                                     basis={V: 'Ni', W: 'Nj'})))
# ...


# ...
def test_gelatize_2d_2():
    print('============ test_gelatize_2d_2 =============')

    V = FemSpace('V', ldim=2, is_vector=True, shape=2)

    v = VectorTestFunction(V, name='v')

    assert(gelatize(rot(v)) == -dx(v[1]) + dy(v[0]))
    assert(gelatize(div(v)) == dx(v[0]) + dy(v[1]))

#    expr = div(v)
#    print('> input         >>> {0}'.format(expr))
#    print('> gelatized     >>> {0}'.format(gelatize(expr)))
# ...

# ...
def test_normalize_2d_2():
    print('============ test_normalize_2d_2 =============')

    V = FemSpace('V', ldim=2, is_vector=True, shape=2)

    v = VectorTestFunction(V, name='v')

    v0_x = Symbol('v_x[0]')
    v1_x = Symbol('v_x[1]')
    v0_y = Symbol('v_y[0]')
    v1_y = Symbol('v_y[1]')

#    assert(normalize_weak_from(rot(v)) == -v1_x + v0_y)
#    assert(normalize_weak_from(div(v)) == v0_x + v1_y)

#    expr = div(v)
#    print('> input         >>> {0}'.format(expr))
#    print('> normal form   >>> {0}'.format(normalize_weak_from(expr)))
# ...

# ...
def test_bilinear_form_2d_1():
    print('============ test_bilinear_form_2d_1 =============')

    W = FemSpace('W', ldim=2)
    V = FemSpace('V', ldim=2)

    w = TestFunction(W, name='w')
    v = TrialFunction(V, name='v')

    expr = inner(grad(w), grad(v)) + w*v

    a = BilinearForm(expr, trial_space=V, test_space=W)
    print('> input         >>> {0}'.format(a))
    print('> gelatized     >>> {0}'.format(gelatize(a)))
    print('> normal form   >>> {0}'.format(normalize_weak_from(a)))

    a_expr = normalize_weak_from(a, basis={V: 'Nj', W: 'Ni'})
    print('> basis  form   >>> {0}'.format(a_expr))
    print('')
# ...

# ...
def test_bilinear_form_2d_2():
    print('============ test_bilinear_form_2d_2 =============')

    W = FemSpace('W', ldim=2)
    V = FemSpace('V', ldim=2)

    w = TestFunction(W, name='w')
    v = TrialFunction(V, name='v')

    expr = cross(curl(w), curl(v)) + 0.2 * w * v

    a = BilinearForm(expr, trial_space=V, test_space=W)
    print('> input         >>> {0}'.format(a))
    print('> gelatized     >>> {0}'.format(gelatize(a)))
    print('> normal form   >>> {0}'.format(normalize_weak_from(a)))

    a_expr = normalize_weak_from(a, basis={V: 'Nj', W: 'Ni'})
    print('> basis  form   >>> {0}'.format(a_expr))
    print('')
# ...

# ...
def test_bilinear_form_2d_3():
    print('============ test_bilinear_form_2d_3 =============')

    W = FemSpace('W', ldim=2)
    V = FemSpace('V', ldim=2)

    w = TestFunction(W, name='w')
    v = TrialFunction(V, name='v')

    bx = Constant('bx')
    by = Constant('by')
    b = Tuple(bx, by)

    expr = 0.2 * w * v + dot(b, grad(v)) * w

    a = BilinearForm(expr, trial_space=V, test_space=W)
    print('> input         >>> {0}'.format(a))
    print('> gelatized     >>> {0}'.format(gelatize(a)))
    print('> normal form   >>> {0}'.format(normalize_weak_from(a)))

    a_expr = normalize_weak_from(a, basis={V: 'Nj', W: 'Ni'})
    print('> basis  form   >>> {0}'.format(a_expr))
    print('')
# ...

# ...
def test_bilinear_form_2d_4():
    print('============ test_bilinear_form_2d_4 =============')

    W = FemSpace('W', ldim=2, is_block=True, shape=2)
    V = FemSpace('V', ldim=2, is_block=True, shape=2)

    w = VectorTestFunction(W, name='w')
    v = VectorTrialFunction(V, name='v')

    expr = rot(w) * rot(v) + div(w) * div(v) + 0.2 * dot(w, v)

    a = BilinearForm(expr, trial_space=V, test_space=W)
    print('> input         >>> {0}'.format(a))
    print('> gelatized     >>> {0}'.format(gelatize(a)))
    print('> normal form   >>> {0}'.format(normalize_weak_from(a)))
    print('')
# ...

# ... TODO debug
def test_bilinear_form_2d_10():
    print('============ test_bilinear_form_2d_10 =============')

    W = FemSpace('W', ldim=2)
    V = FemSpace('V', ldim=2)

    w = TestFunction(W, name='w')
    v = TrialFunction(V, name='v')

    a = BilinearForm(inner(grad(w), grad(v)), trial_space=V, test_space=W)
    b = BilinearForm(w*v, trial_space=V, test_space=W)

    c = a + b
    print('> input         >>> {0}'.format(c))
    print('> gelatized     >>> {0}'.format(gelatize(c)))
    print('> normal form   >>> {0}'.format(normalize_weak_from(c)))
    print('')

    v1 = TestFunction(V, name='v1')
    u1 = TrialFunction(W, name='u1')

    d = a(v1, u1) + b
    print('> input         >>> {0}'.format(d))
    print('> gelatized     >>> {0}'.format(gelatize(d)))
    print('> normal form   >>> {0}'.format(normalize_weak_from(d)))
    print('')
# ...


# .....................................................
if __name__ == '__main__':
    test_gelatize_2d_1()
    test_normalize_2d_1()

#    test_gelatize_2d_2()
#    test_normalize_2d_2()

#    test_bilinear_form_2d_1()
#    test_bilinear_form_2d_2()
#    test_bilinear_form_2d_3()
#    test_bilinear_form_2d_4()
##    test_bilinear_form_2d_10()
