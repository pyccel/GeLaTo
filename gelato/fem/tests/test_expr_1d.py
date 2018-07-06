# coding: utf-8

# TODO split the asserts between algebraic and weak formulations ones

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
from gelato.fem.expr import atomize, normalize
from gelato.fem.expr import normalize_weak_from


# ...
def test_atomize_1d_1():
    print('============ test_atomize_1d_1 =============')

    V = FemSpace('V', ldim=1)

    v = TestFunction(V, name='v')
    w = TestFunction(V, name='w')
    c = Constant('c')
    F = Field('F')

    # ...
    assert(atomize(grad(v)) == dx(v))
    assert(atomize(grad(c*v)) == c*dx(v))
    assert(atomize(grad(F*v)) == F*dx(v) + v*dx(F))

    assert(atomize(dot(grad(v), grad(w))) == dx(v)*dx(w))
    # ...

    # ...
    assert(atomize(grad(v*w)) == w*dx(v) + v*dx(w))
    assert(atomize(div(grad(v*w))) == v*dx(dx(w)) + 2*dx(v)*dx(w) + dx(dx(v))*w)
    # ...

#    expr = dot(grad(v), grad(w))
#    print('> input         >>> {0}'.format(expr))
#    print('> atomized     >>> {0}'.format(atomize(expr)))
# ...

# ...
def test_normalize_1d_1():
    print('============ test_normalize_1d_1 =============')

    V = FemSpace('V', ldim=1)
    W = FemSpace('W', ldim=1)

    v = TestFunction(V, name='v')
    u = TrialFunction(W, name='u')
    c = Constant('c')
    F = Field('F')

    Ni, Ni_x = symbols('Ni Ni_x')
    Nj, Nj_x = symbols('Nj Nj_x')

    # ...
    assert(normalize(grad(v), basis={V: 'Ni'}) == Ni_x)
    assert(normalize(grad(c*v), basis={V: 'Ni'}) == c*Ni_x)
    assert(normalize(grad(v) + c*v, basis={V: 'Ni'}) == Ni_x + c*Ni)
    assert(normalize(dot(grad(v), grad(u)), basis={V: 'Ni', W: 'Nj'}) == Ni_x*Nj_x)
    assert(normalize(dot(grad(v), grad(u)) + c*v*u, basis={V: 'Ni', W: 'Nj'}) == Ni_x*Nj_x + c*Ni*Nj)
    # ...

#    expr = dot(grad(v), grad(u))
#    print('> input         >>> {0}'.format(expr))

#    print('> normal form   >>> {0}'.format(normalize(expr, basis={V: 'Ni'})))
#    print('> normal form   >>> {0}'.format(normalize(dot(grad(v), grad(u)),
#                                                     basis={V: 'Ni', W: 'Nj'})))
# ...

# ...
def test_bilinear_form_1d_1():
    print('============ test_bilinear_form_1d_1 =============')

    W = FemSpace('W', ldim=1)
    V = FemSpace('V', ldim=1)

    w = TestFunction(W, name='w')
    v = TrialFunction(V, name='v')

    expr = inner(grad(w), grad(v))

    a = BilinearForm(expr, trial_space=V, test_space=W)
    print('> input         >>> {0}'.format(a))
    print('> atomized     >>> {0}'.format(atomize(a)))
    print('> normal form   >>> {0}'.format(normalize_weak_from(a)))

    a_expr = normalize_weak_from(a, basis={V: 'Nj', W: 'Ni'})
    print('> basis  form   >>> {0}'.format(a_expr))
    print('')
# ...

# ...
def test_bilinear_form_1d_2():
    print('============ test_bilinear_form_1d_2 =============')

    W = FemSpace('W', ldim=1)
    V = FemSpace('V', ldim=1)

    w = TestFunction(W, name='w')
    v = TrialFunction(V, name='v')

    b = Constant('b')

    expr = inner(grad(b*w), grad(v))

    a = BilinearForm(expr, trial_space=V, test_space=W)
    print('> input         >>> {0}'.format(a))
    print('> atomized     >>> {0}'.format(atomize(a)))
    print('> normal form   >>> {0}'.format(normalize_weak_from(a)))

    a_expr = normalize_weak_from(a, basis={V: 'Nj', W: 'Ni'})
    print('> basis  form   >>> {0}'.format(a_expr))
    print('')
# ...

# ...
def test_bilinear_form_1d_3():
    print('============ test_bilinear_form_1d_3 =============')

    W = FemSpace('W', ldim=1)
    V = FemSpace('V', ldim=1)

    w = TestFunction(W, name='w')
    v = TrialFunction(V, name='v')

    F = Field('F')
    expr = inner(grad(w), grad(v)) + F*w*v

    a = BilinearForm(expr, trial_space=V, test_space=W)
    print('> input         >>> {0}'.format(a))
    print('> atomized     >>> {0}'.format(atomize(a)))
    print('> normal form   >>> {0}'.format(normalize_weak_from(a)))

    a_expr = normalize_weak_from(a, basis={V: 'Nj', W: 'Ni'})
    print('> basis  form   >>> {0}'.format(a_expr))
    print('')
# ...

# ...
def test_bilinear_form_1d_4():
    print('============ test_bilinear_form_1d_4 =============')

    W = FemSpace('W', ldim=1)
    V = FemSpace('V', ldim=1)

    w = TestFunction(W, name='w')
    v = TrialFunction(V, name='v')

    F = Field('F')

    expr = inner(grad(F*w), grad(v))

    a = BilinearForm(expr, trial_space=V, test_space=W)
    print('> input         >>> {0}'.format(a))
    print('> atomized     >>> {0}'.format(atomize(a)))
    print('> normal form   >>> {0}'.format(normalize_weak_from(a)))

    a_expr = normalize_weak_from(a, basis={V: 'Nj', W: 'Ni'})
    print('> basis  form   >>> {0}'.format(a_expr))
    print('')
# ...

# ...
def test_bilinear_form_1d_5():
    print('============ test_bilinear_form_1d_5 =============')

    W = FemSpace('W', ldim=1)
    V = FemSpace('V', ldim=1)

    w = TestFunction(W, name='w')
    v = TrialFunction(V, name='v')

    expr = dx(dx(v))*dx(dx(dx(w))) + w*v

    a = BilinearForm(expr, trial_space=V, test_space=W)
    print('> input         >>> {0}'.format(a))
    print('> atomized     >>> {0}'.format(atomize(a)))
    print('> normal form   >>> {0}'.format(normalize_weak_from(a)))

    a_expr = normalize_weak_from(a, basis={V: 'Nj', W: 'Ni'})
    print('> basis  form   >>> {0}'.format(a_expr))
    print('')
# ...

# ...
def test_bilinear_form_1d_6():
    print('============ test_bilinear_form_1d_6 =============')

    W = FemSpace('W', ldim=1)
    V = FemSpace('V', ldim=1)

    w = TestFunction(W, name='w')
    v = TrialFunction(V, name='v')

    F = Field('F')

    expr = inner(grad(dx(F)*v), grad(w)) + w*v

    a = BilinearForm(expr, trial_space=V, test_space=W)
    print('> input         >>> {0}'.format(a))
    print('> atomized     >>> {0}'.format(atomize(a)))
    print('> normal form   >>> {0}'.format(normalize_weak_from(a)))

    a_expr = normalize_weak_from(a, basis={V: 'Nj', W: 'Ni'})
    print('> basis  form   >>> {0}'.format(a_expr))
    print('')
# ...

# ... TODO debug
def test_bilinear_form_1d_7():
    print('============ test_bilinear_form_1d_7 =============')

    W = FemSpace('W', ldim=1)
    V = FemSpace('V', ldim=1)

    w1 = TestFunction(W,  name='w1')
    w2 = TestFunction(W,  name='w2')
    v1 = TrialFunction(V, name='v1')
    v2 = TrialFunction(V, name='v2')
    t1 = TestFunction(W,  name='t1')
    t2 = TestFunction(W,  name='t2')

    a = BilinearForm(inner(grad(w1), grad(v1)) + w2*v2, trial_space=V, test_space=W)
    b = BilinearForm(w1*v2, trial_space=V, test_space=W)

    ls = [a + b, a((t1,t2), (v1,v2)) + b(t1, v2)]
    for c in ls:
        print('> input         >>> {0}'.format(c))
        print('> atomized     >>> {0}'.format(atomize(c)))
        print('> normal form   >>> {0}'.format(normalize_weak_from(c)))
        print('')
# ...

# ... TODO debug
def test_bilinear_form_1d_8():
    print('============ test_bilinear_form_1d_8 =============')

    W = FemSpace('W', ldim=1)
    V = FemSpace('V', ldim=1)

    w1 = TestFunction(W,  name='w1')
    w2 = TestFunction(W,  name='w2')
    v1 = TrialFunction(V, name='v1')
    v2 = TrialFunction(V, name='v2')
    t1 = TestFunction(W,  name='t1')
    t2 = TestFunction(W,  name='t2')

    a = BilinearForm(w1*v2, trial_space=V, test_space=W)

    ls = [a(t1, v2), a(v2, t1)]
    for c in ls:
        print('> input         >>> {0}'.format(c))
        print('> atomized     >>> {0}'.format(atomize(c)))
        print('> normal form   >>> {0}'.format(normalize_weak_from(c)))
        print('')
# ...

# ... TODO debug
def test_bilinear_form_1d_10():
    print('============ test_bilinear_form_1d_10 =============')

    W = FemSpace('W', ldim=1)
    V = FemSpace('V', ldim=1)

    w = TestFunction(W, name='w')
    v = TrialFunction(V, name='v')

    a = BilinearForm(inner(grad(w), grad(v)), trial_space=V, test_space=W)
    b = BilinearForm(w*v, trial_space=V, test_space=W)
    c = a + b
    print('> input         >>> {0}'.format(c))
    print('> atomized     >>> {0}'.format(atomize(c)))
    print('> normal form   >>> {0}'.format(normalize_weak_from(c)))
    print('')
# ...

# .....................................................
if __name__ == '__main__':
    test_atomize_1d_1()
    test_normalize_1d_1()

#    test_bilinear_form_1d_1()
#    test_bilinear_form_1d_2()
#    test_bilinear_form_1d_3()
#    test_bilinear_form_1d_4()
#    test_bilinear_form_1d_5()
#    test_bilinear_form_1d_6()
##    test_bilinear_form_1d_7()
##    test_bilinear_form_1d_8()
##    test_bilinear_form_1d_10()
