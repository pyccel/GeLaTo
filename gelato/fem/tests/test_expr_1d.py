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
from gelato.fem.expr import gelatize


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
    U = FemSpace('U', ldim=1)

    v = TestFunction(V, name='v')
    u = TestFunction(U, name='u')
    c = Constant('c')
    F = Field('F')

    Ni, Ni_x = symbols('Ni Ni_x')
    Nj, Nj_x = symbols('Nj Nj_x')

    # ...
    assert(normalize(grad(v), basis={V: 'Ni'}) == Ni_x)
    assert(normalize(grad(c*v), basis={V: 'Ni'}) == c*Ni_x)
    assert(normalize(grad(v) + c*v, basis={V: 'Ni'}) == Ni_x + c*Ni)

    assert(normalize(dot(grad(v), grad(u)), basis={V: 'Ni', U: 'Nj'}) == Ni_x*Nj_x)
    assert(normalize(dot(grad(v), grad(u)) + c*v*u, basis={V: 'Ni', U: 'Nj'}) == Ni_x*Nj_x + c*Ni*Nj)
    # ...

#    expr = dot(grad(v), grad(u))
#    print('> input         >>> {0}'.format(expr))

#    print('> normal form   >>> {0}'.format(normalize(expr, basis={V: 'Ni'})))
#    print('> normal form   >>> {0}'.format(normalize(dot(grad(v), grad(u)),
#                                                     basis={V: 'Ni', U: 'Nj'})))
# ...

# ...
def test_bilinear_form_1d_1():
    print('============ test_bilinear_form_1d_1 =============')

    U = FemSpace('U', ldim=1)
    V = FemSpace('V', ldim=1)

    v = TestFunction(V, name='v')
    u = TestFunction(U, name='u')

    expr = inner(grad(v), grad(u))

    a = BilinearForm(expr, trial_space=V, test_space=U)
    print('> input      >>> {0}'.format(a))

    a_expr = gelatize(a, basis={V: 'Nj', U: 'Ni'})
    print('> gelatized  >>> {0}'.format(a_expr))
    print('')
# ...

# ...
def test_bilinear_form_1d_2():
    print('============ test_bilinear_form_1d_2 =============')

    U = FemSpace('U', ldim=1)
    V = FemSpace('V', ldim=1)

    v = TestFunction(V, name='v')
    u = TestFunction(U, name='u')

    b = Constant('b')

    expr = inner(grad(b*u), grad(v))

    a = BilinearForm(expr, trial_space=V, test_space=U)
    print('> input         >>> {0}'.format(a))

    a_expr = gelatize(a, basis={V: 'Nj', U: 'Ni'})
    print('> gelatized  >>> {0}'.format(a_expr))
    print('')
# ...

# ...
def test_bilinear_form_1d_3():
    print('============ test_bilinear_form_1d_3 =============')

    U = FemSpace('U', ldim=1)
    V = FemSpace('V', ldim=1)

    v = TestFunction(V, name='v')
    u = TestFunction(U, name='u')

    F = Field('F')
    expr = inner(grad(u), grad(v)) + F*u*v

    a = BilinearForm(expr, trial_space=V, test_space=U)
    print('> input         >>> {0}'.format(a))

    a_expr = gelatize(a, basis={V: 'Nj', U: 'Ni'})
    print('> basis  form   >>> {0}'.format(a_expr))
    print('')
# ...

## ... TODO debug
#def test_bilinear_form_1d_4():
#    print('============ test_bilinear_form_1d_4 =============')
#
#    U = FemSpace('U', ldim=1)
#    V = FemSpace('V', ldim=1)
#
#    v = TestFunction(V, name='v')
#    u = TestFunction(U, name='u')
#
#    F = Field('F')
#
#    expr = inner(grad(F*u), grad(v))
#
#    a = BilinearForm(expr, trial_space=V, test_space=U)
#    print('> input         >>> {0}'.format(a))
#
#    a_expr = gelatize(a, basis={V: 'Nj', U: 'Ni'})
#    print('> basis  form   >>> {0}'.format(a_expr))
#    print('')
## ...
#
## ... TODO debug
#def test_bilinear_form_1d_5():
#    print('============ test_bilinear_form_1d_5 =============')
#
#    U = FemSpace('U', ldim=1)
#    V = FemSpace('V', ldim=1)
#
#    v = TestFunction(V, name='v')
#    u = TestFunction(U, name='u')
#
#    expr = dx(dx(v))*dx(dx(dx(u))) + u*v
#
#    a = BilinearForm(expr, trial_space=V, test_space=U)
#    print('> input         >>> {0}'.format(a))
#
#    a_expr = gelatize(a, basis={V: 'Nj', U: 'Ni'})
#    print('> basis  form   >>> {0}'.format(a_expr))
#    print('')
## ...
#
## ... TODO debug
#def test_bilinear_form_1d_6():
#    print('============ test_bilinear_form_1d_6 =============')
#
#    U = FemSpace('U', ldim=1)
#    V = FemSpace('V', ldim=1)
#
#    v = TestFunction(V, name='v')
#    u = TestFunction(U, name='u')
#
#    F = Field('F')
#
#    expr = inner(grad(dx(F)*v), grad(u)) + u*v
#
#    a = BilinearForm(expr, trial_space=V, test_space=U)
#    print('> input         >>> {0}'.format(a))
#
#    a_expr = gelatize(a, basis={V: 'Nj', U: 'Ni'})
#    print('> basis  form   >>> {0}'.format(a_expr))
#    print('')
## ...
#
## ... TODO debug
#def test_bilinear_form_1d_7():
#    print('============ test_bilinear_form_1d_7 =============')
#
#    U = FemSpace('U', ldim=1)
#    V = FemSpace('V', ldim=1)
#
#    w1 = TestFunction(U,  name='w1')
#    w2 = TestFunction(U,  name='w2')
#    v1 = TestFunction(V, name='v1')
#    v2 = TestFunction(V, name='v2')
#    t1 = TestFunction(U,  name='t1')
#    t2 = TestFunction(U,  name='t2')
#
#    a = BilinearForm(inner(grad(w1), grad(v1)) + w2*v2, trial_space=V, test_space=U)
#    b = BilinearForm(w1*v2, trial_space=V, test_space=U)
#
#    ls = [a + b, a((t1,t2), (v1,v2)) + b(t1, v2)]
#    for c in ls:
#        print('> input         >>> {0}'.format(c))
#        print('> atomized     >>> {0}'.format(atomize(c)))
#        print('> normal form   >>> {0}'.format(gelatize(c)))
#        print('')
## ...
#
## ... TODO debug
#def test_bilinear_form_1d_8():
#    print('============ test_bilinear_form_1d_8 =============')
#
#    U = FemSpace('U', ldim=1)
#    V = FemSpace('V', ldim=1)
#
#    w1 = TestFunction(U,  name='w1')
#    w2 = TestFunction(U,  name='w2')
#    v1 = TestFunction(V, name='v1')
#    v2 = TestFunction(V, name='v2')
#    t1 = TestFunction(U,  name='t1')
#    t2 = TestFunction(U,  name='t2')
#
#    a = BilinearForm(w1*v2, trial_space=V, test_space=U)
#
#    ls = [a(t1, v2), a(v2, t1)]
#    for c in ls:
#        print('> input         >>> {0}'.format(c))
#        print('> atomized     >>> {0}'.format(atomize(c)))
#        print('> normal form   >>> {0}'.format(gelatize(c)))
#        print('')
## ...
#
## ... TODO debug
#def test_bilinear_form_1d_10():
#    print('============ test_bilinear_form_1d_10 =============')
#
#    U = FemSpace('U', ldim=1)
#    V = FemSpace('V', ldim=1)
#
#    w = TestFunction(U, name='w')
#    v = TestFunction(V, name='v')
#
#    a = BilinearForm(inner(grad(w), grad(v)), trial_space=V, test_space=U)
#    b = BilinearForm(w*v, trial_space=V, test_space=U)
#    c = a + b
#    print('> input         >>> {0}'.format(c))
#    print('> atomized     >>> {0}'.format(atomize(c)))
#    print('> normal form   >>> {0}'.format(gelatize(c)))
#    print('')
## ...

# .....................................................
if __name__ == '__main__':
    test_atomize_1d_1()
    test_normalize_1d_1()

    test_bilinear_form_1d_1()
    test_bilinear_form_1d_2()
    test_bilinear_form_1d_3()

#    test_bilinear_form_1d_4()
#    test_bilinear_form_1d_5()
#    test_bilinear_form_1d_6()
##    test_bilinear_form_1d_7()
##    test_bilinear_form_1d_8()
##    test_bilinear_form_1d_10()
