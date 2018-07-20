# coding: utf-8

# TODO split the asserts between algebraic and weak formulations ones
# TODO: - __call__ examples are not working anymore

from sympy.core.containers import Tuple
from sympy import symbols
from sympy import Symbol
from sympy import Function
from sympy import pi, cos, sin
from sympy import srepr

from gelato.core import dx, dy, dz
from gelato.core import Constant
from gelato.core import Field
from gelato.core import grad, dot, inner, cross, rot, curl, div
from gelato.core import H1Space
from gelato.core import TestFunction
from gelato.core import VectorTestFunction
from gelato.core import BilinearForm, LinearForm, FunctionForm
from gelato.core import atomize, normalize
from gelato.core import gelatize


# ...
def test_atomize_1d_1():
    print('============ test_atomize_1d_1 =============')

    V = H1Space('V', ldim=1)

    v = TestFunction(V, name='v')
    w = TestFunction(V, name='w')
    c = Constant('c')
    F = Field('F', space=V)
    x = Symbol('x')
    f = Function('f')

    # ...
    assert(atomize(grad(v)) == dx(v))
    assert(atomize(grad(c*v)) == c*dx(v))
    assert(atomize(grad(F*v)) == F*dx(v) + v*dx(F))
    assert(atomize(f(x)*grad(v)) == dx(v)*f(x))

    assert(atomize(dot(grad(v), grad(w))) == dx(v)*dx(w))
    # ...

    # ...
    assert(atomize(grad(v*w)) == w*dx(v) + v*dx(w))
    assert(atomize(div(grad(v*w))) == 2*dx(v)*dx(w) + dx(dx(v))*w + dx(dx(w))*v)
    # ...

#    expr = div(grad(v*w))
#    print('> input         >>> {0}'.format(expr))
#    print('> atomized     >>> {0}'.format(atomize(expr)))
#    print(expr.is_commutative)
# ...

# ...
def test_normalize_1d_1():
    print('============ test_normalize_1d_1 =============')

    V = H1Space('V', ldim=1)
    U = H1Space('U', ldim=1)

    v = TestFunction(V, name='v')
    u = TestFunction(U, name='u')
    c = Constant('c')
    F = Field('F', space=V)

    Ni, Ni_x = symbols('Ni Ni_x')
    Nj, Nj_x = symbols('Nj Nj_x')

    # ...
    assert(normalize(grad(v), basis={V: 'Ni'}) == Ni_x)
    assert(normalize(grad(c*v), basis={V: 'Ni'}) == c*Ni_x)
    assert(normalize(grad(v) + c*v, basis={V: 'Ni'}) == Ni_x + c*Ni)
    assert(normalize(grad(F*v), basis={V: 'Ni'}) == F*Ni_x + Ni*dx(F))

    assert(normalize(dot(grad(v), grad(u)), basis={V: 'Ni', U: 'Nj'}) == Ni_x*Nj_x)
    assert(normalize(dot(grad(v), grad(u)) + c*v*u, basis={V: 'Ni', U: 'Nj'}) == Ni_x*Nj_x + c*Ni*Nj)
    assert(normalize(dot(grad(F*v), grad(u)), basis={V: 'Ni', U: 'Nj'}) == Nj_x*(F*Ni_x + Ni*dx(F)))
    # ...

#    expr = dot(grad(F*v), grad(u))
#    print('> input         >>> {0}'.format(expr))

#    print('> normal form   >>> {0}'.format(normalize(expr, basis={V: 'Ni'})))
#    print('> normal form   >>> {0}'.format(normalize(expr, basis={V: 'Ni', U: 'Nj'})))
# ...

# ...
def test_gelatize_1d_1():
    print('============ test_gelatize_1d_1 =============')

    V = H1Space('V', ldim=1)
    U = H1Space('U', ldim=1)

    v = TestFunction(V, name='v')
    u = TestFunction(U, name='u')
    c = Constant('c')
    F = Field('F', space=V)

    Ni, Ni_x = symbols('Ni Ni_x')
    Nj, Nj_x = symbols('Nj Nj_x')

    # ...
    a = BilinearForm((v,u), dot(grad(v), grad(u)))
    assert(gelatize(a, basis={V: 'Ni', U: 'Nj'}) == Ni_x*Nj_x)
    # ...

    # ...
    a = BilinearForm((v,u), dot(grad(v), grad(u)) + c*v*u)
    assert(gelatize(a, basis={V: 'Ni', U: 'Nj'}) == Ni_x*Nj_x + c*Ni*Nj)
    # ...

    # ...
    a = BilinearForm((v,u), dot(grad(v), grad(u)) + F*v*u)
    assert(gelatize(a, basis={V: 'Ni', U: 'Nj'}) == Ni_x*Nj_x + F*Ni*Nj)
    # ...

    # ...
    a = BilinearForm((v,u), dot(grad(F*v), grad(u)))
    assert(gelatize(a, basis={V: 'Ni', U: 'Nj'}) == F*Ni_x*Nj_x + Ni*Nj_x*dx(F))
    # ...
# ...

# ...
def test_bilinear_form_1d_1():
    print('============ test_bilinear_form_1d_1 =============')

    U = H1Space('U', ldim=1)
    V = H1Space('V', ldim=1)

    v = TestFunction(V, name='v')
    u = TestFunction(U, name='u')

    expr = inner(grad(v), grad(u))

    a = BilinearForm((v,u), expr)
    print('> input      >>> {0}'.format(a))

    a_expr = gelatize(a, basis={V: 'Nj', U: 'Ni'})
    print('> gelatized  >>> {0}'.format(a_expr))
    print('')
# ...

# ...
def test_bilinear_form_1d_2():
    print('============ test_bilinear_form_1d_2 =============')

    U = H1Space('U', ldim=1)
    V = H1Space('V', ldim=1)

    v = TestFunction(V, name='v')
    u = TestFunction(U, name='u')

    b = Constant('b')

    expr = inner(grad(b*u), grad(v))

    a = BilinearForm((v,u), expr)
    print('> input         >>> {0}'.format(a))

    a_expr = gelatize(a, basis={V: 'Nj', U: 'Ni'})
    print('> gelatized  >>> {0}'.format(a_expr))
    print('')
# ...

# ...
def test_bilinear_form_1d_3():
    print('============ test_bilinear_form_1d_3 =============')

    U = H1Space('U', ldim=1)
    V = H1Space('V', ldim=1)

    v = TestFunction(V, name='v')
    u = TestFunction(U, name='u')

    F = Field('F', space=V)
    expr = inner(grad(u), grad(v)) + F*u*v

    a = BilinearForm((v,u), expr)
    print('> input         >>> {0}'.format(a))

    a_expr = gelatize(a, basis={V: 'Nj', U: 'Ni'})
    print('> basis  form   >>> {0}'.format(a_expr))
    print('')
# ...

## ... TODO debug
#def test_bilinear_form_1d_4():
#    print('============ test_bilinear_form_1d_4 =============')
#
#    U = H1Space('U', ldim=1)
#    V = H1Space('V', ldim=1)
#
#    v = TestFunction(V, name='v')
#    u = TestFunction(U, name='u')
#
#    F = Field('F', space=V)
#
#    expr = inner(grad(F*u), grad(v))
#
#    a = BilinearForm((v,u), expr)
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
#    U = H1Space('U', ldim=1)
#    V = H1Space('V', ldim=1)
#
#    v = TestFunction(V, name='v')
#    u = TestFunction(U, name='u')
#
#    expr = dx(dx(v))*dx(dx(dx(u))) + u*v
#
#    a = BilinearForm((v,u), expr)
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
#    U = H1Space('U', ldim=1)
#    V = H1Space('V', ldim=1)
#
#    v = TestFunction(V, name='v')
#    u = TestFunction(U, name='u')
#
#    F = Field('F', space=V)
#
#    expr = inner(grad(dx(F)*v), grad(u)) + u*v
#
#    a = BilinearForm((v,u), expr)
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
#    U = H1Space('U', ldim=1)
#    V = H1Space('V', ldim=1)
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
#    U = H1Space('U', ldim=1)
#    V = H1Space('V', ldim=1)
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

# ...
def test_bilinear_form_1d_10():
    print('============ test_bilinear_form_1d_10 =============')

    U = H1Space('U', ldim=1)
    V = H1Space('V', ldim=1)

    u = TestFunction(U, name='u')
    v = TestFunction(V, name='v')

    u1 = TestFunction(U, name='u1')
    v1 = TestFunction(V, name='v1')

    Ni, Ni_x = symbols('Ni Ni_x')
    Nj, Nj_x = symbols('Nj Nj_x')

    c1 = Symbol('c1')
    c2 = Symbol('c2')

    a = BilinearForm((v,u), inner(grad(u), grad(v)))
    b = BilinearForm((v,u), u*v)
    adv = BilinearForm((v,u), dx(u)*v)

    # ...
    expected = Ni*Nj + Ni_x*Nj_x
    assert(gelatize(a + b, basis={V: 'Nj', U: 'Ni'}) == expected)
    # ...

    # ...
    expected = 2*Ni_x*Nj_x
    assert(gelatize(2 * a, basis={V: 'Nj', U: 'Ni'}) == expected)
    # ...

    # ...
    expected = c1*Ni_x*Nj_x
    assert(gelatize(c1*a, basis={V: 'Nj', U: 'Ni'}) == expected)
    # ...

    # ...
    expected = c2*Ni*Nj + c1*Ni_x*Nj_x
    assert(gelatize(c1*a + c2*b, basis={V: 'Nj', U: 'Ni'}) == expected)
    # ...

    # ...
    expected = Ni_x*Nj_x*c1 + c2*(Ni*Nj + Ni_x*Nj)
    assert(gelatize(c1*a  + c2*(b + adv), basis={V: 'Nj', U: 'Ni'}) == expected)
    # ...

    # ...
    assert(gelatize(a(u1, v1), basis={V: 'Nj', U: 'Ni'}) == gelatize(a(v1, u1), basis={V: 'Nj', U: 'Ni'}))
    # ...

#    # ... TODO debug
#    expected = Ni_x*Nj
#    assert(gelatize(adv(v1, u1), basis={V: 'Nj', U: 'Ni'}) == expected)
#
#    expected = Nj_x*Ni
#    assert(gelatize(adv(u1, v1), basis={V: 'Nj', U: 'Ni'}) == expected)
#    # ...

#    expr = adv(v1, u1)
#    print('> input      >>> {0}'.format(expr))
#    print('> gelatized  >>> {0}'.format(gelatize(expr, basis={V: 'Nj', U: 'Ni'}) ))
#    print('')
# ...

# ...
def test_linear_form_1d_1():
    print('============ test_linear_form_1d_1 =============')

    V = H1Space('V', ldim=1)

    v = TestFunction(V, name='v')

    x = V.coordinates
    f = Function('f')

    expr = cos(2*pi*x)*v
#    expr = f(x)*v

    a = LinearForm(v, expr)
    print('> input      >>> {0}'.format(a))

    a_expr = gelatize(a, basis={V: 'Nj'})
    print('> gelatized  >>> {0}'.format(a_expr))
    print('')
# ...

# ...
def test_linear_form_1d_10():
    print('============ test_linear_form_1d_10 =============')

    V = H1Space('V', ldim=1)

    v = TestFunction(V, name='v')

    x = V.coordinates
    f = Function('f')

    Ni, Ni_x, Ni_xx = symbols('Ni Ni_x Ni_xx')

    c1 = Symbol('c1')
    c2 = Symbol('c2')

    # ...
    expected = cos(2*pi*x)*Ni
    assert(gelatize(LinearForm(v, cos(2*pi*x)*v),
                    basis={V: 'Ni'}) == expected)
    # ...

    # ...
    expected = f(x)*Ni
    assert(gelatize(LinearForm(v, f(x)*v),
                    basis={V: 'Ni'}) == expected)
    # ...

    # ...
    expected = cos(2*pi*x)*Ni_x
    assert(gelatize(LinearForm(v, cos(2*pi*x)*dx(v)),
                    basis={V: 'Ni'}) == expected)
    # ...

    # ...
    expected = f(x)*Ni_xx
    assert(gelatize(LinearForm(v, f(x)*dx(dx(v))),
                    basis={V: 'Ni'}) == expected)
    # ...

    expr = LinearForm(v, cos(2*pi*x)*dx(v))
    print('> input      >>> {0}'.format(expr))
    print('> gelatized  >>> {0}'.format(gelatize(expr, basis={V: 'Ni'}) ))
    print('')
# ...

# ...
def test_function_form_1d_10():
    print('============ test_function_form_1d_10 =============')

    V = H1Space('V', ldim=1)

    F = Field('F', space=V)

    x = V.coordinates
    f = Function('f')

    Ni, Ni_x, Ni_xx = symbols('Ni Ni_x Ni_xx')

    c1 = Symbol('c1')
    c2 = Symbol('c2')

    # ...
    expected = -2*pi*sin(2*pi*x)
    assert(gelatize(FunctionForm(grad(cos(2*pi*x)), coordinates=[x])) == expected)
    # ...

    # ...
    expected = F-cos(2*pi*x)
    assert(gelatize(FunctionForm(F-cos(2*pi*x))) == expected)
    # ...

    # ...
    expected = (F-cos(2*pi*x))**2
    assert(gelatize(FunctionForm((F-cos(2*pi*x))**2)) == expected)
    # ...

    # ...
    expected = dx(F) + 2*pi*sin(2*pi*x)
    assert(gelatize(FunctionForm(grad(F-cos(2*pi*x)))) == expected)
    # ...

    # ...
    expected = (dx(F) + 2*pi*sin(2*pi*x))**2
    assert(gelatize(FunctionForm((grad(F-cos(2*pi*x)))**2)) == expected)
    # ...

#    expr = FunctionForm()
#    print('> input      >>> {0}'.format(expr))
#    print('> gelatized  >>> {0}'.format(gelatize(expr) ))
#    print('')
# ...

# .....................................................
if __name__ == '__main__':
    test_atomize_1d_1()
    test_normalize_1d_1()
    test_gelatize_1d_1()

    test_bilinear_form_1d_1()
    test_bilinear_form_1d_2()
    test_bilinear_form_1d_3()

##    test_bilinear_form_1d_4()
##    test_bilinear_form_1d_5()
##    test_bilinear_form_1d_6()
##    test_bilinear_form_1d_7()
##    test_bilinear_form_1d_8()

    test_bilinear_form_1d_10()
    test_linear_form_1d_10()
    test_function_form_1d_10()

    test_linear_form_1d_1()
