# coding: utf-8

# TODO - split the asserts between algebraic and weak formulations ones
#      - add assert for grad in vector case

from sympy import Symbol
from sympy.core.containers import Tuple
from sympy import symbols
from sympy import IndexedBase
from sympy import Matrix

from gelato.calculus import dx, dy, dz
from gelato.calculus import Constant
from gelato.calculus import Field
from gelato.calculus import grad, dot, inner, cross, rot, curl, div

from gelato.fem.core import BasicSobolevSpace
from gelato.fem.core import TestFunction
from gelato.fem.core import VectorTestFunction
from gelato.fem.expr import BilinearForm
from gelato.fem.expr import atomize, normalize, matricize
from gelato.fem.expr import gelatize


# ...
def test_atomize_2d_1():
    print('============ test_atomize_2d_1 =============')

    V = BasicSobolevSpace('V', ldim=2)

    v = TestFunction(V, name='v')
    w = TestFunction(V, name='w')
    c = Constant('c')
    F = Field('F')

    # ... expressions that can be normalized (valid for a weak formulation)
    assert(atomize(grad(v)) == Tuple(dx(v),
                                      dy(v)))
    assert(atomize(grad(c*v)) == Tuple(c*dx(v),
                                        c*dy(v)))
    assert(atomize(grad(F*v)) == Tuple(F*dx(v) + v*dx(F),
                                        F*dy(v) + v*dy(F)))

    assert(atomize(dot(grad(v), grad(w))) == dx(v)*dx(w) + dy(v)*dy(w))
    # ...

#    expr = grad(F*v)
#    print('> input         >>> {0}'.format(expr))
#    print('> atomized     >>> {0}'.format(atomize(expr)))
# ...

# ...
def test_normalize_2d_1():
    print('============ test_normalize_2d_1 =============')

    V = BasicSobolevSpace('V', ldim=2)
    U = BasicSobolevSpace('U', ldim=2)

    v = TestFunction(V, name='v')
    u = TestFunction(U, name='u')
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
    assert(normalize(dot(grad(v), grad(u)), basis={V: 'Ni', U: 'Nj'}) == Ni_x*Nj_x + Ni_y*Nj_y)
    assert(normalize(dot(grad(v), grad(u)) + c*v*u, basis={V: 'Ni', U: 'Nj'}) == Ni_x*Nj_x + Ni_y*Nj_y + c*Ni*Nj)
    # ...

#    expr = dot(grad(v), b)
#    print('> input         >>> {0}'.format(expr))
#
#    print('> normal form   >>> {0}'.format(normalize(expr, basis={V: 'Ni'})))
#    print('> normal form   >>> {0}'.format(normalize(dot(grad(v), grad(u)),
#                                                     basis={V: 'Ni', U: 'Nj'})))
# ...


# ...
def test_atomize_2d_2():
    print('============ test_atomize_2d_2 =============')

    V = BasicSobolevSpace('V', ldim=2, is_block=True, shape=2)

    v = VectorTestFunction(V, name='v')

    assert(atomize(rot(v)) == -dx(v[1]) + dy(v[0]))
    assert(atomize(div(v)) == dx(v[0]) + dy(v[1]))

#    expr = div(v)
#    print('> input         >>> {0}'.format(expr))
#    print('> atomized     >>> {0}'.format(atomize(expr)))
# ...

# ...
def test_normalize_2d_2():
    print('============ test_normalize_2d_2 =============')

    V = BasicSobolevSpace('V', ldim=2, is_block=True, shape=2)
    U = BasicSobolevSpace('U', ldim=2, is_block=True, shape=2)

    v = VectorTestFunction(V, name='v')
    u = VectorTestFunction(U, name='u')

    Ni = IndexedBase('Ni', shape=2)
    Ni_x = IndexedBase('Ni_x', shape=2)
    Ni_y = IndexedBase('Ni_y', shape=2)

    Nj = IndexedBase('Nj', shape=2)
    Nj_x = IndexedBase('Nj_x', shape=2)
    Nj_y = IndexedBase('Nj_y', shape=2)

    assert(normalize(v[0], basis={V: 'Ni'}) == Ni[0])
    assert(normalize(dx(v[0]), basis={V: 'Ni'}) == Ni_x[0])
    assert(normalize(div(v), basis={V: 'Ni'}) == Ni_x[0] + Ni_y[1])
    assert(normalize(rot(v), basis={V: 'Ni'}) == -Ni_x[1] + Ni_y[0])

    assert(normalize(v[0]*u[0], basis={V: 'Ni', U: 'Nj'}) == Ni[0]*Nj[0])
    assert(normalize(v[1]*dx(u[0]), basis={V: 'Ni', U: 'Nj'}) == Ni[1]*Nj_x[0])
    assert(normalize(dy(v[0])*u[1], basis={V: 'Ni', U: 'Nj'}) == Ni_y[0]*Nj[1])
    assert(normalize(dx(v[1])*dy(u[1]), basis={V: 'Ni', U: 'Nj'}) == Ni_x[1]*Nj_y[1])

    expected = (Ni_x[0] + Ni_y[1]) * (Nj_x[0] + Nj_y[1])
    assert(normalize(div(v) * div(u), basis={V: 'Ni', U: 'Nj'}) == expected)

    expected = (-Ni_x[1] + Ni_y[0]) * (-Nj_x[1] + Nj_y[0])
    assert(normalize(rot(v) * rot(u), basis={V: 'Ni', U: 'Nj'}) == expected)

#    expr = dx(v[0])
#    print('> input         >>> {0}'.format(expr))

#    print('> normal form   >>> {0}'.format(normalize(expr, basis={V: 'Ni'})))
##    print('> normal form   >>> {0}'.format(normalize(expr, basis={V: 'Ni', U: 'Nj'})))
# ...

# ...
def test_matricize_2d_2():
    print('============ test_matricize_2d_2 =============')

    V = BasicSobolevSpace('V', ldim=2, is_block=True, shape=2)
    U = BasicSobolevSpace('U', ldim=2, is_block=True, shape=2)

    v = VectorTestFunction(V, name='v')
    u = VectorTestFunction(U, name='u')

    Ni, Ni_x, Ni_y = symbols('Ni Ni_x Ni_y')
    Nj, Nj_x, Nj_y = symbols('Nj Nj_x Nj_y')

    c1 = Symbol('c1')
    c2 = Symbol('c2')

    # ...
    expr = v[0]*u[0]
    expr = normalize(expr, basis={V: 'Ni', U: 'Nj'})
    expected = Matrix([[Ni*Nj, 0], [0, 0]])
    assert(matricize(expr) == expected)
    # ...

    # ...
    expr = v[1]*dx(u[0])
    expr = normalize(expr, basis={V: 'Ni', U: 'Nj'})
    expected = Matrix([[0, Ni*Nj_x], [0, 0]])
    assert(matricize(expr) == expected)
    # ...

    # ...
    expr = dy(v[0])*u[1]
    expr = normalize(expr, basis={V: 'Ni', U: 'Nj'})
    expected = Matrix([[0, 0], [Ni_y*Nj, 0]])
    assert(matricize(expr) == expected)
    # ...

    # ...
    expr = dx(v[1])*dy(u[1])
    expr = normalize(expr, basis={V: 'Ni', U: 'Nj'})
    expected = Matrix([[0, 0], [0, Ni_x*Nj_y]])
    assert(matricize(expr) == expected)
    # ...

    # ...
    expr = div(v) * div(u)
    expr = normalize(expr, basis={V: 'Ni', U: 'Nj'})
    expected = Matrix([[Ni_x*Nj_x, Ni_y*Nj_x], [Ni_x*Nj_y, Ni_y*Nj_y]])
    assert(matricize(expr) == expected)
    # ...

    # ...
    expr = rot(v) * rot(u)
    expr = normalize(expr, basis={V: 'Ni', U: 'Nj'})
    expected = Matrix([[Ni_y*Nj_y, -Ni_x*Nj_y], [-Ni_y*Nj_x, Ni_x*Nj_x]])
    assert(matricize(expr) == expected)
    # ...

    # ...
    expr = div(v) * div(u) + rot(v) * rot(u)
    expr = normalize(expr, basis={V: 'Ni', U: 'Nj'})
    expected = Matrix([[Ni_x*Nj_x + Ni_y*Nj_y, -Ni_x*Nj_y + Ni_y*Nj_x],
                       [Ni_x*Nj_y - Ni_y*Nj_x, Ni_x*Nj_x + Ni_y*Nj_y]])
    assert(matricize(expr) == expected)
    # ...

    # ...
    expr = c1 * div(v) * div(u) + rot(v) * rot(u)
    expr = normalize(expr, basis={V: 'Ni', U: 'Nj'})
    expected = Matrix([[Ni_x*Nj_x*c1 + Ni_y*Nj_y, -Ni_x*Nj_y + Ni_y*Nj_x*c1],
                       [Ni_x*Nj_y*c1 - Ni_y*Nj_x, Ni_x*Nj_x + Ni_y*Nj_y*c1]])
    assert(matricize(expr) == expected)
    # ...

    # ...
    expr = c1 * div(v) * div(u) + c2 * rot(v) * rot(u)
    expr = normalize(expr, basis={V: 'Ni', U: 'Nj'})
    expected = Matrix([[Ni_x*Nj_x*c1 + Ni_y*Nj_y*c2, -Ni_x*Nj_y*c2 + Ni_y*Nj_x*c1],
                       [Ni_x*Nj_y*c1 - Ni_y*Nj_x*c2, Ni_x*Nj_x*c2 + Ni_y*Nj_y*c1]])
    assert(matricize(expr) == expected)
    # ...

#    expr = c1 * div(v) * div(u) + rot(v) * rot(u)
#    print('> input         >>> {0}'.format(expr))
#    expr = normalize(expr, basis={V: 'Ni', U: 'Nj'})
#    print('> matricize     >>> {0}'.format(matricize(expr)))
# ...

# ...
def test_bilinear_form_2d_1():
    print('============ test_bilinear_form_2d_1 =============')

    U = BasicSobolevSpace('U', ldim=2)
    V = BasicSobolevSpace('V', ldim=2)

    u = TestFunction(U, name='u')
    v = TestFunction(V, name='v')

    expr = inner(grad(u), grad(v)) + u*v

    a = BilinearForm((v,u), expr)
    print('> input      >>> {0}'.format(a))

    a_expr = gelatize(a, basis={V: 'Nj', U: 'Ni'})
    print('> gelatized  >>> {0}'.format(a_expr))
    print('')
# ...

# ...
def test_bilinear_form_2d_2():
    print('============ test_bilinear_form_2d_2 =============')

    U = BasicSobolevSpace('U', ldim=2)
    V = BasicSobolevSpace('V', ldim=2)

    u = TestFunction(U, name='u')
    v = TestFunction(V, name='v')

    expr = cross(curl(u), curl(v)) + 0.2 * u * v

    a = BilinearForm((v,u), expr)
    print('> input      >>> {0}'.format(a))

    a_expr = gelatize(a, basis={V: 'Nj', U: 'Ni'})
    print('> gelatized  >>> {0}'.format(a_expr))
    print('')
# ...

# ...
def test_bilinear_form_2d_3():
    print('============ test_bilinear_form_2d_3 =============')

    U = BasicSobolevSpace('U', ldim=2)
    V = BasicSobolevSpace('V', ldim=2)

    u = TestFunction(U, name='u')
    v = TestFunction(V, name='v')

    bx = Constant('bx')
    by = Constant('by')
    b = Tuple(bx, by)

    expr = 0.2 * u * v + dot(b, grad(v)) * u

    a = BilinearForm((v,u), expr)
    print('> input      >>> {0}'.format(a))

    a_expr = gelatize(a, basis={V: 'Nj', U: 'Ni'})
    print('> gelatized  >>> {0}'.format(a_expr))
    print('')
# ...

# ...
def test_bilinear_form_2d_4():
    print('============ test_bilinear_form_2d_4 =============')

    U = BasicSobolevSpace('U', ldim=2, is_block=True, shape=2)
    V = BasicSobolevSpace('V', ldim=2, is_block=True, shape=2)

    u = VectorTestFunction(U, name='u')
    v = VectorTestFunction(V, name='v')

    expr = rot(u) * rot(v) + div(u) * div(v) + 0.2 * dot(u, v)

    a = BilinearForm((v,u), expr)
    print('> input      >>> {0}'.format(a))

    a_expr = gelatize(a, basis={V: 'Nj', U: 'Ni'})
    print('> gelatized  >>> {0}'.format(a_expr))
    print('')
# ...

# ...
def test_bilinear_form_2d_10():
    print('============ test_bilinear_form_2d_10 =============')

    U = BasicSobolevSpace('U', ldim=2)
    V = BasicSobolevSpace('V', ldim=2)

    u = TestFunction(U, name='u')
    v = TestFunction(V, name='v')

    u1 = TestFunction(U, name='u1')
    v1 = TestFunction(V, name='v1')

    Ni, Ni_x, Ni_y = symbols('Ni Ni_x Ni_y')
    Nj, Nj_x, Nj_y = symbols('Nj Nj_x Nj_y')

    c1 = Symbol('c1')
    c2 = Symbol('c2')

    a = BilinearForm((v,u), inner(grad(u), grad(v)))
    b = BilinearForm((v,u), u*v)
    adv = BilinearForm((v,u), dx(u)*v)

    # ...
    expected = Ni*Nj + Ni_x*Nj_x + Ni_y*Nj_y
    assert(gelatize(a + b, basis={V: 'Nj', U: 'Ni'}) == expected)
    # ...

    # ...
    expected = 2*Ni_x*Nj_x + 2*Ni_y*Nj_y
    assert(gelatize(2 * a, basis={V: 'Nj', U: 'Ni'}) == expected)
    # ...

    # ...
    expected = c1*(Ni_x*Nj_x + Ni_y*Nj_y)
    assert(gelatize(c1*a, basis={V: 'Nj', U: 'Ni'}) == expected)
    # ...

    # ...
    expected = Ni*Nj*c2 + c1*(Ni_x*Nj_x + Ni_y*Nj_y)
    assert(gelatize(c1*a + c2*b, basis={V: 'Nj', U: 'Ni'}) == expected)
    # ...

    # ...
    expected = c1*(Ni_x*Nj_x + Ni_y*Nj_y) + c2*(Ni*Nj + Ni_x*Nj)
    assert(gelatize(c1*a  + c2*(b + adv), basis={V: 'Nj', U: 'Ni'}) == expected)
    # ...

    # ...
    assert(gelatize(a(u1, v1), basis={V: 'Nj', U: 'Ni'}) == gelatize(a(v1, u1), basis={V: 'Nj', U: 'Ni'}))
    # ...

    # ...
    expected = Ni_x*Nj
    assert(gelatize(adv(v1, u1), basis={V: 'Nj', U: 'Ni'}) == expected)

    expected = Nj_x*Ni
    assert(gelatize(adv(u1, v1), basis={V: 'Nj', U: 'Ni'}) == expected)
    # ...

#    expr = c1*a  + c2*(b + adv)
#    print('> input      >>> {0}'.format(expr))
#    print('> gelatized  >>> {0}'.format(gelatize(expr, basis={V: 'Nj', U: 'Ni'}) ))
#    print('')
# ...


# .....................................................
if __name__ == '__main__':
    test_atomize_2d_1()
    test_normalize_2d_1()

    test_atomize_2d_2()
    test_normalize_2d_2()
    test_matricize_2d_2()

    test_bilinear_form_2d_1()
    test_bilinear_form_2d_2()
    test_bilinear_form_2d_3()
    test_bilinear_form_2d_4()
    test_bilinear_form_2d_10()
