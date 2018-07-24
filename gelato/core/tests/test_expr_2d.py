# coding: utf-8

# TODO - split the asserts between algebraic and weak formulations ones
#      - add assert for grad in vector case
# TODO: - __call__ examples are not working anymore

from sympy import Symbol
from sympy.core.containers import Tuple
from sympy import symbols
from sympy import IndexedBase
from sympy import Matrix
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
from gelato.core import atomize, normalize, matricize
from gelato.core import gelatize
from gelato.core import tensorize
from gelato.core import Mass, Stiffness, Advection, AdvectionT


# ...
def test_atomize_2d_1():
    print('============ test_atomize_2d_1 =============')

    V = H1Space('V', ldim=2)

    v = TestFunction(V, name='v')
    w = TestFunction(V, name='w')
    c = Constant('c')
    F = Field('F', space=V)

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

    V = H1Space('V', ldim=2)
    U = H1Space('U', ldim=2)

    v = TestFunction(V, name='v')
    u = TestFunction(U, name='u')

    x,y = V.coordinates

    c = Constant('c')
    F = Field('F', space=V)
    f1 = Function('f1')
    f2 = Function('f2')

    Ni, Ni_x, Ni_y = symbols('Ni Ni_x Ni_y')
    Nj, Nj_x, Nj_y = symbols('Nj Nj_x Nj_y')

    bx, by = symbols('bx by')
    b = Tuple(bx, by)

    f = Tuple(f1(x,y), f2(x,y))

    a00 = Constant('a00')
    a10 = Constant('a10')
    a01 = Constant('a01')
    a11 = Constant('a11')
    A = Matrix([[a00, a01], [a10, a11]])

    # ...
    assert(normalize(grad(v), basis={V: 'Ni'}) == Tuple(Ni_x, Ni_y))
    assert(normalize(grad(c*v), basis={V: 'Ni'}) == Tuple(c*Ni_x, c*Ni_y))
    assert(normalize(dot(b, grad(v)), basis={V: 'Ni'}) == Ni_x*bx + Ni_y*by)
    assert(normalize(dot(b, grad(v)) + c*v, basis={V: 'Ni'}) == Ni_x*bx + Ni_y*by + c*Ni)
    assert(normalize(dot(f, grad(v)), basis={V: 'Ni'}) == Ni_x*f1(x,y) + Ni_y*f2(x,y))
    assert(normalize(dot(Tuple(2, 3), grad(v)), basis={V: 'Ni'}) == 2*Ni_x + 3*Ni_y)
    assert(normalize(grad(F*v), basis={V: 'Ni'}) == Tuple(F*Ni_x + Ni*dx(F),
                                                          F*Ni_y + Ni*dy(F)))
    # TODO debug
#    assert(normalize(A*grad(v), basis={V: 'Ni'}) == 2*Ni_x + 3*Ni_y)

    assert(normalize(dot(grad(v), grad(u)), basis={V: 'Ni', U: 'Nj'}) == Ni_x*Nj_x + Ni_y*Nj_y)
    assert(normalize(dot(grad(v), grad(u)) + c*v*u, basis={V: 'Ni', U: 'Nj'}) == Ni_x*Nj_x + Ni_y*Nj_y + c*Ni*Nj)
    assert(normalize(dot(grad(F*v), grad(u)), basis={V: 'Ni', U: 'Nj'}) == Nj_x*(F*Ni_x + Ni*dx(F)) + Nj_y*(F*Ni_y + Ni*dy(F)))
    # ...

#    expr = dot(A, grad(v))
#    expr = div(dot(A, grad(v)))
#    print('> input         >>> {0}'.format(expr))

#    print('> normal form   >>> {0}'.format(normalize(expr, basis={V: 'Ni'})))
#    print('> normal form   >>> {0}'.format(normalize(expr, basis={V: 'Ni', U: 'Nj'})))
# ...


# ...
def test_gelatize_2d_1():
    print('============ test_gelatize_2d_1 =============')

    V = H1Space('V', ldim=2)
    U = H1Space('U', ldim=2)

    v = TestFunction(V, name='v')
    u = TestFunction(U, name='u')

    x,y = V.coordinates

    c = Constant('c')
    F = Field('F', space=V)
    f1 = Function('f1')
    f2 = Function('f2')

    Ni, Ni_x, Ni_y = symbols('Ni Ni_x Ni_y')
    Nj, Nj_x, Nj_y = symbols('Nj Nj_x Nj_y')

    bx, by = symbols('bx by')
    b = Tuple(bx, by)

    f = Tuple(f1(x,y), f2(x,y))

    a00 = Constant('a00')
    a10 = Constant('a10')
    a01 = Constant('a01')
    a11 = Constant('a11')
    A = Matrix([[a00, a01], [a10, a11]])

    # ...
    a = BilinearForm((v,u), dot(grad(v), grad(u)))
    assert(gelatize(a, basis={V: 'Ni', U: 'Nj'}) == Ni_x*Nj_x + Ni_y*Nj_y)
    # ...

    # ...
    a = BilinearForm((v,u), dot(grad(v), grad(u)) + c*v*u)
    assert(gelatize(a, basis={V: 'Ni', U: 'Nj'}) == Ni_x*Nj_x + Ni_y*Nj_y + c*Ni*Nj)
    # ...

    # ...
    a = BilinearForm((v,u), dot(grad(v), grad(u)) + F*v*u)
    assert(gelatize(a, basis={V: 'Ni', U: 'Nj'}) == Ni_x*Nj_x + Ni_y*Nj_y + F*Ni*Nj)
    # ...

    # ...
    a = BilinearForm((v,u), dot(grad(F*v), grad(u)))
    assert(gelatize(a, basis={V: 'Ni', U: 'Nj'}) == F*Ni_x*Nj_x + F*Ni_y*Nj_y + Ni*Nj_x*dx(F) + Ni*Nj_y*dy(F))
    # ...

# ...

# ...
def test_atomize_2d_2():
    print('============ test_atomize_2d_2 =============')

    V = H1Space('V', ldim=2, is_block=True, shape=2)

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

    V = H1Space('V', ldim=2, is_block=True, shape=2)
    U = H1Space('U', ldim=2, is_block=True, shape=2)

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

    V = H1Space('V', ldim=2, is_block=True, shape=2)
    U = H1Space('U', ldim=2, is_block=True, shape=2)

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

    U = H1Space('U', ldim=2)
    V = H1Space('V', ldim=2)

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

    U = H1Space('U', ldim=2)
    V = H1Space('V', ldim=2)

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

    U = H1Space('U', ldim=2)
    V = H1Space('V', ldim=2)

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

    U = H1Space('U', ldim=2, is_block=True, shape=2)
    V = H1Space('V', ldim=2, is_block=True, shape=2)

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

    U = H1Space('U', ldim=2)
    V = H1Space('V', ldim=2)

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

#    # ... TODO debug
#    expected = Ni_x*Nj
#    assert(gelatize(adv(v1, u1), basis={V: 'Nj', U: 'Ni'}) == expected)
#
#    expected = Nj_x*Ni
#    assert(gelatize(adv(u1, v1), basis={V: 'Nj', U: 'Ni'}) == expected)
#    # ...

#    expr = c1*a  + c2*(b + adv)
#    print('> input      >>> {0}'.format(expr))
#    print('> gelatized  >>> {0}'.format(gelatize(expr, basis={V: 'Nj', U: 'Ni'}) ))
#    print('')
# ...

# ...
def test_linear_form_2d_10():
    print('============ test_linear_form_2d_10 =============')

    V = H1Space('V', ldim=2)

    v = TestFunction(V, name='v')

    x,y = V.coordinates
    f = Function('f')
    g = Function('g')

    Ni, Ni_x, Ni_y = symbols('Ni Ni_x Ni_y')
    Nj, Nj_x, Nj_y = symbols('Nj Nj_x Nj_y')

    c1 = Symbol('c1')
    c2 = Symbol('c2')

    bx, by = symbols('bx by')
    b = Tuple(bx, by)
    fg = Tuple(f(x,y), g(x,y))

    a = LinearForm(v, cos(2*pi*x)*cos(4*pi*y)*v)

    # ...
    expected = cos(2*pi*x)*cos(4*pi*y)*Ni
    assert(gelatize(LinearForm(v, cos(2*pi*x)*cos(4*pi*y)*v),
                    basis={V: 'Ni'}) == expected)
    # ...

    # ...
    expected = f(x,y)*Ni
    assert(gelatize(LinearForm(v, f(x,y)*v),
                    basis={V: 'Ni'}) == expected)
    # ...

    # ...
    expected = bx*Ni_x + by*Ni_y + f(x,y)*Ni
    assert(gelatize(LinearForm(v, dot(b, grad(v)) + f(x,y)*v),
                    basis={V: 'Ni'}) == expected)
    # ...

    # ...
    expected = f(x,y)*Ni_x + g(x,y)*Ni_y
    assert(gelatize(LinearForm(v, dot(fg, grad(v))),
                    basis={V: 'Ni'}) == expected)
    # ...

#    expr =
#    print('> input      >>> {0}'.format(expr))
#    print('> gelatized  >>> {0}'.format(gelatize(expr, basis={V: 'Ni'}) ))
#    print('')
# ...

# ...
def test_function_form_2d_10():
    print('============ test_function_form_2d_10 =============')

    V = H1Space('V', ldim=2)

    F = Field('F', space=V)

    x,y = V.coordinates

    f = Function('f')
    g = Function('g')

    Ni, Ni_x, Ni_y = symbols('Ni Ni_x Ni_y')
    Nj, Nj_x, Nj_y = symbols('Nj Nj_x Nj_y')

    c1 = Symbol('c1')
    c2 = Symbol('c2')

    bx, by = symbols('bx by')
    b = Tuple(bx, by)
    fg = Tuple(f(x,y), g(x,y))

    # ...
    expected = 4*pi**2*sin(2*pi*x)**2*cos(3*pi*y)**2 + 9*pi**2*sin(3*pi*y)**2*cos(2*pi*x)**2
    e = cos(2*pi*x)*cos(3*pi*y)
    assert(gelatize(FunctionForm(dot(grad(e), grad(e)), coordinates=[x,y])) == expected)
    # ...

    # ...
    expected = F - cos(2*pi*x)*cos(3*pi*y)
    assert(gelatize(FunctionForm(F-cos(2*pi*x)*cos(3*pi*y))) == expected)
    # ...

    # ...
    expected = (F - cos(2*pi*x)*cos(3*pi*y))**2
    assert(gelatize(FunctionForm((F - cos(2*pi*x)*cos(3*pi*y))**2)) == expected)
    # ...

    # ...
    expected = (dx(F) + 2*pi*sin(2*pi*x)*cos(3*pi*y))**2 + (dy(F) + 3*pi*sin(3*pi*y)*cos(2*pi*x))**2
    e = F -cos(2*pi*x)*cos(3*pi*y)
    assert(gelatize(FunctionForm(dot(grad(e), grad(e)))) == expected)
    # ...

#    e = F -cos(2*pi*x)*cos(3*pi*y)
#    expr = FunctionForm(dot(grad(e), grad(e)))
#    print('> input      >>> {0}'.format(expr))
#    print('> gelatized  >>> {0}'.format(gelatize(expr) ))
#    print('')
# ...

# ...
def test_linear_form_2d_1():
    print('============ test_linear_form_2d_1 =============')

    V = H1Space('V', ldim=2)

    v = TestFunction(V, name='v')

    x,y = V.coordinates
    f = Function('f')

    expr = cos(2*pi*x)*cos(4*pi*y)*v
#    expr = f(x)*v

    a = LinearForm(v, expr)
    print('> input      >>> {0}'.format(a))

    a_expr = gelatize(a, basis={V: 'Nj'})
    print('> gelatized  >>> {0}'.format(a_expr))
    print('')
# ...

# ...
def test():
    print('============ test =============')

    V = H1Space('V', ldim=2, is_block=True, shape=2)
    U = H1Space('U', ldim=2)

    v = VectorTestFunction(V, name='v')
    u = TestFunction(U, name='u')

    Ni, Ni_x, Ni_y = symbols('Ni Ni_x Ni_y')
    Nj, Nj_x, Nj_y = symbols('Nj Nj_x Nj_y')

    # ...
    expr = v[0]*u
    expr = normalize(expr, basis={V: 'Ni', U: 'Nj'})
    expected = Matrix([[Ni*Nj], [0]])
    assert(matricize(expr) == expected)
    # ...

    # ...
    expr = dot(v, grad(u))
    expr = normalize(expr, basis={V: 'Ni', U: 'Nj'})
    expected = Matrix([[Ni*Nj_x], [Ni*Nj_y]])
    assert(matricize(expr) == expected)
    # ...

#    expr = v[0]*u
#    print('> input         >>> {0}'.format(expr))
#    expr = normalize(expr, basis={V: 'Ni', U: 'Nj'})
#    print('> matricize     >>> {0}'.format(matricize(expr)))

# ...
def test_tensorize_2d_1():
    print('============ test_tensorize_2d_1 =============')

    V = H1Space('V', ldim=2)
    V_0 = H1Space('V_0', ldim=1, coordinates=['x'])
    V_1 = H1Space('V_1', ldim=1, coordinates=['y'])

    v = TestFunction(V, name='v')
    u = TestFunction(V, name='u')

    v0 = TestFunction(V_0, name='v0')
    u0 = TestFunction(V_0, name='u0')

    v1 = TestFunction(V_1, name='v1')
    u1 = TestFunction(V_1, name='u1')

    c = Constant('c')

    # ...
    expected = Mass(v1, u1)*Mass(v0, u0)
    assert(tensorize(BilinearForm((v,u), u*v)) == expected)
    # ...

    # ...
    expected = Mass(v1, u1)*Stiffness(v0, u0)
    assert(tensorize(BilinearForm((v,u), dx(u)*dx(v))) == expected)
    # ...

    # ...
    expected = Advection(v1, u1)*Mass(v0, u0)
    assert(tensorize(BilinearForm((v,u), dy(u) * v)) == expected)
    # ...

    # ...
    expected =  Mass(v1,u1)*Advection(v0,u0)
    assert(tensorize(BilinearForm((v,u), dx(u) * v)) == expected)
    # ...

    # ...
    expected = Mass(v1,u1)*Stiffness(v0,u0) + Stiffness(v1,u1)*Mass(v0,u0)
    assert(tensorize(BilinearForm((v,u), dot(grad(v), grad(u)))) == expected)
    # ...

    # ...
    expected = Advection(v1,u1)*Mass(v0,u0) + Mass(v1,u1)*Advection(v0,u0) + Mass(v1,u1)*Stiffness(v0,u0) + Stiffness(v1,u1)*Mass(v0,u0)
    assert(tensorize(BilinearForm((v,u), dot(grad(v), grad(u)) + dx(u)*v + dy(u)*v)) == expected)
    # ...

#    expr = dot(grad(v), grad(u)) + dx(u)*v + dy(u)*v
#    expr = BilinearForm((v,u), expr)
#
#    print('> input         >>> {0}'.format(expr))
#    print('> tensorized    >>> {0}'.format(tensorize(expr)))
# ...


# .....................................................
if __name__ == '__main__':
#    test_atomize_2d_1()
#    test_normalize_2d_1()
#    test_gelatize_2d_1()
#
#    test_atomize_2d_2()
#    test_normalize_2d_2()
#    test_matricize_2d_2()
#
###    test_bilinear_form_2d_1()
###    test_bilinear_form_2d_2()
###    test_bilinear_form_2d_3()
###    test_bilinear_form_2d_4()
#
#    test_bilinear_form_2d_10()
#    test_linear_form_2d_10()
#    test_function_form_2d_10()
#
#    test_linear_form_2d_1()
#
#    test()

    test_tensorize_2d_1()
