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
from sympy import pi, cos

from gelato.core import dx, dy, dz
from gelato.core import Constant
from gelato.core import Field
from gelato.core import grad, dot, inner, cross, rot, curl, div
from gelato.core import H1Space
from gelato.core import TestFunction
from gelato.core import VectorTestFunction
from gelato.core import BilinearForm, LinearForm
from gelato.core import atomize, normalize, matricize
from gelato.core import gelatize


# ...
def test_atomize_3d_1():
    print('============ test_atomize_3d_1 =============')

    V = H1Space('V', ldim=3)

    v = TestFunction(V, name='v')
    w = TestFunction(V, name='w')
    c = Constant('c')
    F = Field('F', space=V)

    # ... expressions that can be normalized (valid for a weak formulation)
    assert(atomize(grad(v)) == Tuple(dx(v),
                                      dy(v),
                                      dz(v)))
    assert(atomize(grad(c*v)) == Tuple(c*dx(v),
                                        c*dy(v),
                                        c*dz(v)))
    assert(atomize(grad(F*v)) == Tuple(F*dx(v) + v*dx(F),
                                        F*dy(v) + v*dy(F),
                                        F*dz(v) + v*dz(F)))

    assert(atomize(dot(grad(v), grad(w))) == dx(v)*dx(w) + dy(v)*dy(w) + dz(v)*dz(w))
    # ...

#    expr = grad(F*v)
#    print('> input         >>> {0}'.format(expr))
#    print('> atomized     >>> {0}'.format(atomize(expr)))
# ...

# ...
def test_normalize_3d_1():
    print('============ test_normalize_3d_1 =============')

    V = H1Space('V', ldim=3)
    U = H1Space('U', ldim=3)

    v = TestFunction(V, name='v')
    u = TestFunction(U, name='u')

    x,y,z = V.coordinates

    c = Constant('c')
    F = Field('F', space=V)
    f1 = Function('f1')
    f2 = Function('f2')
    f3 = Function('f3')

    Ni, Ni_x, Ni_y, Ni_z = symbols('Ni Ni_x Ni_y Ni_z')
    Nj, Nj_x, Nj_y, Nj_z = symbols('Nj Nj_x Nj_y Nj_z')

    bx, by, bz = symbols('bx by bz')
    b = Tuple(bx, by, bz)

    f = Tuple(f1(x,y,z), f2(x,y,z), f3(x,y,z))

    a00 = Constant('a00')
    a10 = Constant('a10')
    a20 = Constant('a20')
    a01 = Constant('a01')
    a11 = Constant('a11')
    a21 = Constant('a21')
    a02 = Constant('a02')
    a12 = Constant('a12')
    a22 = Constant('a22')
    A = Matrix([[a00, a01, a02], [a10, a11, a12], [a20, a21, a22]])

    # ...
    assert(normalize(grad(v), basis={V: 'Ni'}) == Tuple(Ni_x, Ni_y, Ni_z))
    assert(normalize(grad(c*v), basis={V: 'Ni'}) == Tuple(c*Ni_x, c*Ni_y, c*Ni_z))
    assert(normalize(dot(b, grad(v)), basis={V: 'Ni'}) == Ni_x*bx + Ni_y*by + Ni_z*bz)
    assert(normalize(dot(b, grad(v)) + c*v, basis={V: 'Ni'}) == Ni_x*bx + Ni_y*by + Ni_z*bz + c*Ni)
    assert(normalize(dot(f, grad(v)), basis={V: 'Ni'}) == Ni_x*f1(x,y,z) + Ni_y*f2(x,y,z) + Ni_z*f3(x,y,z))
    assert(normalize(dot(Tuple(2, 3, 4), grad(v)), basis={V: 'Ni'}) == 2*Ni_x + 3*Ni_y + 4*Ni_z)
    assert(normalize(grad(F*v), basis={V: 'Ni'}) == Tuple(F*Ni_x + Ni*dx(F),
                                                          F*Ni_y + Ni*dy(F),
                                                          F*Ni_z + Ni*dz(F)))

    assert(normalize(dot(grad(v), grad(u)), basis={V: 'Ni', U: 'Nj'}) == Ni_x*Nj_x + Ni_y*Nj_y + Ni_z*Nj_z)
    assert(normalize(dot(grad(v), grad(u)) + c*v*u, basis={V: 'Ni', U: 'Nj'}) == Ni_x*Nj_x + Ni_y*Nj_y + Ni_z*Nj_z + c*Ni*Nj)
    # ...

#    expr = dot(A, grad(v))
#    expr = div(dot(A, grad(v)))
#    print('> input         >>> {0}'.format(expr))
#
#    print('> normal form   >>> {0}'.format(normalize(expr, basis={V: 'Ni'})))
#    print('> normal form   >>> {0}'.format(normalize(dot(grad(v), grad(u)),
#                                                     basis={V: 'Ni', U: 'Nj'})))
# ...

# ...
def test_atomize_3d_2():
    print('============ test_atomize_3d_2 =============')

    V = H1Space('V', ldim=3, is_vector=True, shape=3)

    v = VectorTestFunction(V, name='v')

    assert(atomize(curl(v)) == Tuple( dy(v[2]) - dz(v[1]),
                                      -dx(v[2]) + dz(v[0]),
                                       dx(v[1]) - dy(v[0])))
    assert(atomize(div(v)) == dx(v[0]) + dy(v[1]) + dz(v[2]))

#    expr = curl(v)
#    print('> input         >>> {0}'.format(expr))
#    print('> atomized     >>> {0}'.format(atomize(expr)))
# ...

# ...
def test_normalize_3d_2():
    print('============ test_normalize_3d_2 =============')

    V = H1Space('V', ldim=3, is_block=True, shape=3)
    U = H1Space('U', ldim=3, is_block=True, shape=3)

    v = VectorTestFunction(V, name='v')
    u = VectorTestFunction(U, name='u')

    Ni = IndexedBase('Ni', shape=3)
    Ni_x = IndexedBase('Ni_x', shape=3)
    Ni_y = IndexedBase('Ni_y', shape=3)
    Ni_z = IndexedBase('Ni_z', shape=3)

    Nj = IndexedBase('Nj', shape=3)
    Nj_x = IndexedBase('Nj_x', shape=3)
    Nj_y = IndexedBase('Nj_y', shape=3)
    Nj_z = IndexedBase('Nj_z', shape=3)

    assert(normalize(v[0], basis={V: 'Ni'}) == Ni[0])
    assert(normalize(dx(v[0]), basis={V: 'Ni'}) == Ni_x[0])
    assert(normalize(div(v), basis={V: 'Ni'}) == Ni_x[0] + Ni_y[1] + Ni_z[2])
    assert(normalize(curl(v), basis={V: 'Ni'}) == Tuple( Ni_y[2] - Ni_z[1],
                                                        -Ni_x[2] + Ni_z[0],
                                                         Ni_x[1] - Ni_y[0]))

    assert(normalize(v[0]*u[0], basis={V: 'Ni', U: 'Nj'}) == Ni[0]*Nj[0])
    assert(normalize(v[1]*dx(u[0]), basis={V: 'Ni', U: 'Nj'}) == Ni[1]*Nj_x[0])
    assert(normalize(dy(v[0])*u[1], basis={V: 'Ni', U: 'Nj'}) == Ni_y[0]*Nj[1])
    assert(normalize(dx(v[1])*dy(u[1]), basis={V: 'Ni', U: 'Nj'}) == Ni_x[1]*Nj_y[1])
    assert(normalize(dx(v[1])*dz(u[2]), basis={V: 'Ni', U: 'Nj'}) == Ni_x[1]*Nj_z[2])

    expected = (Ni_x[0] + Ni_y[1] + Ni_z[2])*(Nj_x[0] + Nj_y[1] + Nj_z[2])
    assert(normalize(div(v) * div(u), basis={V: 'Ni', U: 'Nj'}) == expected)

    expected = ((Ni_x[1] - Ni_y[0])*(Nj_x[1] - Nj_y[0])
             + (-Ni_x[2] + Ni_z[0])*(-Nj_x[2] + Nj_z[0])
             + (Ni_y[2] - Ni_z[1])*(Nj_y[2] - Nj_z[1]))
    assert(normalize(dot(curl(v), curl(u)), basis={V: 'Ni', U: 'Nj'}) == expected)

#    expr = dot(curl(v), curl(u))
#    print('> input         >>> {0}'.format(expr))
#
##    print('> normal form   >>> {0}'.format(normalize(expr, basis={V: 'Ni'})))
#    print('> normal form   >>> {0}'.format(normalize(expr, basis={V: 'Ni', U: 'Nj'})))
# ...

# ...
def test_matricize_3d_2():
    print('============ test_matricize_3d_2 =============')

    V = H1Space('V', ldim=3, is_block=True, shape=3)
    U = H1Space('U', ldim=3, is_block=True, shape=3)

    v = VectorTestFunction(V, name='v')
    u = VectorTestFunction(U, name='u')

    Ni, Ni_x, Ni_y, Ni_z = symbols('Ni Ni_x Ni_y Ni_z')
    Nj, Nj_x, Nj_y, Nj_z = symbols('Nj Nj_x Nj_y Nj_z')

    c1 = Symbol('c1')
    c2 = Symbol('c2')

    # ...
    expr = v[0]*u[0]
    expr = normalize(expr, basis={V: 'Ni', U: 'Nj'})
    expected = Matrix([[Ni*Nj, 0, 0], [0, 0, 0], [0, 0, 0]])
    assert(matricize(expr) == expected)
    # ...

    # ...
    expr = v[1]*dx(u[0])
    expr = normalize(expr, basis={V: 'Ni', U: 'Nj'})
    expected = Matrix([[0, Ni*Nj_x, 0], [0, 0, 0], [0, 0, 0]])
    assert(matricize(expr) == expected)
    # ...

    # ...
    expr = dy(v[0])*u[1]
    expr = normalize(expr, basis={V: 'Ni', U: 'Nj'})
    expected = Matrix([[0, 0, 0], [Ni_y*Nj, 0, 0], [0, 0, 0]])
    assert(matricize(expr) == expected)
    # ...

    # ...
    expr = dx(v[1])*dy(u[1])
    expr = normalize(expr, basis={V: 'Ni', U: 'Nj'})
    expected = Matrix([[0, 0, 0], [0, Ni_x*Nj_y, 0], [0, 0, 0]])
    assert(matricize(expr) == expected)
    # ...

    # ...
    expr = div(v) * div(u)
    expr = normalize(expr, basis={V: 'Ni', U: 'Nj'})
    expected = Matrix([[Ni_x*Nj_x, Ni_y*Nj_x, Ni_z*Nj_x],
                       [Ni_x*Nj_y, Ni_y*Nj_y, Ni_z*Nj_y],
                       [Ni_x*Nj_z, Ni_y*Nj_z, Ni_z*Nj_z]])
    assert(matricize(expr) == expected)
    # ...

    # ...
    expr = dot(curl(v), curl(u))
    expr = normalize(expr, basis={V: 'Ni', U: 'Nj'})
    expected = Matrix([[Ni_y*Nj_y + Ni_z*Nj_z, -Ni_x*Nj_y, -Ni_x*Nj_z],
                       [-Ni_y*Nj_x, Ni_x*Nj_x + Ni_z*Nj_z, -Ni_y*Nj_z],
                       [-Ni_z*Nj_x, -Ni_z*Nj_y, Ni_x*Nj_x + Ni_y*Nj_y]])
    assert(matricize(expr) == expected)
    # ...

    # ...
    expr = div(v) * div(u) + dot(curl(v), curl(u))
    expr = normalize(expr, basis={V: 'Ni', U: 'Nj'})
    expected = Matrix([[Ni_x*Nj_x + Ni_y*Nj_y + Ni_z*Nj_z, -Ni_x*Nj_y + Ni_y*Nj_x, -Ni_x*Nj_z + Ni_z*Nj_x],
                       [Ni_x*Nj_y - Ni_y*Nj_x, Ni_x*Nj_x + Ni_y*Nj_y + Ni_z*Nj_z, -Ni_y*Nj_z + Ni_z*Nj_y],
                       [Ni_x*Nj_z - Ni_z*Nj_x, Ni_y*Nj_z - Ni_z*Nj_y, Ni_x*Nj_x + Ni_y*Nj_y + Ni_z*Nj_z]])
    assert(matricize(expr) == expected)
    # ...

    # ...
    expr = c1 * div(v) * div(u) + dot(curl(v), curl(u))
    expr = normalize(expr, basis={V: 'Ni', U: 'Nj'})
    expected =  Matrix([[Ni_x*Nj_x*c1 + Ni_y*Nj_y + Ni_z*Nj_z,
                         -Ni_x*Nj_y + Ni_y*Nj_x*c1,
                         -Ni_x*Nj_z + Ni_z*Nj_x*c1],
                        [Ni_x*Nj_y*c1 - Ni_y*Nj_x,
                         Ni_x*Nj_x + Ni_y*Nj_y*c1 + Ni_z*Nj_z,
                         -Ni_y*Nj_z + Ni_z*Nj_y*c1],
                        [Ni_x*Nj_z*c1 - Ni_z*Nj_x,
                         Ni_y*Nj_z*c1 - Ni_z*Nj_y,
                         Ni_x*Nj_x + Ni_y*Nj_y + Ni_z*Nj_z*c1]])
    assert(matricize(expr) == expected)
    # ...

    # ...
    expr = c1 * div(v) * div(u) + c2 * dot(curl(v), curl(u))
    expr = normalize(expr, basis={V: 'Ni', U: 'Nj'})
    expected = Matrix([[Ni_x*Nj_x*c1 + Ni_y*Nj_y*c2 + Ni_z*Nj_z*c2,
                        -Ni_x*Nj_y*c2 + Ni_y*Nj_x*c1,
                        -Ni_x*Nj_z*c2 + Ni_z*Nj_x*c1],
                       [Ni_x*Nj_y*c1 - Ni_y*Nj_x*c2,
                        Ni_x*Nj_x*c2 + Ni_y*Nj_y*c1 + Ni_z*Nj_z*c2,
                        -Ni_y*Nj_z*c2 + Ni_z*Nj_y*c1],
                       [Ni_x*Nj_z*c1 - Ni_z*Nj_x*c2,
                        Ni_y*Nj_z*c1 - Ni_z*Nj_y*c2,
                        Ni_x*Nj_x*c2 + Ni_y*Nj_y*c2 + Ni_z*Nj_z*c1]])
    assert(matricize(expr) == expected)
    # ...

#    expr = c1 * div(v) * div(u) + c2 * dot(curl(v), curl(u))
#    print('> input         >>> {0}'.format(expr))
#    expr = normalize(expr, basis={V: 'Ni', U: 'Nj'})
#    print('> matricize     >>> {0}'.format(matricize(expr)))
# ...

# ...
def test_bilinear_form_3d_1():
    print('============ test_bilinear_form_3d_1 =============')

    U = H1Space('U', ldim=3)
    V = H1Space('V', ldim=3)

    u = TestFunction(U, name='u')
    v = TestFunction(V, name='v')

    expr = inner(grad(u), grad(v))

    a = BilinearForm((v,u), expr)
    print('> input      >>> {0}'.format(a))

    a_expr = gelatize(a, basis={V: 'Nj', U: 'Ni'})
    print('> gelatized  >>> {0}'.format(a_expr))
    print('')
# ...

# ...
def test_bilinear_form_3d_2():
    print('============ test_bilinear_form_3d_2 =============')

    U = H1Space('U', ldim=3, is_vector=True, shape=3)
    V = H1Space('V', ldim=3, is_vector=True, shape=3)

    u = VectorTestFunction(U, name='u')
    v = VectorTestFunction(V, name='v')

    expr = div(u) * div(v) + 0.2 * dot(u, v)

    a = BilinearForm((v,u), expr)
    print('> input      >>> {0}'.format(a))

    a_expr = gelatize(a, basis={V: 'Nj', U: 'Ni'})
    print('> gelatized  >>> {0}'.format(a_expr))
    print('')
# ...

# ...
def test_bilinear_form_3d_3():
    print('============ test_bilinear_form_3d_3 =============')

    U = H1Space('U', ldim=3, is_vector=True, shape=3)
    V = H1Space('V', ldim=3, is_vector=True, shape=3)

    u = VectorTestFunction(U, name='u')
    v = VectorTestFunction(V, name='v')

    expr = dot(curl(u), curl(v)) + 0.2 * dot(u, v)

    a = BilinearForm((v,u), expr)
    print('> input      >>> {0}'.format(a))

    a_expr = gelatize(a, basis={V: 'Nj', U: 'Ni'})
    print('> gelatized  >>> {0}'.format(a_expr))
    print('')
# ...

## ... TODO debug
#def test_bilinear_form_3d_4():
#    print('============ test_bilinear_form_3d_4 =============')
#
#    U = H1Space('U', ldim=3, is_vector=True, shape=3)
#    V = H1Space('V', ldim=3, is_vector=True, shape=3)
#
#    u = VectorTestFunction(U, name='u')
#    v = VectorTestFunction(V, name='v')
#
#    bx = Constant('bx')
#    by = Constant('by')
#    bz = Constant('bz')
#    b = Tuple(bx, by, bz)
#
#    expr = dot(curl(cross(b,u)), curl(cross(b,v))) + 0.2 * dot(u, v)
#
#    a = BilinearForm((v,u), expr)
#    print('> input      >>> {0}'.format(a))
#
#    a_expr = gelatize(a, basis={V: 'Nj', U: 'Ni'})
#    print('> gelatized  >>> {0}'.format(a_expr))
#    print('')
## ...
#
## ... TODO debug
#def test_bilinear_form_3d_5():
#    print('============ test_bilinear_form_3d_5 =============')
#
#    U = H1Space('U', ldim=3, is_vector=True, shape=3)
#    V = H1Space('V', ldim=3, is_vector=True, shape=3)
#
#    u = VectorTestFunction(U, name='u')
#    v = VectorTestFunction(V, name='v')
#
#    bx = Constant('bx')
#    by = Constant('by')
#    bz = Constant('bz')
#    b = Tuple(bx, by, bz)
#
#    c0,c1,c2 = symbols('c0 c1 c2')
#
#    expr = c0 * dot(u, v) - c1 * div(u) * div(v) + c2 *dot(curl(cross(b,u)), curl(cross(b,v)))
#
#    a = BilinearForm((v,u), expr)
#    print('> input      >>> {0}'.format(a))
#
#    a_expr = gelatize(a, basis={V: 'Nj', U: 'Ni'})
#    print('> gelatized  >>> {0}'.format(a_expr))
#    print('')
## ...

# ...
def test_bilinear_form_3d_10():
    print('============ test_bilinear_form_3d_10 =============')

    U = H1Space('U', ldim=3)
    V = H1Space('V', ldim=3)

    u = TestFunction(U, name='u')
    v = TestFunction(V, name='v')

    u1 = TestFunction(U, name='u1')
    v1 = TestFunction(V, name='v1')

    Ni, Ni_x, Ni_y, Ni_z = symbols('Ni Ni_x Ni_y Ni_z')
    Nj, Nj_x, Nj_y, Nj_z = symbols('Nj Nj_x Nj_y Nj_z')

    c1 = Symbol('c1')
    c2 = Symbol('c2')

    a = BilinearForm((v,u), inner(grad(u), grad(v)))
    b = BilinearForm((v,u), u*v)
    adv = BilinearForm((v,u), dx(u)*v)

    # ...
    expected = Ni*Nj + Ni_x*Nj_x + Ni_y*Nj_y + Ni_z*Nj_z
    assert(gelatize(a + b, basis={V: 'Nj', U: 'Ni'}) == expected)
    # ...

    # ...
    expected = Ni*Nj + Ni_x*Nj_x + Ni_y*Nj_y + Ni_z*Nj_z
    assert(gelatize(a + b, basis={V: 'Nj', U: 'Ni'}) == expected)
    # ...

    # ...
    expected = 2*Ni_x*Nj_x + 2*Ni_y*Nj_y + 2*Ni_z*Nj_z
    assert(gelatize(2 * a, basis={V: 'Nj', U: 'Ni'}) == expected)
    # ...

    # ...
    expected = c1*(Ni_x*Nj_x + Ni_y*Nj_y + Ni_z*Nj_z)
    assert(gelatize(c1*a, basis={V: 'Nj', U: 'Ni'}) == expected)
    # ...

    # ...
    expected = Ni*Nj*c2 + c1*(Ni_x*Nj_x + Ni_y*Nj_y + Ni_z*Nj_z)
    assert(gelatize(c1*a + c2*b, basis={V: 'Nj', U: 'Ni'}) == expected)
    # ...

    # ...
    expected = c1*(Ni_x*Nj_x + Ni_y*Nj_y + Ni_z*Nj_z) + c2*(Ni*Nj + Ni_x*Nj)
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
def test_linear_form_3d_10():
    print('============ test_linear_form_3d_10 =============')

    V = H1Space('V', ldim=3)

    v = TestFunction(V, name='v')
    x,y,z = V.coordinates

    f = Function('f')
    g = Function('g')
    r = Function('r')

    Ni, Ni_x, Ni_y, Ni_z = symbols('Ni Ni_x Ni_y Ni_z')
    Nj, Nj_x, Nj_y, Nj_z = symbols('Nj Nj_x Nj_y Nj_z')

    c1 = Symbol('c1')
    c2 = Symbol('c2')

    bx, by, bz = symbols('bx by bz')
    b = Tuple(bx, by, bz)
    fgr = Tuple(f(x,y,z), g(x,y,z), r(x,y,z))

    a = LinearForm(v, cos(2*pi*x)*cos(4*pi*y)*cos(5*pi*z)*v)

    # ...
    expected = cos(2*pi*x)*cos(4*pi*y)*cos(5*pi*z)*Ni
    assert(gelatize(LinearForm(v, cos(2*pi*x)*cos(4*pi*y)*cos(5*pi*z)*v),
                    basis={V: 'Ni'}) == expected)
    # ...

    # ...
    expected = f(x,y,z)*Ni
    assert(gelatize(LinearForm(v, f(x,y,z)*v),
                    basis={V: 'Ni'}) == expected)
    # ...

    # ...
    expected = bx*Ni_x + by*Ni_y + bz*Ni_z + f(x,y,z)*Ni
    assert(gelatize(LinearForm(v, dot(b, grad(v)) + f(x,y,z)*v),
                    basis={V: 'Ni'}) == expected)
    # ...

    # ...
    expected = f(x,y,z)*Ni_x + g(x,y,z)*Ni_y + r(x,y,z)*Ni_z
    assert(gelatize(LinearForm(v, dot(fgr, grad(v))),
                    basis={V: 'Ni'}) == expected)
    # ...

#    expr = c1*a  + c2*(b + adv)
#    print('> input      >>> {0}'.format(expr))
#    print('> gelatized  >>> {0}'.format(gelatize(expr, basis={V: 'Ni'}) ))
#    print('')
# ...

# ...
def test_linear_form_3d_1():
    print('============ test_linear_form_3d_1 =============')

    V = H1Space('V', ldim=3)

    v = TestFunction(V, name='v')
    x,y,z = V.coordinates

    f = Function('f')

    expr = cos(2*pi*x)*cos(4*pi*y)*cos(5*pi*z)*v
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

    V = H1Space('V', ldim=3, is_block=True, shape=3)
    U = H1Space('U', ldim=3)

    v = VectorTestFunction(V, name='v')
    u = TestFunction(U, name='u')

    Ni, Ni_x, Ni_y, Ni_z = symbols('Ni Ni_x Ni_y Ni_z')
    Nj, Nj_x, Nj_y, Nj_z = symbols('Nj Nj_x Nj_y Nj_z')

    # ...
    expr = v[0]*u
    expr = normalize(expr, basis={V: 'Ni', U: 'Nj'})
    expected = Matrix([[Ni*Nj], [0], [0]])
    assert(matricize(expr) == expected)
    # ...

    # ...
    expr = dot(v, grad(u))
    expr = normalize(expr, basis={V: 'Ni', U: 'Nj'})
    expected = Matrix([[Ni*Nj_x], [Ni*Nj_y], [Ni*Nj_z]])
    assert(matricize(expr) == expected)
    # ...

#    expr = v[0]*u
#    print('> input         >>> {0}'.format(expr))
#    expr = normalize(expr, basis={V: 'Ni', U: 'Nj'})
#    print('> matricize     >>> {0}'.format(matricize(expr)))


# .....................................................
if __name__ == '__main__':
    test_atomize_3d_1()
    test_normalize_3d_1()

    test_atomize_3d_2()
    test_normalize_3d_2()
    test_matricize_3d_2()

#    test_bilinear_form_3d_1()
#    test_bilinear_form_3d_2()
#    test_bilinear_form_3d_3()
##    test_bilinear_form_3d_4()
##    test_bilinear_form_3d_5()

    test_bilinear_form_3d_10()
    test_linear_form_3d_10()

    test_linear_form_3d_1()

    test()
