# coding: utf-8

from sympy import Symbol
from sympy.core.containers import Tuple
from sympy import symbols
from sympy import pi, cos, sin
from sympy import srepr

from symfe.core import dx, dy, dz
from symfe.core import Constant
from symfe.core import Field
from symfe.core import grad, dot, inner, cross, rot, curl, div
from symfe.core import H1Space
from symfe.core import TestFunction
from symfe.core import VectorTestFunction
from symfe.core import BilinearForm

from gelato.core import gelatize


# ...
def test_tensorize_2d_1():
    print('============ test_tensorize_2d_1 =============')

    V = H1Space('V', ldim=2)
    V_0 = H1Space('V_0', ldim=1, coordinates=['x'])
    V_1 = H1Space('V_1', ldim=1, coordinates=['y'])

    v = TestFunction(V, name='v')
    u = TestFunction(V, name='u')

    nx, ny = symbols('nx ny', integer=True)

    m_0, s_0, a_0 = symbols('m_0 s_0 a_0', real=True)
    m_1, s_1, a_1 = symbols('m_1 s_1 a_1', real=True)

    c = Constant('c')

    bx = Constant('bx')
    by = Constant('by')
    b = Tuple(bx, by)

    # ...
    expected = m_0*m_1/(nx*ny)
    assert(gelatize(BilinearForm((v,u), u*v), evaluate=False) == expected)
    # ...

    # ...
    expected = m_1*nx*s_0/ny
    assert(gelatize(BilinearForm((v,u), dx(u)*dx(v)), evaluate=False) == expected)
    # ...

    # ...
    expected = a_1*m_0/nx
    assert(gelatize(BilinearForm((v,u), dy(u) * v), evaluate=False) == expected)
    # ...

    # ...
    expected =  a_0*m_1/ny
    assert(gelatize(BilinearForm((v,u), dx(u) * v), evaluate=False) == expected)
    # ...

    # ...
    expected = m_0*ny*s_1/nx + m_1*nx*s_0/ny
    assert(gelatize(BilinearForm((v,u), dot(grad(v), grad(u))), evaluate=False) == expected)
    # ...

    # ...
    expected = a_0*m_1/ny + a_1*m_0/nx + m_0*ny*s_1/nx + m_1*nx*s_0/ny
    assert(gelatize(BilinearForm((v,u), dot(grad(v), grad(u)) + dx(u)*v + dy(u)*v), evaluate=False) == expected)
    # ...

    # ...
    expected = -bx*a_0*m_1/ny - by*a_1*m_0/nx
    assert(gelatize(BilinearForm((v,u), dot(b, grad(v)) * u), evaluate=False) == expected)
    # ...

    # ...
    expected = bx**2*m_1*nx*s_0/ny - 2*bx*by*a_0*a_1 + by**2*m_0*ny*s_1/nx
    assert(gelatize(BilinearForm((v,u), dot(b, grad(v)) * dot(b, grad(u))), evaluate=False) == expected)
    # ...

    degrees = None
#    degrees = [2, 1]

#    evaluate = True
    evaluate = False

    expr = dot(b, grad(v)) * dot(b, grad(u))

#    expr = BilinearForm((v,u), expr)
#    print('> input     >>> {0}'.format(expr))
#    print('> gelatized >>> {0}'.format(gelatize(expr, degrees, evaluate=evaluate)))
# ...


# .....................................................
if __name__ == '__main__':
    test_tensorize_2d_1()
