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
from sympy import S

from symfe.core import dx, dy, dz
from symfe.core import Constant
from symfe.core import Field
from symfe.core import grad, dot, inner, cross, rot, curl, div
from symfe.core import H1Space
from symfe.core import TestFunction
from symfe.core import VectorTestFunction
from symfe.core import BilinearForm
from symfe.core import tensorize
from symfe.core import Mass, Stiffness, Advection, AdvectionT
from symfe.core.basic import _coeffs_registery

# ...
from sympy.core import Add, Mul
from sympy import I as sympy_I

from symfe.core import BilinearForm, BilinearAtomicForm

from gelato.core import (glt_symbol_m,
                         glt_symbol_s,
                         glt_symbol_a,
                         glt_symbol_b)


def gelatize(a, degrees=None, evaluate=True, verbose=False):
    if isinstance(a, BilinearForm) and not(isinstance(a, BilinearAtomicForm)):
        expr = tensorize(a)
        if verbose:
            print('> tensorized = ', expr)
    else:
        expr = a

    if isinstance(expr, Add):
        args = [gelatize(i, degrees=degrees, evaluate=evaluate) for i in expr.args]
        return Add(*args)

    elif isinstance(expr, Mul):
        coeffs  = [i for i in expr.args if isinstance(i, _coeffs_registery)]
        vectors = [i for i in expr.args if not(i in coeffs)]

        i = S.One
        if coeffs:
            i = Mul(*coeffs)

        j = S.One
        if vectors:
            args = [gelatize(i, degrees=degrees, evaluate=evaluate) for i in vectors]
            j = Mul(*args)

        return Mul(i, j)

    elif isinstance(expr, BilinearAtomicForm):

        coord = expr.trial_spaces[0].coordinates

        # ... construct the fourier variable and the number of elements
        t_name = 't{}'.format(coord)
        n_name = 'n{}'.format(coord)

        t = Symbol(t_name)
        n = Symbol(n_name, integer=True)
        # ...

        # ...
        _coordinates = ['x', 'y', 'z']
        index = _coordinates.index(str(coord.name))
        # ...

        if evaluate and ( degrees is None ):
            raise ValueError('> degrees must be provided')

        # ... get the degree
        if not( degrees is None ):
            if not isinstance(degrees, (tuple, list, Tuple)):
                degrees = [degrees]

            p = degrees[index]
        # ...

        if isinstance(expr, Mass):
            if evaluate:
                return glt_symbol_m(p, t, n)

            else:
                name = 'm_{}'.format(index)
                symbol = Symbol(name, real=True)
                return symbol / n

        elif isinstance(expr, Stiffness):
            if evaluate:
                return glt_symbol_s(p, t, n)

            else:
                name = 's_{}'.format(index)
                symbol = Symbol(name, real=True)
                return symbol * n

        elif isinstance(expr, Advection):
            if evaluate:
                return   sympy_I * glt_symbol_a(p, t, n)

            else:
                name = 'a_{}'.format(index)
                symbol = Symbol(name, real=True)
                return symbol

        elif isinstance(expr, AdvectionT):
            if evaluate:
                return - sympy_I * glt_symbol_a(p, t, n)

            else:
                name = 'a_{}'.format(index)
                symbol = - Symbol(name, real=True)
                return symbol

        else:
            raise NotImplementedError('TODO')

    return expr
# ...


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
