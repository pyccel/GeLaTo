# coding: utf-8


from sympy.core import Add, Mul
from sympy import I as sympy_I
from sympy import Symbol
from sympy.core.containers import Tuple
from sympy import S
from sympy.core import Expr, Basic, AtomicExpr

from symfe.core import BilinearForm, BilinearAtomicForm
from symfe.core import tensorize
from symfe.core import Mass, Stiffness, Advection, AdvectionT
from symfe.core.basic import _coeffs_registery

from .glt import (glt_symbol_m,
                  glt_symbol_s,
                  glt_symbol_a,
                  glt_symbol_b)


# ...
def gelatize(a, degrees=None, evaluate=False, verbose=False):
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
            evaluate = True
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


class Glt(AtomicExpr):
    """

    Examples

    """
    _bilinear_form = None
    def __new__(cls, a, degrees=None, evaluate=False):

        if not isinstance(a, BilinearForm):
            raise TypeError('> Expecting a BilinearForm')

        expr = gelatize(a, degrees=degrees, evaluate=evaluate)

        obj = Basic.__new__(cls, expr)
        obj._bilinear_form = a

        return obj

    @property
    def expr(self):
        return self._args[0]

    @property
    def bilinear_form(self):
        return self._bilinear_form
