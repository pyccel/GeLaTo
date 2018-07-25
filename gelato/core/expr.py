# coding: utf-8

from sympy import Function
from sympy.core import Add, Mul
from sympy import I as sympy_I
from sympy import Symbol
from sympy.core.containers import Tuple
from sympy import S
from sympy.core import Expr, Basic, AtomicExpr

from symfe.core import BilinearForm, BilinearAtomicForm
from symfe.core import tensorize
from symfe.core import Mass as MassForm
from symfe.core import Stiffness as StiffnessForm
from symfe.core import Advection as AdvectionForm
from symfe.core import AdvectionT as AdvectionTForm
from symfe.core.basic import _coeffs_registery

from .glt import (Mass,
                  Stiffness,
                  Advection,
                  Bilaplacian)


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
        else:
            p_name = 'p{}'.format(coord)
            p = Symbol(p_name, integer=True)
        # ...

        if isinstance(expr, MassForm):
            symbol = Mass(p, t, evaluate=evaluate)
            return symbol / n

        elif isinstance(expr, StiffnessForm):
            symbol = Stiffness(p, t, evaluate=evaluate)
            return symbol * n

        elif isinstance(expr, AdvectionForm):
            symbol = sympy_I * Advection(p, t, evaluate=evaluate)
            return symbol

        elif isinstance(expr, AdvectionTForm):
            symbol = - sympy_I * Advection(p, t, evaluate=evaluate)
            return symbol

        else:
            raise NotImplementedError('TODO')

    return expr
# ...


class Glt(Function):
    """

    Examples

    """
    _bilinear_form = None
    nargs = None

    def __new__(cls, *args, **options):
        # (Try to) sympify args first

        if options.pop('evaluate', True):
            r = cls.eval(*args)
        else:
            r = None

        if r is None:
            return Basic.__new__(cls, *args, **options)
        else:
            return r

    @classmethod
    def eval(cls, a, degrees=None, evaluate=False):

        if not isinstance(a, BilinearForm):
            raise TypeError('> Expecting a BilinearForm')

        expr = gelatize(a, degrees=degrees, evaluate=evaluate)

        cls._bilinear_form = a

        return expr

    @property
    def bilinear_form(self):
        return self._bilinear_form
