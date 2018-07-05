# coding: utf-8

# TODO - transpose of BilinearForm
#      - check that a BilinearForm is bilinear (using Constant)
#      - check that a LinearForm is linear
#      - add is_symmetric property for BilinearForm


from sympy.core import Basic
from sympy.core import Symbol
from sympy.core import Expr, Add, Mul
from sympy import S
from sympy.core.containers import Tuple
from sympy import preorder_traversal
from sympy import Indexed

from gelato.calculus import _partial_derivatives
from gelato.calculus import _calculus_operators
from gelato.calculus import partial_derivative_as_symbol
from gelato.calculus import sort_partial_derivatives
from gelato.calculus import get_atom_derivatives
from gelato.calculus import dx, dy, dz
from gelato.calculus import Field, Constant
from gelato.calculus import Dot, Inner, Cross
from gelato.calculus import Grad, Rot, Curl, Div
from gelato.calculus import _generic_ops
from gelato.calculus import _coeffs_registery
from gelato.calculus import (Dot_1d, Grad_1d, Div_1d)
from gelato.calculus import (Dot_2d, Cross_2d, Grad_2d, Curl_2d, Rot_2d, Div_2d)
from gelato.calculus import (Dot_3d, Cross_3d, Grad_3d, Curl_3d, Div_3d)


from gelato.fem.core import FemSpace
from gelato.fem.core import TestFunction
from gelato.fem.core import TrialFunction
from gelato.fem.core import VectorTestFunction
from gelato.fem.core import VectorTrialFunction



# TODO once BilinearForm is stable
class LinearForm(Expr):
    pass


class BilinearForm(Expr):
    """

    Examples

    """
    def __new__(cls, expr, trial_space=None, test_space=None):
        if not(trial_space is None) and not(test_space is None):
            assert(trial_space.ldim == test_space.ldim)

        return Basic.__new__(cls, expr, trial_space, test_space)

    @property
    def expr(self):
        return self._args[0]

    @property
    def trial_space(self):
        return self._args[1]

    @property
    def test_space(self):
        return self._args[2]

    @property
    def trial_functions(self):
        ls = [a for a in self.expr.free_symbols if isinstance(a, TrialFunction)]
        # no redanduncy
        ls = list(set(ls))

        # ... reorder symbols by name
        # TODO can we do better?
        d = {}
        for i in ls:
            d[i.name] = i
        names = [i.name for i in ls]
        names.sort()
        ls = [d[name] for name in names]
        # ...

        return ls

    @property
    def test_functions(self):
        ls = [a for a in self.expr.free_symbols if isinstance(a, TestFunction)]
        # no redanduncy
        ls = list(set(ls))

        # ... reorder symbols by name
        # TODO can we do better?
        d = {}
        for i in ls:
            d[i.name] = i
        names = [i.name for i in ls]
        names.sort()
        ls = [d[name] for name in names]
        # ...

        return ls

    @property
    def fields(self):
        ls = [a for a in self.expr.free_symbols if isinstance(a, Field)]
        # no redanduncy
        return list(set(ls))

    def _sympystr(self, printer):
        sstr = printer.doprint
        expr = self.expr
        return sstr(expr)

    def __call__(self, *args):
        if not(len(args) == 2):
            raise ValueError('Expecting exactly two arguments')

        # ...
        tests = args[0]
        if isinstance(tests, TestFunction):
            tests = [tests]
            tests = Tuple(*tests)
        elif isinstance(tests, (tuple, list, Tuple)):
            tests = Tuple(*tests)
        else:
            raise TypeError('Wrong type for test functions')
        # ...

        # ...
        trials = args[1]
        if isinstance(trials, TrialFunction):
            trials = [trials]
            trials = Tuple(*trials)
        elif isinstance(trials, (tuple, list, Tuple)):
            trials = Tuple(*trials)
        else:
            raise TypeError('Wrong type for trial functions')
        # ...

        expr = self.expr

        # ... replacing test functions
        d = {}
        for k,v in zip(self.test_functions, tests):
            d[k] = v
        expr = expr.subs(d)
        # ...

        # ... replacing trial functions
        d = {}
        for k,v in zip(self.trial_functions, trials):
            d[k] = v
        expr = expr.subs(d)
        # ...

        return BilinearForm(expr,
                            trial_space=self.trial_space,
                            test_space=self.test_space)


# ...
def gelatize(expr, dim=None):
    """
    """
    if not isinstance(expr, (BilinearForm, LinearForm, Add, Mul,
                             _partial_derivatives, _calculus_operators,
                             TestFunction, VectorTestFunction, Indexed,
                             Field, Constant, Symbol,
                             list, tuple, Tuple)):
        msg = ('> Wrong input type.')

        raise TypeError(msg, ', given ', type(expr))

    # ... compute dim if None
    if dim is None:
        if isinstance(expr, (BilinearForm, LinearForm)):
            dim = expr.test_space.ldim
        else:
            ls = [i for i in expr.free_symbols if isinstance(i, (TestFunction, VectorTestFunction))]
            if len(ls) == 0:
                raise ValueError('Unable to compute dim')

            atom = ls[0]
            dim = atom.space.ldim
    # ...

    if isinstance(expr, (list, tuple, Tuple)):
        args = [gelatize(i, dim=dim) for i in expr]
        return Tuple(*args)

    elif isinstance(expr, Add):
        args = [gelatize(i, dim=dim) for i in expr.args]
        return Add(*args)

    elif isinstance(expr, Mul):
        coeffs  = [i for i in expr.args if isinstance(i, _coeffs_registery)]
        vectors = [i for i in expr.args if not(i in coeffs)]

        i = S.One
        if coeffs:
            i = Mul(*coeffs)

        j = S.One
        if vectors:
            args = [gelatize(i, dim=dim) for i in vectors]
            j = Mul(*args)

        return Mul(i, j)

    elif isinstance(expr, (Dot, Inner, Cross, Grad, Rot, Curl, Div)):
        # if i = Dot(...) then type(i) is Grad
        op = type(expr)
        new  = eval('{0}_{1}d'.format(op, dim))

        args = [gelatize(i, dim=dim) for i in expr.args]
        return new(*args)

    elif isinstance(expr, (BilinearForm, LinearForm)):
        e = gelatize(expr.expr, dim=dim)

        return BilinearForm(e,
                            trial_space=expr.trial_space,
                            test_space=expr.test_space)

    return expr
# ...

# ...
def normalize_weak_from(a, basis=None):
    """
    basis: dict
        for every space we give the name of the basis function symbol
    """
    # ...
    return a
    if not isinstance(a, (BilinearForm, LinearForm, Add, Mul,
                          _partial_derivatives, _calculus_operators, Indexed)):
        msg = ('> Wrong input type.'
               '  Expecting BilinearForm, LinearForm, Add, Mul,'
               '  partial derivatives, calculus operators, Indexed')

        raise TypeError(msg, ', given ', type(a))

    if isinstance(a, Add):
        args = [normalize_weak_from(i) for i in a.args]
        return Add(*args)

    elif isinstance(a, Mul):
        coeffs  = [i for i in a.args if isinstance(i, _coeffs_registery)]
        vectors = [i for i in a.args if not(i in coeffs)]

        i = S.One
        if coeffs:
            i = Mul(*coeffs)

        j = S.One
        if vectors:
            args = [normalize_weak_from(i) for i in vectors]
            j = Mul(*args)

        return Mul(i, j)
    # ...

    # ...
    a = gelatize(a)
    # ...

    # ...
    if isinstance(a, (BilinearForm, LinearForm)):
        ops = sort_partial_derivatives(a.expr)
        expr = a.expr

        trials = a.trial_functions
        tests = a.test_functions
    else:
        ops = sort_partial_derivatives(a)
        expr = a
        trials = [i for i in expr.free_symbols if isinstance(expr, TrialFunction)]
        tests = [i for i in expr.free_symbols if isinstance(expr, TestFunction)]
    # ...

    # ...
    for i in ops:

        if not(len(i.args) == 1):
            raise ValueError('expecting only one argument for partial derivatives')

        arg = i.args[0]

        name = None
        if not(basis is None):
            atom = get_atom_derivatives(i)
            if isinstance(atom, (TestFunction, VectorTestFunction)):
                if atom.space in list(basis.keys()):
                    name = basis[atom.space]
            elif isinstance(atom, Symbol) and atom.is_Indexed:
                base = atom.base
                if base.space in list(basis.keys()):
                    name = basis[base.space]

        # terms like dx(..)
        new = partial_derivative_as_symbol(i, name=name)
        expr = expr.subs({i: new})
    # ...

    # ...
    if not(basis is None):
        for a in trials + tests:
            if isinstance(a, (TestFunction, VectorTestFunction)):
                if a.space in list(basis.keys()):
                    name = basis[a.space]

            elif isinstance(a, Symbol) and a.is_Indexed:
                base = a.base
                if base.space in list(basis.keys()):
                    name = basis[base.space]

            expr = expr.subs({a: Symbol(name)})
    # ...

    if isinstance(a, BilinearForm):
        return BilinearForm(expr,
                            trial_space=a.trial_space,
                            test_space=a.test_space)

    elif isinstance(a, LinearForm):
        raise NotImplementedError('TODO')

    else:
        return expr
# ...
