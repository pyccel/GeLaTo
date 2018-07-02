# coding: utf-8

import numpy as np

from collections import OrderedDict

from sympy.core import Basic
from sympy import preorder_traversal
from sympy.core.containers import Tuple
from sympy import symbols
from sympy import Symbol
from sympy import Lambda
from sympy import Expr

from gelato.calculus import get_atom_derivatives
from gelato.calculus import get_index_derivatives
from gelato.calculus import sort_partial_derivatives
from gelato.calculus import dx, dy, dz
from gelato.calculus import LinearOperator
from gelato.calculus import Constant
from gelato.calculus import Field
from gelato.calculus import grad, dot, inner
from gelato.calculus import _generic_ops, _partial_derivatives
from gelato.calculus import (Dot_1d, Grad_1d, Div_1d)
from gelato.calculus import (Dot_2d, Cross_2d, Grad_2d, Curl_2d, Rot_2d, Div_2d)
from gelato.calculus import (Dot_3d, Cross_3d, Grad_3d, Curl_3d, Div_3d)


from gelato.fem.core import FemSpace
from gelato.fem.core import TestFunction
from gelato.fem.core import TrialFunction

# TODO once BilinearForm is stable
class LinearForm(Basic):
    pass


class BilinearForm(Basic):
    """

    Examples

    """
    def __new__(cls, expr, name=None, trial_space=None, test_space=None):
        if not(trial_space is None) and not(test_space is None):
            assert(trial_space.ldim == test_space.ldim)

        return Basic.__new__(cls, expr, name, trial_space, test_space)

    @property
    def expr(self):
        return self._args[0]

    @property
    def name(self):
        return self._args[1]

    @property
    def trial_space(self):
        return self._args[2]

    @property
    def test_space(self):
        return self._args[3]

    @property
    def trial_functions(self):
        ls = [a for a in self.expr.free_symbols if isinstance(a, TrialFunction)]
        # no redanduncy
        return list(set(ls))

    @property
    def test_functions(self):
        ls = [a for a in self.expr.free_symbols if isinstance(a, TestFunction)]
        # no redanduncy
        return list(set(ls))

    @property
    def fields(self):
        ls = [a for a in self.expr.free_symbols if isinstance(a, Field)]
        # no redanduncy
        return list(set(ls))

    def _sympystr(self, printer):
        sstr = printer.doprint
        name = self.name
        expr = self.expr
        trials = self.trial_functions
        tests = self.test_functions
        if len(trials) == 1: trials = trials[0]
        if len(tests) == 1: tests = tests[0]

        expr_str = '({tests}, {trials}) -> {expr}'.format(tests=sstr(tests),
                                                          trials=sstr(trials),
                                                          expr=sstr(expr))
        if name is None:
            return expr_str
        else:
            return '{name} := {expr}'.format(name=sstr(self.name), expr=expr_str)


# ...
def gelatize(expr):
    """
    """
    if not isinstance(expr, (BilinearForm, LinearForm)):
        raise TypeError('Expecting a BilinearForm or LinearForm.')

    dim = expr.test_space.ldim

    # ... we first need to find the ordered list of generic operators
    ops = [a for a in preorder_traversal(expr) if isinstance(a, _generic_ops)]
    # ...

    # ...
    for i in ops:
        # if i = Grad(u) then type(i) is Grad
        op = type(i)

        new  = eval('{0}_{1}d'.format(op, dim))
        expr = expr.subs(op, new)
    # ...

    return expr
# ...

# ...
class DerivativeSymbol(Symbol):
    """

    Examples

    """
    def __new__(cls, expr):
        assert(isinstance(expr, _partial_derivatives))

        index = get_index_derivatives(expr)
        var = get_atom_derivatives(expr)

        if not isinstance(var, (TestFunction, TrialFunction, Field)):
            raise TypeError('Expecting TestFunction, TrialFunction, Field')

        return Basic.__new__(cls, var, index)

    @property
    def var(self):
        return self._args[0]

    @property
    def index(self):
        return self._args[1]

    @property
    def name(self):
        code = ''
        for k,n in list(self.index.items()):
            code += k*n

        return '{var}_{code}'.format(var=self.var, code=code)
# ...

# ...
def normalize_weak_from(a):
    """
    """
    # ...
    if type(a) == dict:
        d_expr = {}
        for key, g in list(a.items()):
            # ...
            d_expr[key] = normalize_weak_from(g)
            # ...

        return dict_to_matrix(d_expr)
    # ...

    # ...
    if not isinstance(a, (BilinearForm, LinearForm)):
        raise TypeError('Expecting a BilinearForm or LinearForm.')
    # ...

    a = gelatize(a)

    ops = sort_partial_derivatives(a.expr)
    trials = a.trial_functions
    tests = a.test_functions

    expr = a.expr
    for i in ops:

        if not(len(i.args) == 1):
            raise ValueError('expecting only one argument for partial derivatives')

        arg = i.args[0]

        # terms like dx(..)
        expr = expr.subs({i: DerivativeSymbol(i)})

    return expr
# ...


# ...
def test_bilinear_form_1d_1():
    print('============ test_bilinear_form_1d_1 =============')

    W = FemSpace('W', ldim=1)
    V = FemSpace('V', ldim=1)

    w = TestFunction(W, name='w')
    v = TrialFunction(V, name='v')

    expr = inner(grad(w), grad(v))

    a = BilinearForm(expr, name='a', trial_space=V, test_space=W)
    print('> input         >>> {0}'.format(a))
    print('> gelatized     >>> {0}'.format(gelatize(a)))
    print('> normal form   >>> {0}'.format(normalize_weak_from(a)))
    print('')
# ...

# ...
def test_bilinear_form_1d_2():
    print('============ test_bilinear_form_1d_2 =============')

    W = FemSpace('W', ldim=1)
    V = FemSpace('V', ldim=1)

    w = TestFunction(W, name='w')
    v = TrialFunction(V, name='v')

    b = Constant('b')

    expr = inner(grad(b*w), grad(v))

    a = BilinearForm(expr, name='a', trial_space=V, test_space=W)
    print('> input         >>> {0}'.format(a))
    print('> gelatized     >>> {0}'.format(gelatize(a)))
    print('> normal form   >>> {0}'.format(normalize_weak_from(a)))
    print('')
# ...

# ...
def test_bilinear_form_1d_3():
    print('============ test_bilinear_form_1d_3 =============')

    W = FemSpace('W', ldim=1)
    V = FemSpace('V', ldim=1)

    w = TestFunction(W, name='w')
    v = TrialFunction(V, name='v')

    expr = inner(grad(w), grad(v))

    a = BilinearForm(expr, name='a', trial_space=V, test_space=W)
    print('> input         >>> {0}'.format(a))
    print('> gelatized     >>> {0}'.format(gelatize(a)))
    print('> normal form   >>> {0}'.format(normalize_weak_from(a)))
    print('')
# ...

# ...
def test_bilinear_form_1d_3():
    print('============ test_bilinear_form_1d_3 =============')

    W = FemSpace('W', ldim=1)
    V = FemSpace('V', ldim=1)

    w = TestFunction(W, name='w')
    v = TrialFunction(V, name='v')

    F = Field('F')
    expr = inner(grad(w), grad(v)) + F*w*v

    a = BilinearForm(expr, name='a', trial_space=V, test_space=W)
    print('> input         >>> {0}'.format(a))
    print('> gelatized     >>> {0}'.format(gelatize(a)))
    print('> normal form   >>> {0}'.format(normalize_weak_from(a)))
    print('')
# ...

# ...
def test_bilinear_form_1d_4():
    print('============ test_bilinear_form_1d_4 =============')

    W = FemSpace('W', ldim=1)
    V = FemSpace('V', ldim=1)

    w = TestFunction(W, name='w')
    v = TrialFunction(V, name='v')

    F = Field('F')

    expr = inner(grad(F*w), grad(v))

    a = BilinearForm(expr, name='a', trial_space=V, test_space=W)
    print('> input         >>> {0}'.format(a))
    print('> gelatized     >>> {0}'.format(gelatize(a)))
    print('> normal form   >>> {0}'.format(normalize_weak_from(a)))
    print('')
# ...

# ...
def test_bilinear_form_1d_5():
    print('============ test_bilinear_form_1d_5 =============')

    W = FemSpace('W', ldim=1)
    V = FemSpace('V', ldim=1)

    w = TestFunction(W, name='w')
    v = TrialFunction(V, name='v')

    expr = dx(dx(v))*dx(dx(dx(w))) + w*v

    a = BilinearForm(expr, name='a', trial_space=V, test_space=W)
    print('> input         >>> {0}'.format(a))
    print('> gelatized     >>> {0}'.format(gelatize(a)))
    print('> normal form   >>> {0}'.format(normalize_weak_from(a)))
    print('')
# ...

# ...
def test_bilinear_form_1d_6():
    print('============ test_bilinear_form_1d_6 =============')

    W = FemSpace('W', ldim=1)
    V = FemSpace('V', ldim=1)

    w = TestFunction(W, name='w')
    v = TrialFunction(V, name='v')

    F = Field('F')

    expr = inner(grad(dx(F)*v), grad(w)) + w*v

    a = BilinearForm(expr, name='a', trial_space=V, test_space=W)
    print('> input         >>> {0}'.format(a))
    print('> gelatized     >>> {0}'.format(gelatize(a)))
    print('> normal form   >>> {0}'.format(normalize_weak_from(a)))
    print('')
# ...

# .....................................................
if __name__ == '__main__':
    test_bilinear_form_1d_1()
    test_bilinear_form_1d_2()
    test_bilinear_form_1d_3()
    test_bilinear_form_1d_4()
    test_bilinear_form_1d_5()
    test_bilinear_form_1d_6()
