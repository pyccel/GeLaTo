# coding: utf-8

# TODO transpose of BilinearForm:



from sympy.core import Basic
from sympy.core import Expr, Add, Mul
from sympy import S
from sympy.core.containers import Tuple
from sympy import preorder_traversal

from gelato.calculus import sort_partial_derivatives
from gelato.calculus import dx, dy, dz
from gelato.calculus import Field
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
from gelato.fem.core import DerivativeSymbol



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
def gelatize(expr):
    """
    """
    if not isinstance(expr, (BilinearForm, LinearForm, Add, Mul)):
        msg = ('> Wrong input type.'
               '  Expecting BilinearForm, LinearForm, Add, Mul')
        raise TypeError(msg)

    if isinstance(expr, Add):
        args = [gelatize(i) for i in expr.args]
        return Add(*args)

    elif isinstance(expr, Mul):
        coeffs  = [i for i in expr.args if isinstance(i, _coeffs_registery)]
        vectors = [i for i in expr.args if not(i in coeffs)]

        i = S.One
        if coeffs:
            i = Mul(*coeffs)

        j = S.One
        if vectors:
            j = gelatize(Mul(*vectors), evaluate=False)

        return Mul(i, j)

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
def normalize_weak_from(a):
    """
    """
    # ...
    if not isinstance(a, (BilinearForm, LinearForm, Add, Mul)):
        msg = ('> Wrong input type.'
               '  Expecting BilinearForm, LinearForm, Add, Mul')
        raise TypeError(msg)

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
            j = normalize_weak_from(Mul(*vectors), evaluate=False)

        return Mul(i, j)
    # ...

    a = gelatize(a)

    ops = sort_partial_derivatives(a.expr)
    trials = a.trial_functions
    tests = a.test_functions
#    print('* trials = ', trials)
#    print('* tests = ', tests)
#    v1 = trials[0]
#    x = DerivativeSymbol(dx(v1))
#    print(type(v1))
#    print(x)
#    import sys; sys.exit(0)

    expr = a.expr
#    print(len(ops))
#    print('% ops = ', ops)
#    n = 0
    for i in ops:
#        n += 1

        if not(len(i.args) == 1):
            raise ValueError('expecting only one argument for partial derivatives')

        arg = i.args[0]

        # terms like dx(..)
        expr = expr.subs({i: DerivativeSymbol(i)})

#        x = DerivativeSymbol(i)
#        print('****** ', i, x, type(x))
#        if n == 2:
#            expr = expr.subs({i: DerivativeSymbol(i)})
#            import sys; sys.exit(0)
#        print(expr)

    a = BilinearForm(expr,
                     trial_space=a.trial_space,
                     test_space=a.test_space)
    return a
# ...

