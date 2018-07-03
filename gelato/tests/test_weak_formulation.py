# coding: utf-8

# TODO transpose of BilinearForm:

import numpy as np

from collections import OrderedDict

from sympy.core import Basic
from sympy.core import Expr, Add, Mul
from sympy import S
from sympy.core.containers import Tuple
from sympy import symbols
from sympy import Symbol
from sympy import Lambda
from sympy import preorder_traversal

from gelato.calculus import get_atom_derivatives
from gelato.calculus import get_index_derivatives
from gelato.calculus import sort_partial_derivatives
from gelato.calculus import dx, dy, dz
from gelato.calculus import LinearOperator
from gelato.calculus import Constant
from gelato.calculus import Field
from gelato.calculus import grad, dot, inner, cross, rot, curl, div
from gelato.calculus import _generic_ops, _partial_derivatives
from gelato.calculus import _coeffs_registery
from gelato.calculus import (Dot_1d, Grad_1d, Div_1d)
from gelato.calculus import (Dot_2d, Cross_2d, Grad_2d, Curl_2d, Rot_2d, Div_2d)
from gelato.calculus import (Dot_3d, Cross_3d, Grad_3d, Curl_3d, Div_3d)


from gelato.fem.core import FemSpace
from gelato.fem.core import TestFunction
from gelato.fem.core import TrialFunction

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

    expr = a.expr
    for i in ops:

        if not(len(i.args) == 1):
            raise ValueError('expecting only one argument for partial derivatives')

        arg = i.args[0]

        # terms like dx(..)
        expr = expr.subs({i: DerivativeSymbol(i)})

    a = BilinearForm(expr,
                     trial_space=a.trial_space,
                     test_space=a.test_space)
    return a
# ...


# ...
def test_bilinear_form_1d_0():
    print('============ test_bilinear_form_1d_0 =============')

    W = FemSpace('W', ldim=1)
    V = FemSpace('V', ldim=1)

    w = TestFunction(W, name='w')
    v = TrialFunction(V, name='v')

    a = BilinearForm(inner(grad(w), grad(v)), trial_space=V, test_space=W)
    b = BilinearForm(w*v, trial_space=V, test_space=W)
    c = a + b
    print('> input         >>> {0}'.format(c))
    print('> gelatized     >>> {0}'.format(gelatize(c)))
    print('> normal form   >>> {0}'.format(normalize_weak_from(c)))
    print('')
# ...

# ...
def test_bilinear_form_1d_1():
    print('============ test_bilinear_form_1d_1 =============')

    W = FemSpace('W', ldim=1)
    V = FemSpace('V', ldim=1)

    w = TestFunction(W, name='w')
    v = TrialFunction(V, name='v')

    expr = inner(grad(w), grad(v))

    a = BilinearForm(expr, trial_space=V, test_space=W)
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

    a = BilinearForm(expr, trial_space=V, test_space=W)
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

    a = BilinearForm(expr, trial_space=V, test_space=W)
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

    a = BilinearForm(expr, trial_space=V, test_space=W)
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

    a = BilinearForm(expr, trial_space=V, test_space=W)
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

    a = BilinearForm(expr, trial_space=V, test_space=W)
    print('> input         >>> {0}'.format(a))
    print('> gelatized     >>> {0}'.format(gelatize(a)))
    print('> normal form   >>> {0}'.format(normalize_weak_from(a)))
    print('')
# ...

# ...
def test_bilinear_form_1d_7():
    print('============ test_bilinear_form_1d_7 =============')

    W = FemSpace('W', ldim=1)
    V = FemSpace('V', ldim=1)

    w1 = TestFunction(W,  name='w1')
    w2 = TestFunction(W,  name='w2')
    v1 = TrialFunction(V, name='v1')
    v2 = TrialFunction(V, name='v2')
    t1 = TestFunction(W,  name='t1')
    t2 = TestFunction(W,  name='t2')

    a = BilinearForm(inner(grad(w1), grad(v1)) + w2*v2, trial_space=V, test_space=W)
    b = BilinearForm(w1*v2, trial_space=V, test_space=W)

    ls = [a + b, a((t1,t2), (v1,v2)) + b(t1, v2)]
    for c in ls:
        print('> input         >>> {0}'.format(c))
        print('> gelatized     >>> {0}'.format(gelatize(c)))
        print('> normal form   >>> {0}'.format(normalize_weak_from(c)))
        print('')
# ...

# ... TODO not wokring yet
def test_bilinear_form_1d_8():
    print('============ test_bilinear_form_1d_8 =============')

    W = FemSpace('W', ldim=1)
    V = FemSpace('V', ldim=1)

    w1 = TestFunction(W,  name='w1')
    w2 = TestFunction(W,  name='w2')
    v1 = TrialFunction(V, name='v1')
    v2 = TrialFunction(V, name='v2')
    t1 = TestFunction(W,  name='t1')
    t2 = TestFunction(W,  name='t2')

    a = BilinearForm(w1*v2, trial_space=V, test_space=W)

    ls = [a(t1, v2), a(v2, t1)]
    for c in ls:
        print('> input         >>> {0}'.format(c))
        print('> gelatized     >>> {0}'.format(gelatize(c)))
        print('> normal form   >>> {0}'.format(normalize_weak_from(c)))
        print('')
# ...

# ...
def test_bilinear_form_2d_0():
    print('============ test_bilinear_form_2d_0 =============')

    W = FemSpace('W', ldim=2)
    V = FemSpace('V', ldim=2)

    w = TestFunction(W, name='w')
    v = TrialFunction(V, name='v')

    a = BilinearForm(inner(grad(w), grad(v)), trial_space=V, test_space=W)
    b = BilinearForm(w*v, trial_space=V, test_space=W)

    c = a + b
    print('> input         >>> {0}'.format(c))
    print('> gelatized     >>> {0}'.format(gelatize(c)))
    print('> normal form   >>> {0}'.format(normalize_weak_from(c)))
    print('')

    v1 = TestFunction(V, name='v1')
    u1 = TrialFunction(W, name='u1')

    d = a(v1, u1) + b
    print('> input         >>> {0}'.format(d))
    print('> gelatized     >>> {0}'.format(gelatize(d)))
    print('> normal form   >>> {0}'.format(normalize_weak_from(d)))
    print('')
# ...

# ...
def test_bilinear_form_2d_1():
    print('============ test_bilinear_form_2d_1 =============')

    W = FemSpace('W', ldim=2)
    V = FemSpace('V', ldim=2)

    w = TestFunction(W, name='w')
    v = TrialFunction(V, name='v')

    expr = inner(grad(w), grad(v))

    a = BilinearForm(expr, trial_space=V, test_space=W)
    print('> input         >>> {0}'.format(a))
    print('> gelatized     >>> {0}'.format(gelatize(a)))
    print('> normal form   >>> {0}'.format(normalize_weak_from(a)))
    print('')
# ...

# ...
def test_bilinear_form_2d_2():
    print('============ test_bilinear_form_2d_2 =============')

    W = FemSpace('W', ldim=2)
    V = FemSpace('V', ldim=2)

    w = TestFunction(W, name='w')
    v = TrialFunction(V, name='v')

    expr = cross(curl(w), curl(v)) + 0.2 * w * v

    a = BilinearForm(expr, trial_space=V, test_space=W)
    print('> input         >>> {0}'.format(a))
    print('> gelatized     >>> {0}'.format(gelatize(a)))
    print('> normal form   >>> {0}'.format(normalize_weak_from(a)))
    print('')
# ...

# ...
def test_bilinear_form_2d_3():
    print('============ test_bilinear_form_2d_3 =============')

    W = FemSpace('W', ldim=2)
    V = FemSpace('V', ldim=2)

    w = TestFunction(W, name='w')
    v = TrialFunction(V, name='v')

    bx = Constant('bx')
    by = Constant('by')
    b = Tuple(bx, by)

    expr = 0.2 * w * v + dot(b, grad(v)) * w

    a = BilinearForm(expr, trial_space=V, test_space=W)
    print('> input         >>> {0}'.format(a))
    print('> gelatized     >>> {0}'.format(gelatize(a)))
    print('> normal form   >>> {0}'.format(normalize_weak_from(a)))
    print('')
# ...

# ... TODO debug
def test_bilinear_form_2d_4():
    print('============ test_bilinear_form_2d_4 =============')

    W = FemSpace('W', ldim=2)
    V = FemSpace('V', ldim=2)

    w = TestFunction(W, name='w')
    v = TrialFunction(V, name='v')

    expr = rot(w) * rot(v) + div(w) * div(v) + 0.2 * dot(w, v)

    a = BilinearForm(expr, trial_space=V, test_space=W)
    print('> input         >>> {0}'.format(a))
    print('> gelatized     >>> {0}'.format(gelatize(a)))
    print('> normal form   >>> {0}'.format(normalize_weak_from(a)))
    print('')
# ...

# ...
def test_bilinear_form_3d_0():
    print('============ test_bilinear_form_3d_0 =============')

    W = FemSpace('W', ldim=3)
    V = FemSpace('V', ldim=3)

    w = TestFunction(W, name='w')
    v = TrialFunction(V, name='v')

    a = BilinearForm(inner(grad(w), grad(v)), trial_space=V, test_space=W)
    b = BilinearForm(w*v, trial_space=V, test_space=W)
    c = a + b
    print('> input         >>> {0}'.format(c))
    print('> gelatized     >>> {0}'.format(gelatize(c)))
    print('> normal form   >>> {0}'.format(normalize_weak_from(c)))
    print('')
# ...

# ...
def test_bilinear_form_3d_1():
    print('============ test_bilinear_form_3d_1 =============')

    W = FemSpace('W', ldim=3)
    V = FemSpace('V', ldim=3)

    w = TestFunction(W, name='w')
    v = TrialFunction(V, name='v')

    expr = inner(grad(w), grad(v))

    a = BilinearForm(expr, trial_space=V, test_space=W)
    print('> input         >>> {0}'.format(a))
    print('> gelatized     >>> {0}'.format(gelatize(a)))
    print('> normal form   >>> {0}'.format(normalize_weak_from(a)))
    print('')
# ...

# ... TODO debug
def test_bilinear_form_3d_2():
    print('============ test_bilinear_form_3d_2 =============')

    W = FemSpace('W', ldim=3)
    V = FemSpace('V', ldim=3)

    w = TestFunction(W, name='w')
    v = TrialFunction(V, name='v')

    expr = div(w) * div(v) + 0.2 * dot(w, v)

    a = BilinearForm(expr, trial_space=V, test_space=W)
    print('> input         >>> {0}'.format(a))
    print('> gelatized     >>> {0}'.format(gelatize(a)))
    print('> normal form   >>> {0}'.format(normalize_weak_from(a)))
    print('')
# ...

# ... TODO debug
def test_bilinear_form_3d_3():
    print('============ test_bilinear_form_3d_3 =============')

    W = FemSpace('W', ldim=3)
    V = FemSpace('V', ldim=3)

    w = TestFunction(W, name='w')
    v = TrialFunction(V, name='v')

    expr = dot(curl(w), curl(v)) + 0.2 * dot(w, v)

    a = BilinearForm(expr, trial_space=V, test_space=W)
    print('> input         >>> {0}'.format(a))
    print('> gelatized     >>> {0}'.format(gelatize(a)))
    print('> normal form   >>> {0}'.format(normalize_weak_from(a)))
    print('')
# ...

# ... TODO debug
def test_bilinear_form_3d_4():
    print('============ test_bilinear_form_3d_4 =============')

    W = FemSpace('W', ldim=3)
    V = FemSpace('V', ldim=3)

    w = TestFunction(W, name='w')
    v = TrialFunction(V, name='v')

    expr = dot(curl(cross(b,w)), curl(cross(b,v))) + 0.2 * dot(w, v)

    a = BilinearForm(expr, trial_space=V, test_space=W)
    print('> input         >>> {0}'.format(a))
    print('> gelatized     >>> {0}'.format(gelatize(a)))
    print('> normal form   >>> {0}'.format(normalize_weak_from(a)))
    print('')
# ...

# ... TODO debug
def test_bilinear_form_3d_5():
    print('============ test_bilinear_form_3d_5 =============')

    W = FemSpace('W', ldim=3)
    V = FemSpace('V', ldim=3)

    w = TestFunction(W, name='w')
    v = TrialFunction(V, name='v')

    bx = Constant('bx')
    by = Constant('by')
    bz = Constant('bz')
    b = Tuple(bx, by, bz)

    c0,c1,c2 = symbols('c0 c1 c2')

    expr = c0 * dot(w, v) - c1 * div(w) * div(v) + c2 *dot(curl(cross(b,w)), curl(cross(b,v)))

    a = BilinearForm(expr, trial_space=V, test_space=W)
    print('> input         >>> {0}'.format(a))
    print('> gelatized     >>> {0}'.format(gelatize(a)))
    print('> normal form   >>> {0}'.format(normalize_weak_from(a)))
    print('')
# ...

# .....................................................
if __name__ == '__main__':
    test_bilinear_form_1d_0()
    test_bilinear_form_1d_1()
    test_bilinear_form_1d_2()
    test_bilinear_form_1d_3()
    test_bilinear_form_1d_4()
    test_bilinear_form_1d_5()
    test_bilinear_form_1d_6()
    test_bilinear_form_1d_7()
#    test_bilinear_form_1d_8()

    test_bilinear_form_2d_0()
    test_bilinear_form_2d_1()
    test_bilinear_form_2d_2()
    test_bilinear_form_2d_3()
#    test_bilinear_form_2d_4()

    test_bilinear_form_3d_0()
    test_bilinear_form_3d_1()
#    test_bilinear_form_3d_1()
#    test_bilinear_form_3d_2()
#    test_bilinear_form_3d_3()
#    test_bilinear_form_3d_4()
