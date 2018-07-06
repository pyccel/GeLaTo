# coding: utf-8

# TODO - transpose of BilinearForm
#      - check that a BilinearForm is bilinear (using Constant)
#      - check that a LinearForm is linear
#      - add is_symmetric property for BilinearForm
#      - treat Function atom in atomize, normalize, matricize

from numpy import zeros

from sympy.core import Basic
from sympy.core import Symbol
from sympy.core import Expr, Add, Mul
from sympy import S
from sympy.core.containers import Tuple
from sympy import preorder_traversal
from sympy import Indexed, IndexedBase, Matrix
from sympy.physics.quantum import TensorProduct
from sympy import expand

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
    def __new__(cls, spaces, expr):
        if not isinstance(spaces, (tuple, list, Tuple)):
            raise NotImplementedError('dual space not yet available')

        if not(len(spaces) == 2):
            raise ValueError('Expecting two spaces')

        test_space = spaces[0]
        trial_space = spaces[0]

        if not(trial_space.ldim == test_space.ldim):
            raise ValueError('Incompatible logical dimension between test and trial spaces')

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
def atomize(expr, dim=None):
    """
    """
    if not isinstance(expr, (Add, Mul,
                             _partial_derivatives, _calculus_operators,
                             TestFunction, VectorTestFunction, Indexed,
                             Field, Constant, Symbol,
                             list, tuple, Tuple)):
        msg = ('> Wrong input type.')

        raise TypeError(msg, ', given ', type(expr))

#    print('> expr = ', expr)

    # ... compute dim if None
    if dim is None:
        ls = [i for i in expr.free_symbols if isinstance(i, (TestFunction, VectorTestFunction))]

        if ls:
            atom = ls[0]
            dim = atom.space.ldim
    # ...

    if isinstance(expr, (list, tuple, Tuple)):
        args = [atomize(i, dim=dim) for i in expr]
        return Tuple(*args)

    elif isinstance(expr, Add):
        args = [atomize(i, dim=dim) for i in expr.args]
        return Add(*args)

    elif isinstance(expr, Mul):
        coeffs  = [i for i in expr.args if isinstance(i, _coeffs_registery)]
        vectors = [i for i in expr.args if not(i in coeffs)]

        i = S.One
        if coeffs:
            i = Mul(*coeffs)

        j = S.One
        if vectors:
            args = [atomize(i, dim=dim) for i in vectors]
            j = Mul(*args)

        return Mul(i, j)

    elif isinstance(expr, (Dot, Inner, Cross, Grad, Rot, Curl, Div)):
        # if i = Dot(...) then type(i) is Grad
        op = type(expr)
        new  = eval('{0}_{1}d'.format(op, dim))

        args = [atomize(i, dim=dim) for i in expr.args]
        return new(*args)

    return expr
# ...

# ... TODO remove this function
def normalize_weak_from(a, basis=None):
    raise NotImplementedError('TODO')
# ...


# ...
def normalize(expr, basis=None):
    """
    must be applied after calling atomize

    basis: dict
        for every space we give the name of the basis function symbol
    """
    # ... compute dim
    ls = [i for i in expr.free_symbols if isinstance(i, (TestFunction, VectorTestFunction))]

    if ls:
        atom = ls[0]
        dim = atom.space.ldim
    # ...

    # ...
    expr = atomize(expr)
    # ...

#    print('> expr = ', expr, type(expr))

    if isinstance(expr, (list, tuple, Tuple)):
        args = [normalize(i, basis=basis) for i in expr]
        return Tuple(*args)

    elif isinstance(expr, Add):
        args = [normalize(i, basis=basis) for i in expr.args]
        return Add(*args)

    elif isinstance(expr, Mul):
        coeffs  = [i for i in expr.args if isinstance(i, _coeffs_registery)]
        vectors = [i for i in expr.args if not(i in coeffs)]

        i = S.One
        if coeffs:
            i = Mul(*coeffs)

        j = S.One
        if vectors:
            args = [normalize(i, basis=basis) for i in vectors]
            j = Mul(*args)

        return Mul(i, j)

    elif isinstance(expr, _partial_derivatives):
        ops = sort_partial_derivatives(expr)

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

                elif isinstance(atom, Indexed):
                    base = atom.base
                    if base.space in list(basis.keys()):
                        name = basis[base.space]

            # terms like dx(..)
            new = partial_derivative_as_symbol(i, name=name, dim=dim)
            expr = expr.subs({i: new})
        # ...

        return expr

    elif isinstance(expr, (TestFunction, VectorTestFunction)):
        if expr.space in list(basis.keys()):
            name = basis[expr.space]
            return Symbol(name)

#    elif isinstance(expr, Symbol) and expr.is_Indexed:
#        base = expr.base
#        if base.space in list(basis.keys()):
#            name = basis[base.space]
#            return Symbol(name)

    elif isinstance(expr, Indexed):
        base = expr.base
        if base.space in list(basis.keys()):
            name = basis[base.space]
            indices = expr.indices
            return IndexedBase(name, shape=dim)[indices]

    return expr
# ...

# ...
def matricize(expr):
    """
    must be applied after calling normalize
    """

    # ... we need first to expand the expression
    expr = expand(expr)
    # ...

#    print('> expr = ', expr, type(expr))

    if isinstance(expr, (list, tuple, Tuple)):
        args = [matricize(i) for i in expr]
        return Tuple(*args)

    elif isinstance(expr, Add):
        args = [matricize(i) for i in expr.args]
        # we cannot return Add(*args)
        # since it gives the error:
        # TypeError: cannot add <class 'sympy.matrices.immutable.ImmutableDenseMatrix'> and <class 'sympy.core.numbers.Zero'>
        # when args are of type Matrix
        r = args[0]
        for a in args[1:]:
            r += a
        return r

    elif isinstance(expr, Mul):
        # a coeff can be a symbol, otherwise the expression rot(v) * rot(u) + c * div(v) * div(u)
        # raises an error
        coeffs  = [i for i in expr.args if isinstance(i, _coeffs_registery) or isinstance(i, Symbol)]
        vectors = [i for i in expr.args if not(i in coeffs)]

        i = S.One
        if coeffs:
            i = Mul(*coeffs)

        j = S.One
        if vectors:
            args = [matricize(i) for i in vectors]
            if not(len(args) == 2):
                print(args)
                raise ValueError('Expecting exactly 2 arguments')

            # TODO how to be sure about who is left/right? test/trial?
            left = args[0]
            right = args[1]

            if isinstance(left, Matrix) and isinstance(right, Matrix):
                j = TensorProduct(left.transpose(), right)
            else:
                j = Mul(*args)

        # we cannot return Mul(i, j)
        # since it gives the error:
        # TypeError: cannot add <class 'sympy.matrices.immutable.ImmutableDenseMatrix'> and <class 'sympy.core.mul.Mul'>
        # when args are of type Matrix
        return i * j

    elif isinstance(expr, Indexed):
        base = expr.base
        if not(len(expr.indices) == 1):
            raise ValueError('Expecting exactly one index')

        index = expr.indices[0]
        name = '{}'.format(base)

        dim = base.shape[0]

        M = Matrix(zeros(dim))
        M[index] = Symbol(name)

        return M

    return expr
# ...

# ... TODO compute basis if not given
def gelatize(a, basis=None, verbose=True):

    if not isinstance(a, (BilinearForm, LinearForm)):
        raise TypeError('Expecting a BilinearForm or LinearForm')

    dim = a.test_space.ldim
    expr = a.expr

    expr = atomize(expr, dim=dim)
    if verbose:
        print('> atomized   >>> {0}'.format(expr))

    expr = normalize(expr, basis=basis)
    if verbose:
        print('> normalized >>> {0}'.format(expr))

    expr = matricize(expr)
    if verbose:
        print('> matricized >>> {0}'.format(expr))

    return expr

#    if isinstance(a, BilinearForm):
#        return BilinearForm(expr, trial_space=a.trial_space, test_space=a.test_space)
#
#    elif isinstance(a, LinearForm):
#        raise NotImplementedError('not implemented yet for LinearForm')
##        return LinearForm(expr, test_space=a.test_space)
# ...
