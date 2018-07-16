# coding: utf-8

# TODO - transpose of BilinearForm
#      - add unknown status if only one space is given to the BilinearForm
#        => can not be gelatized except if it is called on a couple test/trial
#      - check that a BilinearForm is bilinear (using Constant)
#      - check that a LinearForm is linear
#      - add is_symmetric property for BilinearForm
#      - treat Function atom in atomize, normalize, matricize
#      - treat Field atom in atomize, normalize, matricize

from numpy import zeros

from sympy.core import Basic
from sympy.core import Symbol
from sympy.core import Function
from sympy.core import Expr, Add, Mul
from sympy import S
from sympy.core.containers import Tuple
from sympy import preorder_traversal
from sympy import Indexed, IndexedBase, Matrix, ImmutableDenseMatrix
from sympy.physics.quantum import TensorProduct
from sympy import expand
from sympy import Integer, Float

from gelato.core import _partial_derivatives
from gelato.core import _calculus_operators
from gelato.core import CalculusFunction
from gelato.core import partial_derivative_as_symbol
from gelato.core import sort_partial_derivatives
from gelato.core import get_atom_derivatives
from gelato.core import dx, dy, dz
from gelato.core import Field, Constant
from gelato.core import Dot, Inner, Cross
from gelato.core import Grad, Rot, Curl, Div
from gelato.core import _generic_ops
from gelato.core import _coeffs_registery
from gelato.core import (Dot_1d, Grad_1d, Div_1d)
from gelato.core import (Dot_2d, Cross_2d, Grad_2d, Curl_2d, Rot_2d, Div_2d)
from gelato.core import (Dot_3d, Cross_3d, Grad_3d, Curl_3d, Div_3d)
from gelato.core import BasicSobolevSpace
from gelato.core import TestFunction
from gelato.core import VectorTestFunction



class LinearForm(Expr):
    """

    Examples

    """
    def __new__(cls, test_functions, expr):
        # ...
        if isinstance(test_functions, (TestFunction, VectorTestFunction)):
            test_functions = [test_functions]
            test_functions = Tuple(*test_functions)

        elif not isinstance(test_functions, (tuple, list, Tuple)):
            raise TypeError('Wrong type for test function(s)')
        # ...

        return Basic.__new__(cls, expr, test_functions)

    @property
    def expr(self):
        return self._args[0]

    @property
    def test_functions(self):
        return self._args[1]

    @property
    def ldim(self):
        return self.test_spaces[0].ldim

    @property
    def test_spaces(self):
        return [u.space for u in self.test_functions]

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
        # ...
        if not(len(args) == 1):
            raise ValueError('Expecting one argument')

        test_functions = args[0]
        if isinstance(test_functions, (TestFunction, VectorTestFunction)):
            test_functions = [test_functions]
            test_functions = Tuple(*test_functions)

        elif not isinstance(test_functions, (tuple, list, Tuple)):
            raise TypeError('Wrong type for test function(s)')
        # ...

        # ...
        expr = self.expr
        # ...

        # ... replacing test functions
        d = {}
        for k,v in zip(self.test_functions, test_functions):
            d[k] = v
        expr = expr.subs(d)
        # ...

        # ... replacing trial functions from tmp symbols
        expr = expr.subs(d_tmp)
        # ...

        return LinearForm(test_functions, expr)


class BilinearForm(Expr):
    """

    Examples

    """
    def __new__(cls, test_trial, expr):
        # ...
        # TODO put this in a private function since it is used in __new__ and __call__
        if not isinstance(test_trial, (tuple, list, Tuple)):
            raise TypeError('(test, trial) must be a tuple, list or Tuple')

        if not(len(test_trial) == 2):
            raise ValueError('Expecting a couple (test, trial)')

        test_functions = test_trial[0]
        if isinstance(test_functions, (TestFunction, VectorTestFunction)):
            test_functions = [test_functions]
            test_functions = Tuple(*test_functions)

        elif not isinstance(test_functions, (tuple, list, Tuple)):
            raise TypeError('Wrong type for test function(s)')

        trial_functions = test_trial[1]
        if isinstance(trial_functions, (TestFunction, VectorTestFunction)):
            trial_functions = [trial_functions]
            trial_functions = Tuple(*trial_functions)

        elif not isinstance(trial_functions, (tuple, list, Tuple)):
            raise TypeError('Wrong type for trial function(s)')
        # ...

        return Basic.__new__(cls, expr, test_functions, trial_functions)

    @property
    def expr(self):
        return self._args[0]

    @property
    def test_functions(self):
        return self._args[1]

    @property
    def trial_functions(self):
        return self._args[2]

    @property
    def ldim(self):
        return self.test_spaces[0].ldim

    @property
    def test_spaces(self):
        return [u.space for u in self.test_functions]

    @property
    def trial_spaces(self):
        return [u.space for u in self.trial_functions]

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
        # ...
        # TODO put this in a private function since it is used in __new__ and __call__
        test_trial = args
        if not isinstance(test_trial, (tuple, list, Tuple)):
            raise TypeError('(test, trial) must be a tuple, list or Tuple')

        if not(len(test_trial) == 2):
            raise ValueError('Expecting a couple (test, trial)')

        test_functions = test_trial[0]
        if isinstance(test_functions, (TestFunction, VectorTestFunction)):
            test_functions = [test_functions]
            test_functions = Tuple(*test_functions)

        elif not isinstance(test_functions, (tuple, list, Tuple)):
            raise TypeError('Wrong type for test function(s)')

        trial_functions = test_trial[1]
        if isinstance(trial_functions, (TestFunction, VectorTestFunction)):
            trial_functions = [trial_functions]
            trial_functions = Tuple(*trial_functions)

        elif not isinstance(trial_functions, (tuple, list, Tuple)):
            raise TypeError('Wrong type for trial function(s)')
        # ...

        # ...
        test_trial = [test_functions, trial_functions]
        test_trial = Tuple(*test_trial)

        expr = self.expr
        # ...

        # in order to avoid problems when swapping indices, we need to create
        # temp symbols

        # ...
        d_tmp = {}
        for x in trial_functions:
            name = '{name}_{hash}'.format(name=x.name, hash=abs(hash(x)))
            if isinstance(x, TestFunction):
                X = Symbol(name)
            elif isinstance(x, VectorTestFunction):
                X = IndexedBase(name, shape=x.shape)
            else:
                raise TypeError('Expecting a TestFunction or VectorTestFunction')

            d_tmp[X] = x
        # ...

        # ... replacing trial functions by tmp symbols
        d = {}
        for k,v in zip(self.trial_functions, d_tmp):
            d[k] = v
        expr = expr.subs(d)
        # ...

        # ... replacing test functions
        d = {}
        for k,v in zip(self.test_functions, test_functions):
            d[k] = v
        expr = expr.subs(d)
        # ...

        # ... replacing trial functions from tmp symbols
        expr = expr.subs(d_tmp)
        # ...

        return BilinearForm(test_trial, expr)


# ...
def atomize(expr, dim=None):
    """
    """
    if not isinstance(expr, (Add, Mul,
                             _partial_derivatives, _calculus_operators,
                             TestFunction, VectorTestFunction, Indexed,
                             Field, Constant, Symbol, Function,
                             Integer, Float, Matrix, ImmutableDenseMatrix,
                             list, tuple, Tuple)):
        msg = ('> Wrong input type.')

        raise TypeError(msg, ', given ', expr, type(expr))

#    print('> expr = ', expr, type(expr))

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
#        print('expr    = ', expr)
#        print('vectors = ', vectors)
#        print('coeffs  = ', coeffs )

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

    elif isinstance(expr, Matrix):
        n,m = expr.shape
        lines = []
        for i in range(0, n):
            line = []
            for j in range(0, m):
                line.append(atomize(expr[i,j], dim=dim))
            lines.append(line)
        return Matrix(lines)

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
        others = [i for i in expr.args if not(i in coeffs)]
        for i in others:
            if isinstance(i, Function) and not(isinstance(i, CalculusFunction)):
                coeffs.append(i)
        vectors = [i for i in expr.args if not(i in coeffs)]

        i = S.One
        if coeffs:
            i = Mul(*coeffs)

        j = S.One
        if vectors:
            args = [matricize(i) for i in vectors]

            if len(args) == 1:
                j = args[0]

            elif len(args) == 2:
                # TODO how to be sure about who is left/right? test/trial?
                left = args[0]
                right = args[1]

                if isinstance(left, Matrix) and isinstance(right, Matrix):
                    j = TensorProduct(left.transpose(), right)
                else:
                    j = Mul(*args)

            else:
                raise ValueError('Expecting one or two arguments')


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
def gelatize(a, basis=None, verbose=False):

    if not isinstance(a, (BilinearForm, LinearForm, Add, Mul)):
        raise TypeError('Expecting a BilinearForm, LinearForm, Add or Mul')

    if isinstance(a, Add):
        args = [gelatize(i, basis=basis, verbose=verbose) for i in a.args]
        return Add(*args)

    elif isinstance(a, Mul):
        # a coeff can be a symbol, otherwise the expression c1 * a
        # raises an error
        coeffs  = [i for i in a.args if isinstance(i, _coeffs_registery) or isinstance(i, Symbol)]
        vectors = [i for i in a.args if not(i in coeffs)]

        i = S.One
        if coeffs:
            i = Mul(*coeffs)

        j = S.One
        if vectors:
            args = [gelatize(i, basis=basis, verbose=verbose) for i in vectors]
            j = Mul(*args)

        return Mul(i, j)

    dim = a.ldim
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
