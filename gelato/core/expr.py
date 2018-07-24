# coding: utf-8

# TODO - transpose of BilinearForm
#      - add unknown status if only one space is given to the BilinearForm
#        => can not be gelatized except if it is called on a couple test/trial
#      - check that a BilinearForm is bilinear (using Constant)
#      - check that a LinearForm is linear
#      - add is_symmetric property for BilinearForm
#      - treat Function atom in atomize, normalize, matricize
#      - treat Field atom in atomize, normalize, matricize
#      - shall we matricize a FunctionForm or not?

from numpy import zeros

from sympy.core import Basic
from sympy.core import Symbol
from sympy.core import Function
from sympy.core import Expr, Add, Mul, Pow
from sympy import S
from sympy.core.containers import Tuple
from sympy import preorder_traversal
from sympy import Indexed, IndexedBase, Matrix, ImmutableDenseMatrix
from sympy.physics.quantum import TensorProduct
from sympy import expand
from sympy import Integer, Float
from sympy.core.expr import AtomicExpr

from .derivatives import _partial_derivatives
from .derivatives import partial_derivative_as_symbol
from .derivatives import sort_partial_derivatives
from .derivatives import get_atom_derivatives
from .derivatives import dx, dy, dz
from .derivatives import (Grad_1d, Div_1d,
                          Grad_2d, Curl_2d, Rot_2d, Div_2d,
                          Grad_3d, Curl_3d, Div_3d)

from .basic import _coeffs_registery
from .basic import CalculusFunction
from .basic import Field, Constant

from .generic import Dot, Inner, Cross
from .generic import Grad, Rot, Curl, Div
from .generic import _generic_ops

from .algebra import (Dot_1d, Dot_2d, Cross_2d, Dot_3d, Cross_3d)

from .space import BasicSobolevSpace
from .space import TestFunction
from .space import VectorTestFunction


class BasicForm(Expr):
    pass


# TODO we should check that the only free symbols are fields, constants or coordinates
class FunctionForm(BasicForm):
    """

    Examples

    """
    _ldim = None
    _coordinates = None
    def __new__(cls, expr, coordinates=None):

        # ... check that there are no test functions in the expression
        ls = [a for a in expr.free_symbols if isinstance(a, (TestFunction, VectorTestFunction))]
        if not(len(ls) == 0):
            raise TypeError('Cannot use test functions in FunctionForm')
        # ...

        # ... compute dim from fields if available
        ls = [a for a in expr.free_symbols if isinstance(a, Field)]
        if ls:
            F = ls[0]
            ldim = F.space.ldim

            if coordinates is None:
                coordinates = F.space.coordinates

        else:
            if coordinates is None:
                raise ValueError('> Coordinates must be provided if the expression has no fields')

            ldim = len(coordinates)
            if ldim == 1:
                coordinates = coordinates[0]
        # ...

        obj = Basic.__new__(cls, expr)
        obj._ldim = ldim
        obj._coordinates = coordinates

        return obj

    @property
    def expr(self):
        return self._args[0]

    @property
    def ldim(self):
        return self._ldim

    @property
    def coordinates(self):
        return self._coordinates

    @property
    def fields(self):
        ls = [a for a in self.expr.free_symbols if isinstance(a, Field)]
        # no redanduncy
        return sorted(list(set(ls)))

    @property
    def constants(self):
        ls = [a for a in self.expr.free_symbols if isinstance(a, Constant)]
        # no redanduncy
        return list(set(ls))

    def _sympystr(self, printer):
        sstr = printer.doprint
        expr = self.expr
        return sstr(expr)

    # TODO how to implement this?
    def __call__(self, *args):
        raise NotImplementedError('')


class LinearForm(BasicForm):
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
        return sorted(list(set(ls)))

    @property
    def constants(self):
        ls = [a for a in self.expr.free_symbols if isinstance(a, Constant)]
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


class BilinearForm(BasicForm):
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
        return sorted(list(set(ls)))

    @property
    def constants(self):
        ls = [a for a in self.expr.free_symbols if isinstance(a, Constant)]
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


class BilinearAtomicForm(BilinearForm, AtomicExpr):
    """

    Examples

    """
    _name = None

    @property
    def name(self):
        return self._name

    def _sympystr(self, printer):
        sstr = printer.doprint
        name = sstr(self.name)

        test = [sstr(i) for i in self.test_functions]
        test = ','.join(i for i in test)

        trial = [sstr(i) for i in self.trial_functions]
        trial = ','.join(i for i in trial)

        return '{name}({test},{trial})'.format(name=name, trial=trial, test=test)

class Mass(BilinearAtomicForm):
    """

    Examples

    """
    _name = 'Mass'
    def __new__(cls, test, trial):

        test_trial = [test, trial]
        expr = test * trial

        return BilinearForm.__new__(cls, test_trial, expr)

class Stiffness(BilinearAtomicForm):
    """

    Examples

    """
    _name = 'Stiffness'
    def __new__(cls, test, trial):

        test_trial = [test, trial]

        coordl = test.space.coordinates.name
        coordr = trial.space.coordinates.name
        if not(coordl == coordr):
            raise ValueError('> Incompatible coordinates')

        ops = {'x': dx, 'y': dy, 'z': dz}
        d = ops[coordl]

        expr = d(test) * d(trial)

        return BilinearForm.__new__(cls, test_trial, expr)

class Advection(BilinearAtomicForm):
    """

    Examples

    """
    _name = 'Advection'
    def __new__(cls, test, trial):

        test_trial = [test, trial]

        coordl = test.space.coordinates.name
        coordr = trial.space.coordinates.name
        if not(coordl == coordr):
            raise ValueError('> Incompatible coordinates')

        ops = {'x': dx, 'y': dy, 'z': dz}
        d = ops[coordl]

        expr = test * d(trial)

        return BilinearForm.__new__(cls, test_trial, expr)

class AdvectionT(BilinearAtomicForm):
    """

    Examples

    """
    _name = 'Advection'
    def __new__(cls, test, trial):

        test_trial = [test, trial]

        coordl = test.space.coordinates.name
        coordr = trial.space.coordinates.name
        if not(coordl == coordr):
            raise ValueError('> Incompatible coordinates')

        ops = {'x': dx, 'y': dy, 'z': dz}
        d = ops[coordl]

        expr = d(test) * trial

        return BilinearForm.__new__(cls, test_trial, expr)

# ...
def atomize(expr, dim=None):
    """
    """
#    if not isinstance(expr, (Add, Mul,
#                             _partial_derivatives, _generic_ops,
#                             TestFunction, VectorTestFunction, Indexed,
#                             Field, Constant, Symbol, Function,
#                             Integer, Float, Matrix, ImmutableDenseMatrix,
#                             list, tuple, Tuple)):
#        msg = ('> Wrong input type.')
#
#        raise TypeError(msg, ', given ', expr, type(expr))

    if not isinstance(expr, (Expr,
                             _partial_derivatives, _generic_ops,
                             TestFunction, VectorTestFunction, Indexed,
                             Field, Constant, Symbol, Function,
                             Integer, Float, Matrix, ImmutableDenseMatrix,
                             list, tuple, Tuple)):
        msg = ('> Wrong input type.')

        raise TypeError(msg, ', given ', expr, type(expr))

#    print('> expr [atomize] = ', expr, type(expr))

    # ... compute dim if None
    if dim is None:
        ls = [i for i in expr.free_symbols if isinstance(i, (TestFunction,
                                                             VectorTestFunction,
                                                             Field))]

        if ls:
            atom = ls[0]
            if atom.space is None:
                raise ValueError('Expecting atom to be associated to a space')

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

    elif isinstance(expr, Pow):

        b = atomize(expr.base, dim=dim)
        e = expr.exp

        return Pow(b, e)

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
def normalize(expr, basis=None, enable_fields=False):
    """
    must be applied after calling atomize

    basis: dict
        for every space we give the name of the basis function symbol
    """
    # ... compute dim
    ls = [i for i in expr.free_symbols if isinstance(i, (TestFunction,
                                                         VectorTestFunction,
                                                         Field))]

    if ls:
        atom = ls[0]
        if atom.space is None:
            raise ValueError('Expecting atom to be associated to a space')

        dim = atom.space.ldim
    # ...

    # ...
    expr = atomize(expr)
    # ...

#    print('> expr [normalize] = ', expr, type(expr))

    if isinstance(expr, (list, tuple, Tuple)):
        args = [normalize(i, basis=basis, enable_fields=enable_fields) for i in expr]
        return Tuple(*args)

    elif isinstance(expr, Add):
        args = [normalize(i, basis=basis, enable_fields=enable_fields) for i in expr.args]
        return Add(*args)

    elif isinstance(expr, Mul):
        coeffs  = [i for i in expr.args if isinstance(i, _coeffs_registery)]
        vectors = [i for i in expr.args if not(i in coeffs)]

        i = S.One
        if coeffs:
            i = Mul(*coeffs)

        j = S.One
        if vectors:
            args = [normalize(i, basis=basis, enable_fields=enable_fields) for i in vectors]
            j = Mul(*args)

        return Mul(i, j)

    elif isinstance(expr, Pow):

        b = normalize(expr.base, basis=basis, enable_fields=enable_fields)
        e = expr.exp

        return Pow(b, e)

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
            if  enable_fields or not isinstance(arg, Field):
                new = partial_derivative_as_symbol(i, name=name, dim=dim)
                expr = expr.subs({i: new})
        # ...

        return expr

    elif isinstance(expr, (TestFunction, VectorTestFunction)):
        if not(basis is None):
            if expr.space in list(basis.keys()):
                name = basis[expr.space]
                return Symbol(name)

#    elif isinstance(expr, Symbol) and expr.is_Indexed:
#        base = expr.base
#        if base.space in list(basis.keys()):
#            name = basis[base.space]
#            return Symbol(name)

    elif isinstance(expr, Indexed):
        if not(basis is None):
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

    if not isinstance(a, (BasicForm, Add, Mul)):
        raise TypeError('Expecting a BasicForm, Add or Mul')

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

    # TODO is it ok to keep this?
    if isinstance(a, FunctionForm):
        return expr

    expr = matricize(expr)
    if verbose:
        print('> matricized >>> {0}'.format(expr))

    return expr
# ...

# TODO - get dim from atoms
#      - check coefficinets/functions
def _tensorize_core(expr, dim, tests, trials):

    if isinstance(expr, Add):
        args = [_tensorize_core(i, dim, tests, trials) for i in expr.args]
        return Add(*args)

    elif isinstance(expr, Mul):
        args = expr.args

        d_atoms = {}
        _coordinates = ['x', 'y', 'z']
        test_trial = list(tests) + list(trials)
        for a in test_trial:
            d_atoms[a] = []

            new = S.One
            for i in range(0, dim):
                coord = _coordinates[i]
                Vi = BasicSobolevSpace('V_{}'.format(i),
                                       ldim=1,
                                       coordinates=[coord])

                ai = TestFunction(Vi, '{test}{i}'.format(test=a.name, i=i))
                d_atoms[a].append(ai)

                new *= ai

            expr = expr.subs({a: new})

        # ...
        # TODO - improve this later
        #      - must distinguish between test/trial
        assert(len(test_trial) == 2)

        # u :: test
        # v :: trial
        u = tests[0]
        v = trials[0]

        ops = {'x': dx, 'y': dy, 'z': dz}
        for ui,vi in zip(d_atoms[u], d_atoms[v]):
            coord = ui.space.coordinates.name
            d = ops[coord]

            # ... Mass
            old = ui * vi
            new = Mass(ui,vi)

            expr = expr.subs({old: new})
            # ...

            # ... Stiffness
            old = d(ui) * d(vi)
            new = Stiffness(ui,vi)

            expr = expr.subs({old: new})
            # ...

            # ... Advection
            old = ui * d(vi)
            new = Advection(ui,vi)

            expr = expr.subs({old: new})
            # ...

            # ... Transpose of Advection
            old = d(ui) * vi
            new = AdvectionT(ui,vi)

            expr = expr.subs({old: new})
            # ...
        # ...

    return expr

def tensorize(a):

    if not isinstance(a, BilinearForm):
        raise TypeError('Expecting a BilinearForm')

    dim = a.ldim
    expr = a.expr
    tests = a.test_functions
    trials = a.trial_functions

    expr = atomize(expr)
    expr = _tensorize_core(expr, dim, tests, trials)

    return expr
