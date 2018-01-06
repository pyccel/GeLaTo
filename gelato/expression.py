# -*- coding: utf-8 -*-
#
#
# TODO use to_assign and post processing as expression and not latex => helpful
#      for Fortran and Lua (code gen).
"""This module contains different functions to create and treate the GLT symbols."""

import numpy as np

from sympy.core.sympify import sympify
from sympy.simplify.simplify import simplify
from sympy import Symbol
from sympy import Lambda
from sympy import Function
from sympy import lambdify
from sympy import cos
from sympy import sin
from sympy import Rational
from sympy import diff
from sympy import Matrix
from sympy import latex
from sympy import Integral
from sympy import I as sympy_I
from sympy.core import Basic
from sympy.core.singleton import S
from sympy.simplify.simplify import nsimplify
from sympy.utilities.lambdify import implemented_function
from sympy.matrices.dense import MutableDenseMatrix
from sympy import Mul, Add
from sympy import Tuple
from sympy import postorder_traversal
from sympy import preorder_traversal
from sympy import Indexed
from sympy import IndexedBase

from itertools import product

from gelato.calculus import _generic_ops, _partial_derivatives
from gelato.calculus import (Dot_1d, Grad_1d, Div_1d)
from gelato.calculus import (Dot_2d, Cross_2d, Grad_2d, Curl_2d, Rot_2d, Div_2d)
from gelato.calculus import (Dot_3d, Cross_3d, Grad_3d, Curl_3d, Div_3d)

#try:
#    from pyccel.ast.core import Variable
#except:
#    pass


# TODO find a better solution.
#      this code is duplicated in printing.latex
ARGS_x       = ["x", "y", "z"]
ARGS_u       = ["u", "v", "w"]
ARGS_s       = ["s", "ss"]
BASIS_TEST   = "Ni"
BASIS_TRIAL  = "Nj"
BASIS_PREFIX = ["x", "y", "z", "xx", "yy", "zz", "xy", "yz", "xz"]
TOLERANCE    = 1.e-10
#TOLERANCE    = 1.e-4
SETTINGS     = ["glt_integrate", "glt_formatting", "glt_formatting_atoms"]


# ...
_coord_registery = ['x', 'y', 'z']
_basis_registery = ['Ni',
                    'Ni_x', 'Ni_y', 'Ni_z',
                    'Ni_xx', 'Ni_yy', 'Ni_zz',
                    'Ni_xy', 'Ni_yz', 'Ni_zx',
                    'Nj',
                    'Nj_x', 'Nj_y', 'Nj_z',
                    'Ni_xx', 'Ni_yy', 'Ni_zz',
                    'Ni_xy', 'Ni_yz', 'Ni_zx']


# ...
def gelatize(expr, dim):
    # ... in the case of a Lambda expression
    args = None
    if isinstance(expr, Lambda):
        args = expr.variables
        expr = expr.expr
    # ...

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

    if args:
        return Lambda(args, expr)
    else:
        return expr
# ...

# ...
def dict_to_matrix(d, instructions=None, **settings):
    """
    converts a dictionary of expressions to a matrix

    d: dict
        dictionary of expressions

    instructions: list
        a list to keep track of the applied instructions.

    settings: dict
        dictionary for different settings
    """
    # ...
    assert(type(d) == dict)
    # ...

    # ...
    n_rows = 1
    n_cols = 1
    for key, values in list(d.items()):
        if key[0]+1 > n_rows:
            n_rows = key[0] + 1
        if key[1]+1 > n_cols:
            n_cols = key[1] + 1
    # ...

    # ...
    expressions = []
    for i_row in range(0, n_rows):
        row_expr = []
        for i_col in range(0, n_cols):
            _expr = None
            try:
                _expr = d[i_row,i_col]
            except:
                _expr = S.Zero
            row_expr.append(_expr)
        expressions.append(row_expr)
    # ...

    # ...
    expr = Matrix(expressions)
    # ...

    # ... updates the latex expression
    if instructions is not None:
        # ...
        title  = "GLT symbol"
        instructions.append(latex_title_as_paragraph(title))
        # ...

        # ...
        sets = {}
        for key, value in list(settings.items()):
            if not(key == "glt_integrate"):
                sets[key] = value

        instructions.append(glt_latex(expr, **sets))
        # ...
    # ...

    return expr
# ...

# ...
def initialize_weak_form(f, dim):
    if not isinstance(f, Lambda):
        raise TypeError('Expecting a Lambda')

    args = f.variables
    n_args = len(args)
    if (n_args - dim) % 2 == 1:
        raise ValueError('Wrong number of arguments')

    n = (n_args - dim) / 2

    coords = Tuple(*args[:dim])
    tests  = Tuple(*args[dim:dim+n])
    trials = Tuple(*args[dim+n:])

#    print('> coords : {0}'.format(coords))
#    print('> tests  : {0}'.format(tests))
#    print('> trials : {0}'.format(trials))

    test_names  = [str(i) for i in tests]
    trial_names = [str(i) for i in trials]
    coord_names = [str(i) for i in coords]

    d = {}
    d_args = {}
    # TODO must fix the precision for S.Zero?
    for i_test in range(0, n):
        for i_trial in range(0, n):
            d[(i_test, i_trial)] = S.Zero
            d_args[(i_test, i_trial)] = []


    # ...
    def _find_atom(expr, atom):
        """."""
#        if not(isinstance(atom, (Symbol, IndexedBase, Variable))):
#            raise TypeError('Wrong type, given {0}'.format(type(atom)))
        if not(isinstance(atom, (Symbol, IndexedBase))):
            raise TypeError('Wrong type, given {0}'.format(type(atom)))

        if isinstance(expr, (list, tuple, Tuple)):
            ls = [_find_atom(i, atom) for i in expr]
            return np.array(ls).any()

        if isinstance(expr, Add):
            return _find_atom(expr._args, atom)

        if isinstance(expr, Mul):
            return _find_atom(expr._args, atom)

        if isinstance(expr, Function):
            return _find_atom(expr.args, atom)

        if isinstance(expr, Indexed):
            return (str(expr.base) == str(atom))

#        if isinstance(expr, Variable):
#            return (str(expr) == str(atom))

        if isinstance(expr, Symbol):
            return (str(expr) == str(atom))

        return False
    # ...

    # ...
    def _is_vector(expr, atom):
        """."""
#        if not(isinstance(atom, (Symbol, IndexedBase, Variable))):
#            raise TypeError('Wrong type, given {0}'.format(type(atom)))
        if not(isinstance(atom, (Symbol, IndexedBase))):
            raise TypeError('Wrong type, given {0}'.format(type(atom)))

        if isinstance(expr, (list, tuple, Tuple)):
            ls = [_is_vector(i, atom) for i in expr]
            return np.array(ls).any()

        if isinstance(expr, Add):
            return _is_vector(expr._args, atom)

        if isinstance(expr, Mul):
            return _is_vector(expr._args, atom)

        if isinstance(expr, Function):
            return _is_vector(expr.args, atom)

        if isinstance(expr, Indexed):
            return True

        return False
    # ...

    # ... be careful here, we are using side effect on (d, d_args)
    def _decompose(expr):
        if isinstance(expr, Mul):
            for i_test, test in enumerate(tests):
                for i_trial, trial in enumerate(trials):
                    if _find_atom(expr, test) and _find_atom(expr, trial):
                        d[(i_test, i_trial)] += expr
                        d_args[(i_test, i_trial)] = Tuple(test, trial)
        elif isinstance(expr, Add):
            for e in expr._args:
                _decompose(e)
#        else:
#            raise NotImplementedError('given type {0}'.format(type(expr)))

        return d, d_args
    # ...

    expr = f.expr
    expr = expr.expand()
#    expr = expr.subs({Function('Grad'): Grad})
#    expr = expr.subs({Function('Dot'): Dot})

    d, d_args = _decompose(expr)

    d_expr = {}
    for k,expr in d.items():
        args = list(coords)

        found_vector = False
        for u in d_args[k]:
            if _is_vector(expr, u):
                found_vector = True
                for i in range(0, dim):
                    uofi = IndexedBase(str(u))[i]
                    ui = Symbol('{0}{1}'.format(u, i+1))
                    expr = expr.subs(uofi, ui)
                    args += [ui]
            else:
                args += [u]

        d_expr[k] = Lambda(args, expr)
        if found_vector:
            d_expr[k], _infos = initialize_weak_form(d_expr[k], dim)

    if len(d_expr) == 1:
        key = d_expr.keys()[0]
        d_expr = d_expr[key]

    info = {}
    info['coords'] = coords
    info['tests']  = tests
    info['trials'] = trials

    return d_expr, info

# ...
def normalize_weak_from(f):
    """
    Converts an expression using dx, dy, etc to a normal form, where we
    introduce symbols with suffix to define derivatives.

    f: dict, Lambda
        a valid weak formulation in terms of dx, dy etc
    """
    # ...
    if type(f) == dict:
        d_expr = {}
        for key, g in list(f.items()):
            # ...
            d_expr[key] = normalize_weak_from(g)
            # ...

        return dict_to_matrix(d_expr)
    # ...

    # ...
    if not isinstance(f, Lambda):
        raise TypeError('Expecting a Lambda expression')
    # ...

    # ...
    expr = f.expr

    args   = f.variables
    n_args = len(args)
    # ...

    # ...
    coords = [i for i in f.variables if str(i) in _coord_registery]
    dim    = len(coords)

    if (n_args - dim) % 2 == 1:
        raise ValueError('Wrong number of arguments')

    n = (n_args - dim) / 2

    coords = Tuple(*args[:dim])
    tests  = Tuple(*args[dim:dim+n])
    trials = Tuple(*args[dim+n:])

    # ... we first need to find the ordered list of generic operators
    ops = [a for a in preorder_traversal(expr) if isinstance(a, _partial_derivatives)]
    # ...

    # ...
    for i in ops:
        # if i = dx(u) then type(i) is dx
        op = type(i)
        coordinate = op.coordinate
        for a in i.args:
            # ... test functions
            if a in tests:
                expr = expr.subs({i: Symbol('Ni_{0}'.format(coordinate))})

            if isinstance(a, Indexed) and a.base in tests:
                expr = expr.subs({i: Symbol('Ni_{0}'.format(coordinate))})
            # ...

            # ... trial functions
            if a in trials:
                expr = expr.subs({i: Symbol('Nj_{0}'.format(coordinate))})

            if isinstance(a, Indexed) and a.base in trials:
                expr = expr.subs({i: Symbol('Nj_{0}'.format(coordinate))})
            # ...
    # ...

    # ...
    for i in tests:
        expr = expr.subs({i: Symbol('Ni')})

    for i in trials:
        expr = expr.subs({i: Symbol('Nj')})
    # ...

    return expr
# ...

# ...
class weak_formulation(Function):
    """

    Examples
    ========

    """

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
    def eval(cls, *_args):
        """."""

        if not _args:
            return

        # ...
        f   = _args[0]
        dim = _args[1]
        # ...

        # ...
        f, info = initialize_weak_form(f, dim)

        coords = info['coords']
        tests  = info['tests']
        trials = info['trials']

        test_names  = [str(i) for i in tests]
        trial_names = [str(i) for i in trials]
        coord_names = [str(i) for i in coords]
        # ...

        # ...
        expr = normalize_weak_from(f)

        if isinstance(expr, Matrix):
            expressions = []
            nr = expr.shape[0]
            nc = expr.shape[1]
            for ir in range(0, nr):
                for ic in range(0, nc):
                    expressions += [expr[ir,ic]]
            expr = Tuple(*expressions)

            if len(expr) == 1:
                expr = expr[0]
        # ...

        # ... TODO improve
        free_symbols = [str(i) for i in expr.free_symbols]
        free_symbols.sort()

        args  = _coord_registery[:dim]
        args += [i for i in free_symbols if i in _basis_registery]

        args = [Symbol(i) for i in args]
        # ...

        expr = Lambda(args, expr)

        return expr
# ...

# ...
def construct_weak_form(expr, dim, verbose=False, is_block=False):

    if not isinstance(expr, Lambda):
        raise TypeError('Expecting a Lambda expression')

    # ...
    expr = gelatize(expr, dim=dim)

    if verbose:
        print '> gelatized   := {0}'.format(expr)
    # ...

    # ...
    if is_block:
        expr, info = initialize_weak_form(expr, dim=dim)

        if verbose:
            if isinstance(expr, dict):
                print '> temp form   :='
                # for a nice printing, we print the dictionary entries one by one
                for key, value in expr.items():
                    print '\t\t', key, '\t', value
    # ...

    # ...
    expr = normalize_weak_from(expr)

    if verbose:
        print '> normal form := {0}'.format(expr)
    # ...

    return expr
# ...
