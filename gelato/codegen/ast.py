
from collections import OrderedDict
from itertools import groupby
import string
import random
import numpy as np

from sympy import Basic
from sympy import symbols, Symbol, IndexedBase, Indexed, Function
from sympy import Mul, Add, Tuple
from sympy import Matrix, ImmutableDenseMatrix
from sympy import sqrt as sympy_sqrt
from sympy import S as sympy_S

from pyccel.ast.core import For
from pyccel.ast.core import Assign
from pyccel.ast.core import AugAssign
from pyccel.ast.core import Slice
from pyccel.ast.core import Range
from pyccel.ast.core import FunctionDef
from pyccel.ast.core import FunctionCall
from pyccel.ast.core import Import
from pyccel.ast import Zeros
from pyccel.ast import Import
from pyccel.ast import DottedName
from pyccel.ast import Nil
from pyccel.ast import Len
from pyccel.ast import If, Is, Return
from pyccel.ast import String, Print, Shape
from pyccel.ast import Comment, NewLine
from pyccel.parser.parser import _atomic

from sympde.core import grad
from sympde.core import Constant
from sympde.core import Mapping
from sympde.core import Field
from sympde.core import Covariant, Contravariant
from sympde.core import BilinearForm
from sympde.core.derivatives import _partial_derivatives
from sympde.core.derivatives import get_max_partial_derivatives
from sympde.core.space import FunctionSpace
from sympde.core.space import TestFunction
from sympde.core.space import VectorTestFunction
from sympde.core.space import Trace
from sympde.printing.pycode import pycode  # TODO remove from here
from sympde.core.derivatives import print_expression
from sympde.core.derivatives import get_atom_derivatives
from sympde.core.derivatives import get_index_derivatives
from sympde.core.math import math_atoms_as_str

from gelato.core import gelatize

from spl.fem.splines import SplineSpace

def random_string( n ):
    chars    = string.ascii_uppercase + string.ascii_lowercase + string.digits
    selector = random.SystemRandom()
    return ''.join( selector.choice( chars ) for _ in range( n ) )


class GelatoBasic(Basic):

    def __new__(cls, tag, name=None, prefix=None, debug=False, detailed=False):

        if name is None:
            if prefix is None:
                raise ValueError('prefix must be given')

            name = '{prefix}_{tag}'.format(tag=tag, prefix=prefix)

        obj = Basic.__new__(cls)
        obj._name = name
        obj._tag = tag
        obj._dependencies = []
        obj._debug = debug
        obj._detailed = detailed

        return obj

    @property
    def name(self):
        return self._name

    @property
    def tag(self):
        return self._tag

    @property
    def func(self):
        return self._func

    @property
    def basic_args(self):
        return self._basic_args

    @property
    def dependencies(self):
        return self._dependencies

    @property
    def debug(self):
        return self._debug

    @property
    def detailed(self):
        return self._detailed


# target is used when there are multiple expression (domain/boundaries)
class Kernel(GelatoBasic):

    def __new__(cls, weak_form, discrete_space, name=None):

        if not isinstance(weak_form, BilinearForm):
            raise TypeError('> Expecting a weak formulation')

        tag = random_string( 8 )
        obj = GelatoBasic.__new__(cls, tag, name=name, prefix='kernel')

        obj._weak_form = weak_form
        obj._discrete_space = discrete_space
        obj._eval_fields = None
        obj._eval_mapping = None

        obj._func = obj._initialize()

        return obj

    @property
    def weak_form(self):
        return self._weak_form

    @property
    def discrete_space(self):
        return self._discrete_space

    @property
    def n_rows(self):
        return self._n_rows

    @property
    def n_cols(self):
        return self._n_cols

    @property
    def max_nderiv(self):
        return self._max_nderiv

    @property
    def constants(self):
        return self._constants

    @property
    def fields(self):
        return self._fields

    @property
    def fields_coeffs(self):
        return self._fields_coeffs

    @property
    def mapping_coeffs(self):
        if not self.eval_mapping:
            return ()

        return self.eval_mapping.mapping_coeffs

    @property
    def mapping_values(self):
        if not self.eval_mapping:
            return ()

        return self.eval_mapping.mapping_values

    @property
    def eval_fields(self):
        return self._eval_fields

    @property
    def eval_mapping(self):
        return self._eval_mapping

    def build_arguments(self, data):

        other = data

        if self.mapping_values:
            other = self.mapping_values + other

        if self.constants:
            other = other + self.constants

        return self.basic_args + other

    def _initialize(self):
        weak_form = self.weak_form
        dim       = weak_form.ldim

        # ...
        n_rows = 1 ; n_cols = 1
        if isinstance(weak_form, (Matrix, ImmutableDenseMatrix)):
            n_rows = weak_form.shape[0]
            n_cols = weak_form.shape[1]

        self._n_rows = n_rows
        self._n_cols = n_cols
        # ...

        # ... discrete values
        Vh = self.discrete_space

        n_elements = Vh.ncells
        if isinstance(Vh, SplineSpace):
            degrees    = Vh.degree

        else:
            degrees    = Vh.degrees
        # ...

        # ...
        kernel_expr = gelatize(weak_form, degrees=degrees, n_elements=n_elements)
        # ...

        # ...
        constants = tuple(kernel_expr.atoms(Constant))
        self._constants = constants
        # ...

        # ...
        field_atoms = tuple(kernel_expr.atoms(Field))
        # ...

        # ...
        degrees    = symbols('p1:%d'%(dim+1), integer=True)
        n_elements = symbols('n1:%d'%(dim+1), integer=True)
        tis        = symbols('t1:%d'%(dim+1), real=True)
        arr_tis    = symbols('arr_t1:%d'%(dim+1), cls=IndexedBase)
        indices    = symbols('i1:%d'%(dim+1))
        lengths    = symbols('nt1:%d'%(dim+1))
        ranges     = [Range(lengths[i]) for i in range(dim)]
        # ...

        # ...
        d_symbols = {}
        for i in range(0, n_rows):
            for j in range(0, n_cols):
                mat = IndexedBase('symbol_{i}{j}'.format(i=i,j=j))
                d_symbols[i,j] = mat
        # ...

        # ... replace tx/ty/tz by t1/t2/t3
        txs = [Symbol(tx) for tx in ['tx', 'ty', 'tz'][:dim]]
        for ti, tx in zip(tis, txs):
            kernel_expr = kernel_expr.subs(tx, ti)
        # ...

        # ...
        prelude = []
        for l,arr_ti in zip(lengths, arr_tis):
            prelude += [Assign(l, Len(arr_ti))]
        # ...

        # ...
        slices = [Slice(None,None)]*dim
        for i_row in range(0, n_rows):
            for i_col in range(0, n_cols):
                symbol = d_symbols[i_row,i_col]
                prelude += [Assign(symbol[slices], 0.)]
        # ...

        # ...
        body = []
        for i in range(dim-1,-1,-1):
            x = indices[i]
            rx = ranges[i]

            ti = tis[i]
            arr_ti = arr_tis[i]
            body += [Assign(ti, arr_ti[x])]

            for i_row in range(0, n_rows):
                for i_col in range(0, n_cols):
                    symbol = d_symbols[i_row,i_col]
                    symbol = symbol[x]
                    # TODO matrix case
                    body += [Assign(symbol, kernel_expr)]

            body = [For(x, rx, body)]
        # ...

        # ...
        body = prelude + body
#        print(body)
        # ...

        # ...
        self._basic_args = [*arr_tis]
        self._basic_args = tuple(self._basic_args)
        # ...

        # ...
        mats = []
        for i in range(0, n_rows):
            for j in range(0, n_cols):
                mats.append(d_symbols[i,j])
        mats = tuple(mats)
        # ...

        # function args
        func_args = self.build_arguments(mats)

        return FunctionDef(self.name, list(func_args), [], body)
