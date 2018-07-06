# coding: utf-8

from numpy import unique

from sympy.core import Basic
from sympy.tensor import Indexed, IndexedBase
from sympy.core.compatibility import is_sequence
from sympy.core import Symbol
from sympy.core import Expr
from sympy.core.containers import Tuple

from gelato.calculus import Field, Constant
from gelato.calculus import _partial_derivatives
from gelato.calculus import get_index_derivatives, get_atom_derivatives


class DottedName(Basic):
    """
    Represents a dotted variable.

    Examples

    >>> from pyccel.ast.core import DottedName
    >>> DottedName('matrix', 'n_rows')
    matrix.n_rows
    >>> DottedName('pyccel', 'stdlib', 'parallel')
    pyccel.stdlib.parallel
    """
    def __new__(cls, *args):
        return Basic.__new__(cls, *args)

    @property
    def name(self):
        return self._args

    def __str__(self):
        return '.'.join(str(n) for n in self.name)

    def _sympystr(self, printer):
        sstr = printer.doprint
        return '.'.join(sstr(n) for n in self.name)

# ...
class FemSpace(Basic):
    """
    Represents a symbolic finite elements space.

    Examples

    """
    _ldim = None
    _shape = None
    _is_vector = False
    _is_block = False
    def __new__(cls, name, ldim=None, shape=None, is_vector=False, is_block=False):
        if is_vector or is_block:
            if shape is None:
                raise ValueError('shape must be provided for a vector/block space')

        obj = Basic.__new__(cls, name)
        obj._ldim = ldim
        obj._shape = shape
        obj._is_vector = is_vector
        obj._is_block = is_block
        return obj

    @property
    def name(self):
        return self._args[0]

    @property
    def ldim(self):
        return self._ldim

    @property
    def shape(self):
        return self._shape

    @property
    def is_vector(self):
        return self._is_vector

    @property
    def is_block(self):
        return self._is_block

    @property
    def ldim(self):
        return self._ldim

    def _sympystr(self, printer):
        sstr = printer.doprint
        return sstr(self.name)
# ...

class SplineFemSpace(FemSpace):
    """
    Represents a 1D fem space of splines.

    Examples

    >>> from gelato.fem.core import SplineFemSpace
    >>> V = SplineFemSpace('V')
    >>> V.degree
    V.p
    >>> V.n_elements
    V.n
    """
    _ldim = 1
    _degree = None
    _n_elements = None
    def __new__(cls, name):
        return Basic.__new__(cls, name)

    @property
    def name(self):
        return self._args[0]

    @property
    def degree(self):
        if not (self._degree is None):
            return self._degree
        else:
            return DottedName(self.name, 'p')

    @property
    def n_elements(self):
        if not (self._n_elements is None):
            return self._n_elements
        else:
            return DottedName(self.name, 'n')

    def _sympystr(self, printer):
        sstr = printer.doprint
        return sstr(self.name)

class TensorFemSpace(FemSpace):
    """
    Represents a tensor product of fem spaces.

    Examples

    >>> from gelato.fem.core import SplineFemSpace
    >>> from gelato.fem.core import TensorFemSpace
    >>> V1 = SplineFemSpace('V1')
    >>> V2 = SplineFemSpace('V2')
    >>> V  = TensorFemSpace('V', V1, V2)
    >>> V.degree
    (V1.p, V2.p)
    >>> V.n_elements
    (V1.n, V2.n)
    """
    _ldim = 0

    def __new__(cls, name, *args):
        return Basic.__new__(cls, name, *args)

    @property
    def name(self):
        return self._args[0]

    @property
    def spaces(self):
        return self._args[1:]

    @property
    def degree(self):
        args = [V.degree for V in self.spaces]
        return Tuple(*args)

    @property
    def n_elements(self):
        args = [V.n_elements for V in self.spaces]
        return Tuple(*args)

    @property
    def ldim( self ):
        return sum([V.ldim for V in self.spaces])

    def _sympystr(self, printer):
        sstr = printer.doprint
        return sstr(self.name)


class VectorFemSpace(FemSpace):
    """
    Represents a tensor product of fem spaces.

    Examples

    >>> from gelato.fem.core import SplineFemSpace
    >>> from gelato.fem.core import TensorFemSpace
    >>> from gelato.fem.core import VectorFemSpace
    >>> V1 = SplineFemSpace('V1')
    >>> V2 = SplineFemSpace('V2')
    >>> Vx = TensorFemSpace('Vx', V1, V2)
    >>> Vy = TensorFemSpace('Vy', V2, V1)
    >>> V = VectorFemSpace('V', Vx, Vy)
    >>> V.degree
    ((V1.p, V2.p), (V2.p, V1.p))
    >>> V.n_elements
    ((V1.n, V2.n), (V2.n, V1.n))
    """
    _ldim = 0

    def __new__(cls, name, *args):
        return Basic.__new__(cls, name, *args)

    @property
    def name(self):
        return self._args[0]

    @property
    def spaces(self):
        return self._args[1:]

    @property
    def degree(self):
        args = [V.degree for V in self.spaces]
        return Tuple(*args)

    @property
    def n_elements(self):
        args = [V.n_elements for V in self.spaces]
        return Tuple(*args)

    @property
    def ldim( self ):
        # make sure that all spaces have the same parametric dimension
        ldims = [V.ldim for V in self.spaces]
        assert (len(unique(ldims)) == 1)

        return ldims[0]

    def _sympystr(self, printer):
        sstr = printer.doprint
        return sstr(self.name)

class TestFunction(Symbol):
    """
    Represents a test function as an element of a fem space.

    Examples

    >>> from gelato.fem.core import SplineFemSpace
    >>> from gelato.fem.core import TestFunction
    >>> V = SplineFemSpace('V')
    >>> phi = TestFunction(V, 'phi')
    """
    _space = None
    def __new__(cls, space, name=None):
        obj =  Basic.__new__(cls, name)
        obj._space = space
        return obj

    @property
    def space(self):
        return self._space

    @property
    def name(self):
        return self._args[0]

    def _sympystr(self, printer):
        sstr = printer.doprint
        return sstr(self.name)


# this class is needed, otherwise sympy will convert VectorTestFunction to
# IndexedBase
class IndexedTestTrial(Indexed):
    """Represents a mathematical object with indices.

    """
    is_commutative = True
    is_Indexed = True
    is_symbol = True
    is_Atom = True

    def __new__(cls, base, *args, **kw_args):
        assert(isinstance(base, VectorTestFunction))

        if not args:
            raise IndexException("Indexed needs at least one index.")

        return Expr.__new__(cls, base, *args, **kw_args)

    # free_symbols is redefined otherwise an expression u[0].free_symbols will
    # give the error:  AttributeError: 'int' object has no attribute 'free_symbols'
    @property
    def free_symbols(self):
        base_free_symbols = self.base.free_symbols
        symbolic_indices = [i for i in self.indices if isinstance(i, Basic)]
        if len(symbolic_indices) > 0:
            raise ValueError('symbolic indices not yet available')

        return base_free_symbols

        # TODO uncomment if needed
#        indices_free_symbols = {
#            fs for i in symbolic_indices for fs in i.free_symbols}
#        if base_free_symbols:
#            return {self} | base_free_symbols | indices_free_symbols
#        else:
#            return indices_free_symbols


class VectorTestFunction(Symbol, IndexedBase):
    """
    Represents a vector test function as an element of a fem space.

    Examples

    """
    _space = None
    def __new__(cls, space, name=None):
        if not(space.is_vector) and not(space.is_block):
            raise ValueError('Expecting a vector/block space')

        obj = Basic.__new__(cls, name)
        obj._space = space
        return obj

    @property
    def space(self):
        return self._space

    @property
    def name(self):
        return self._args[0]

    @property
    def shape(self):
        # we return a list to make it compatible with IndexedBase sympy object
        return [self.space.shape]

    def __getitem__(self, *args):

        if self.shape and len(self.shape) != len(args):
            raise IndexException("Rank mismatch.")

        if not(len(args) == 1):
            raise ValueError('expecting exactly one argument')

        assumptions ={}
        obj = IndexedTestTrial(self, *args)
        return obj

