# coding: utf-8

from numpy import unique

from sympy.core import Basic
from sympy.core import Symbol
from sympy.core.containers import Tuple


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
    def __new__(cls, name, ldim=None):
        obj = Basic.__new__(cls, name)
        obj._ldim = ldim
        return obj

    @property
    def name(self):
        return self._args[0]

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
    def __new__(cls, space, name=None):
        return Basic.__new__(cls, space, name)

    @property
    def space(self):
        return self._args[0]

    @property
    def name(self):
        return self._args[1]

class TrialFunction(Symbol):
    """
    Represents a trial function as an element of a fem space.

    Examples

    >>> from gelato.fem.core import SplineFemSpace
    >>> from gelato.fem.core import TrialFunction
    >>> V = SplineFemSpace('V')
    >>> u = TrialFunction(V, 'u')
    """
    def __new__(cls, space, name=None):
        return Basic.__new__(cls, space, name)

    @property
    def space(self):
        return self._args[0]

    @property
    def name(self):
        return self._args[1]

