# coding: utf-8

from sympy.core import Basic
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

class SplineFemSpace(Basic):
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

class TensorFemSpace(Basic):
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

    def _sympystr(self, printer):
        sstr = printer.doprint
        return sstr(self.name)

class VectorFemSpace(Basic):
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

    def _sympystr(self, printer):
        sstr = printer.doprint
        return sstr(self.name)
