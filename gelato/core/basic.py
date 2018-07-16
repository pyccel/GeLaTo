# coding: utf-8
# TODO: - Unknown is not used here (mlhipy) remove it?

from numpy import unique

from sympy.core import Basic
from sympy.tensor import Indexed, IndexedBase
from sympy.core import Symbol
from sympy.core import Expr
from sympy.core.containers import Tuple
from sympy import Function
from sympy import Integer, Float

# ...
class Constant(Symbol):
    """
    Represents a constant symbol.

    Examples

    """
    pass
# ...

# ...
class Field(Symbol):
    """
    Represents a Field variable.

    Examples

    """
    pass
# ...

# ...
class Unknown(Symbol):
    """
    Represents an unknown function

    Examples

    """
    pass
# ...

class CalculusFunction(Function):
    """this class is needed to distinguish between functions and calculus
    functions when manipulating our expressions"""
    pass

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



_coeffs_registery = (Integer, Float, Constant)
