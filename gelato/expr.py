# coding: utf-8

from sympy import Function
from sympy.core import Add, Mul
from sympy import I as sympy_I
from sympy import Symbol
from sympy.core.containers import Tuple
from sympy import S
from sympy.core import Expr, Basic, AtomicExpr
from sympy import simplify
from sympy import Matrix, ImmutableDenseMatrix
from sympy.physics.quantum import TensorProduct

from sympde.expr import LinearForm, BilinearForm
from sympde.expr import TensorExpr

from sympde.expr import Mass as MassForm
from sympde.expr import Stiffness as StiffnessForm
from sympde.expr import Advection as AdvectionForm
from sympde.expr import AdvectionT as AdvectionTForm
from sympde.expr import Bilaplacian as BilaplacianForm
from sympde.expr import Basic1dForm
from sympde.topology import SymbolicExpr

from .glt import (Mass,
                  Stiffness,
                  Advection,
                  Bilaplacian)


def gelatize(a, degrees=None, n_elements=None, evaluate=False, mapping=None,
             human=False):

    if not isinstance(a, BilinearForm):
        raise TypeError('> Expecting a BilinearForm')

    # ... compute tensor form
    expr = TensorExpr(a, mapping=mapping)
    # ...

    # ...
    atoms = list(expr.atoms(TensorExpr))
    if atoms:
        for i in atoms:
            expr = expr.subs(i, i._args[0])
    # ...

    # ... coordinates as strings
    coordinates = ['x', 'y', 'z']
    # ...

    # ... get the degree
    if not( degrees is None ):
        if not isinstance(degrees, (tuple, list, Tuple)):
            degrees = [degrees]

    else:
        degrees = [Symbol('p{}'.format(i), integer=True) for i in coordinates]
    # ...

    # ... coordinates as symbols
    coordinates = [Symbol(i) for i in coordinates]
    # ...

    # ...
    forms = list(expr.atoms(Basic1dForm))
    for form in forms:

        p = degrees[form.axis]
        coord = coordinates[form.axis]

        # ... construct the fourier variable and the number of elements
        t_name = 't{}'.format(coord)
        n_name = 'n{}'.format(coord)

        t = Symbol(t_name)
        n = Symbol(n_name, integer=True)
        # ...

        if isinstance(form, MassForm):
            symbol = Mass(p, t, evaluate=evaluate)
            expr = expr.subs(form, symbol / n)

        elif isinstance(form, StiffnessForm):
            symbol = Stiffness(p, t, evaluate=evaluate)
            expr = expr.subs(form, symbol * n)

        elif isinstance(form, AdvectionForm):
            symbol = sympy_I * Advection(p, t, evaluate=evaluate)
            expr = expr.subs(form, symbol)

        elif isinstance(form, AdvectionTForm):
            symbol = - sympy_I * Advection(p, t, evaluate=evaluate)
            expr = expr.subs(form, symbol)

        elif isinstance(form, BilaplacianForm):
            symbol = Bilaplacian(p, t, evaluate=evaluate)
            expr = expr.subs(form, symbol * n**3)

        else:
            raise NotImplementedError('{} not available yet'.format(type(form)))
    # ...

    # ...
    if not( n_elements is None ):
        dim = a.ldim
        ns = ['nx', 'ny', 'nz'][:dim]
        ns = [Symbol(i, integer=True) for i in ns]

        if isinstance(n_elements, int):
            n_elements = [n_elements]*dim

        if not( len(ns) == len(n_elements) ):
            raise ValueError('Wrong size for n_elements')

        for n,v in zip(ns, n_elements):
            expr = expr.subs(n, v)
    # ...

    # ...
    if mapping and human:
        expr = SymbolicExpr(expr)
    # ...

    return expr
