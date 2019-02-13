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

from .glt import (Mass,
                  Stiffness,
                  Advection,
                  Bilaplacian)


# ...
def _gelatize(a, degrees=None, evaluate=False, verbose=False):
    if isinstance(a, BilinearForm) and not(isinstance(a, BilinearAtomicForm)):
        expr = tensorize(a)
        if verbose:
            print('> tensorized = ', expr)
    else:
        expr = a

    if isinstance(expr, Add):
        args = [_gelatize(i, degrees=degrees, evaluate=evaluate) for i in expr.args]
        return Add(*args)

    elif isinstance(expr, Mul):
        coeffs  = [i for i in expr.args if isinstance(i, _coeffs_registery)]
        vectors = [i for i in expr.args if not(i in coeffs)]

        i = S.One
        if coeffs:
            i = Mul(*coeffs)

        j = S.One
        if vectors:
            args = [_gelatize(i, degrees=degrees, evaluate=evaluate) for i in vectors]
            j = Mul(*args)

        return Mul(i, j)

    if isinstance(expr, TensorProduct):
        args = [_gelatize(i, degrees=degrees, evaluate=evaluate) for i in expr.args]
        return Mul(*args)

    elif isinstance(expr, (Matrix, ImmutableDenseMatrix)):

        n_rows, n_cols = expr.shape
        lines = []
        for i in range(0, n_rows):
            line = []
            for j in range(0, n_cols):
                eij = _gelatize(expr[i,j], degrees=degrees, evaluate=evaluate)
                line.append(eij)
            lines.append(line)
        return Matrix(lines)

    elif isinstance(expr, BilinearAtomicForm):

        coord = expr.trial_spaces[0].coordinates

        # ... construct the fourier variable and the number of elements
        t_name = 't{}'.format(coord)
        n_name = 'n{}'.format(coord)

        t = Symbol(t_name)
        n = Symbol(n_name, integer=True)
        # ...

        # ...
        _coordinates = ['x', 'y', 'z']
        index = _coordinates.index(str(coord.name))
        # ...

        if evaluate and ( degrees is None ):
            raise ValueError('> degrees must be provided')

        # ... get the degree
        if not( degrees is None ):
            if not isinstance(degrees, (tuple, list, Tuple)):
                degrees = [degrees]

            p = degrees[index]
        else:
            p_name = 'p{}'.format(coord)
            p = Symbol(p_name, integer=True)
        # ...

        if isinstance(expr, MassForm):
            symbol = Mass(p, t, evaluate=evaluate)
            return symbol / n

        elif isinstance(expr, StiffnessForm):
            symbol = Stiffness(p, t, evaluate=evaluate)
            return symbol * n

        elif isinstance(expr, AdvectionForm):
            symbol = sympy_I * Advection(p, t, evaluate=evaluate)
            return symbol

        elif isinstance(expr, AdvectionTForm):
            symbol = - sympy_I * Advection(p, t, evaluate=evaluate)
            return symbol

        elif isinstance(expr, BilaplacianForm):
            symbol = Bilaplacian(p, t, evaluate=evaluate)
            return symbol * n**3

        else:
            raise NotImplementedError('TODO')

    return expr
# ...


def gelatize(a, degrees=None, n_elements=None, evaluate=False):

    if not isinstance(a, BilinearForm):
        raise TypeError('> Expecting a BilinearForm')

    # ... compute tensor form
    expr = TensorExpr(a)
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

    return expr
