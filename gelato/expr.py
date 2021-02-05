# coding: utf-8

from sympy import I as sympy_I
from sympy import Symbol
from sympy.core.containers import Tuple
from sympy import S
from sympy.core import Expr, Basic
from sympy import simplify, expand

from sympde.expr import BilinearForm
from sympde.expr import TensorExpr

from sympde.expr import Mass as MassForm
from sympde.expr import Stiffness as StiffnessForm
from sympde.expr import Advection as AdvectionForm
from sympde.expr import AdvectionT as AdvectionTForm
from sympde.expr import Bilaplacian as BilaplacianForm
from sympde.expr import Basic1dForm
from sympde.topology import SymbolicExpr
from sympde.calculus.matrices import SymbolicDeterminant

from .glt import (BasicGlt, Mass, Stiffness, Advection, Bilaplacian)

__all__ = ('gelatize', 'GltExpr')

#==============================================================================
def gelatize(a, degrees=None, n_elements=None, evaluate=False, mapping=None,
             human=False, expand=False):

    if not isinstance(a, BilinearForm):
        raise TypeError('> Expecting a BilinearForm')

    dim = a.ldim

    # ... compute tensor form
    expr = TensorExpr(a, mapping=mapping, expand=expand)
    # ...

    # ...
    atoms = list(expr.atoms(TensorExpr))
    if atoms:
        for i in atoms:
            expr = expr.subs(i, i._args[0])
    # ...

    # ... coordinates as strings
    coordinates = ['x', 'y', 'z'][:dim]
    # ...

    # ... symbolic degree, ncells, fourier variables
    ps = [Symbol('p{}'.format(i), integer=True) for i in coordinates]
    ns = [Symbol('n{}'.format(i), integer=True) for i in coordinates]
    ts = [Symbol('t{}'.format(i))               for i in coordinates]
    # ...

    # ...
    forms = list(expr.atoms(Basic1dForm))
    for form in forms:

        p = ps[form.axis]
        n = ns[form.axis]
        t = ts[form.axis]

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

        if isinstance(n_elements, int):
            n_elements = [n_elements]*dim

        if not( len(ns) == len(n_elements) ):
            raise ValueError('Wrong size for n_elements')

        for n,v in zip(ns, n_elements):
            expr = expr.subs(n, v)
    # ...

    # ... get the degree
    if not( degrees is None ):
        if not isinstance(degrees, (tuple, list, Tuple)):
            degrees = [degrees]

        assert(len(degrees) == len(ps))

        atoms = list(expr.atoms(BasicGlt))

        d = {}
        for p,v in zip(ps, degrees):
            d[p] = v

        for atom in atoms:
            p,t = atom.args[:]
            newp = d[p]
            newatom = atom.func(newp, t)
            expr = expr.subs(atom, newatom)
    # ...

    # ...
    if mapping and human:
        expr *= SymbolicDeterminant(mapping)
        expr  = SymbolicExpr(expr)
    # ...

    return expr

#==============================================================================
# TODO add __call__
class GltExpr(Expr):
    is_Function = True

    def __new__(cls, form):

        assert(isinstance(form, BilinearForm))

        # ...
        dim = form.ldim

        fourier_vars = [Symbol(i) for i in ['tx', 'ty', 'tz'][:dim]]
#        space_vars   = [Symbol(i) for i in ['x', 'y', 'z'][:dim]]
        atoms  = form.atoms(Symbol)
        space_vars   = [i for i in atoms if i in form.coordinates]
        # ...

        return Basic.__new__(cls, fourier_vars, space_vars, form)

    @property
    def fourier_variables(self):
        return self._args[0]

    @property
    def space_variables(self):
        return self._args[1]

    @property
    def form(self):
        return self._args[2]

    @property
    def ldim(self):
        return self.form.ldim

    @property
    def coordinates(self):
        return self.form.coordinates

    @property
    def fields(self):
        return self.form.fields

    @property
    def constants(self):
        return self.form.constants

    def __call__(self, *args, **kwargs):

        mapping    = kwargs.pop('mapping',    None)
        human      = kwargs.pop('human',      True)
        degrees    = kwargs.pop('degrees',    None)
        n_elements = kwargs.pop('n_elements', None)

        expr =  gelatize( self.form,
                          degrees = degrees, n_elements = n_elements,
                          mapping = mapping, human = human, evaluate = True)

        dim = self.ldim

        for i in ['tx', 'ty', 'tz'][:dim] +  ['x', 'y', 'z'][:dim]:
            I = kwargs.pop(i, None)
            if not(I is None):
                S = Symbol(i)
                expr = expr.subs(S, I)

        return expr
