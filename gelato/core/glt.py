# -*- coding: utf-8 -*-
#
#
"""This module contains different functions to create and treate the GLT symbols."""

from sympy import Symbol
from sympy import Function
from sympy import bspline_basis
from sympy import lambdify
from sympy import cos
from sympy import sin
from sympy import Rational
from sympy import diff
from sympy import I as sympy_I
from sympy.core import Basic
from sympy.core.singleton import S
from sympy.simplify.simplify import nsimplify
from sympy import Tuple

from itertools import product


# TODO add it to glt_function
TOLERANCE    = 1.e-10

class BasicGlt(Function):
    """

    Examples

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

    @property
    def name(self):
        return self._name

    def _sympystr(self, printer):
        sstr = printer.doprint

        name = sstr(self.name)
        p = sstr(self.args[0])
        t = sstr(self.args[1])

        return '{name}({p},{t})'.format(name=name, p=p, t=t)
#        return '{name}_{p}({t})'.format(name=name, p=p, t=t)

# ...
class Mass(BasicGlt):
    """
    A class for the mass symbol
    """
    nargs = 2
    _name = 'Mass'

    @classmethod
    def eval(cls, p, t):

        if p is S.Infinity:
            raise NotImplementedError('Add symbol limit for p -> oo')

        elif isinstance(p, Symbol):
            return Mass(p, t, evaluate=False)

        elif isinstance(p, int):

            # ...
            r  = Symbol('r')

            pp = 2*p + 1
            N = pp + 1
            L = list(range(0, N + pp + 1))

            b0 = bspline_basis(pp, L, 0, r)
            bsp = lambdify(r, b0)
            # ...

            # ... we use nsimplify to get the rational number
            phi = []
            for i in range(0, p+1):
                y = bsp(p+1-i)
                y = nsimplify(y, tolerance=TOLERANCE, rational=True)
                phi.append(y)
            # ...

            # ...
            m = phi[0] * cos(S.Zero)
            for i in range(1, p+1):
                m += 2 * phi[i] * cos(i * t)
            # ...

            return m
# ...

# ...
class Stiffness(BasicGlt):
    """
    A class for the stiffness symbol
    """
    nargs = 2
    _name = 'Stiffness'

    @classmethod
    def eval(cls, p, t):

        if p is S.Infinity:
            raise NotImplementedError('Add symbol limit for p -> oo')

        elif isinstance(p, Symbol):
            return Stiffness(p, t, evaluate=False)

        elif isinstance(p, int):

            # ...
            r  = Symbol('r')

            pp = 2*p + 1
            N = pp + 1
            L = list(range(0, N + pp + 1))

            b0    = bspline_basis(pp, L, 0, r)
            b0_r  = diff(b0, r)
            b0_rr = diff(b0_r, r)
            bsp   = lambdify(r, b0_rr)
            # ...

            # ... we use nsimplify to get the rational number
            phi = []
            for i in range(0, p+1):
                y = bsp(p+1-i)
                y = nsimplify(y, tolerance=TOLERANCE, rational=True)
                phi.append(y)
            # ...

            # ...
            m = -phi[0] * cos(S.Zero)
            for i in range(1, p+1):
                m += -2 * phi[i] * cos(i * t)
            # ...

            return m
# ...

# ...
class Advection(BasicGlt):
    """
    A class for the advection symbol
    """
    nargs = 2
    _name = 'Advection'

    @classmethod
    def eval(cls, p, t):

        if p is S.Infinity:
            raise NotImplementedError('Add symbol limit for p -> oo')

        elif isinstance(p, Symbol):
            return Advection(p, t, evaluate=False)

        elif isinstance(p, int):

            # ...
            r  = Symbol('r')

            pp = 2*p + 1
            N = pp + 1
            L = list(range(0, N + pp + 1))

            b0   = bspline_basis(pp, L, 0, r)
            b0_r = diff(b0, r)
            bsp  = lambdify(r, b0_r)
            # ...

            # ... we use nsimplify to get the rational number
            phi = []
            for i in range(0, p+1):
                y = bsp(p+1-i)
                y = nsimplify(y, tolerance=TOLERANCE, rational=True)
                phi.append(y)
            # ...

            # ...
            m = -phi[0] * cos(S.Zero)
            for i in range(1, p+1):
                m += -2 * phi[i] * sin(i * t)
            # ...

            return m
# ...

# ...
class Bilaplacian(BasicGlt):
    """
    A class for the bilaplacian symbol
    """
    nargs = 2
    _name = 'Bilaplacian'

    @classmethod
    def eval(cls, p, t):

        if p is S.Infinity:
            raise NotImplementedError('Add symbol limit for p -> oo')

        elif isinstance(p, Symbol):
            return Bilaplacian(p, t, evaluate=False)

        elif isinstance(p, int):

            # ...
            r  = Symbol('r')

            pp = 2*p + 1
            N = pp + 1
            L = list(range(0, N + pp + 1))

            b0    = bspline_basis(pp, L, 0, r)
            b0_r  = diff(b0, r)
            b0_rr = diff(b0_r, r)
            b0_rrr = diff(b0_rr, r)
            b0_rrrr = diff(b0_rrr, r)
            bsp   = lambdify(r, b0_rrrr)
            # ...

            # ... we use nsimplify to get the rational number
            phi = []
            for i in range(0, p+1):
                y = bsp(p+1-i)
                y = nsimplify(y, tolerance=TOLERANCE, rational=True)
                phi.append(y)
            # ...

            # ...
            m = phi[0] * cos(S.Zero)
            for i in range(1, p+1):
                m += 2 * phi[i] * cos(i * t)
            # ...

            return m
# ...
