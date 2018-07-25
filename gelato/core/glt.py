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


# ...
class glt_symbol_m(Function):
    """
    A class for the mass symbol
    """
    nargs = 3

    @classmethod
    def eval(cls, p, t, n=None):

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

        # ... scaling
        if not( n is None ):
            if isinstance(n, Symbol):
                m = m/n

            else:
                m *= Rational(1,n)
        # ...

        return m
# ...

# ...
class glt_symbol_s(Function):
    """
    A class for the stiffness symbol
    """
    nargs = 3

    @classmethod
    def eval(cls, p, t, n=None):

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

        # ... scaling
        if not( n is None ):
            m *= n
        # ...

        return m
# ...

# ...
class glt_symbol_a(Function):
    """
    A class for the advection symbol
    """
    nargs = 3

    @classmethod
    def eval(cls, p, t, n=None):

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
class glt_symbol_b(Function):
    """
    A class for the bilaplacian symbol
    """
    nargs = 3

    @classmethod
    def eval(cls, p, t, n=None):

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

        # ... scaling
        if not( n is None ):
            m *= n**3
        # ...

        return m
# ...
