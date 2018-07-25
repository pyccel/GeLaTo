# -*- coding: utf-8 -*-
#
#
# TODO use to_assign and post processing as expression and not latex => helpful
#      for Fortran and Lua (code gen).
"""This module contains different functions to create and treate the GLT symbols."""

import numpy as np

from sympy.core.sympify import sympify
from sympy.simplify.simplify import simplify
from sympy import Symbol
from sympy import Lambda
from sympy import Function
from sympy import bspline_basis
from sympy import lambdify
from sympy import cos
from sympy import sin
from sympy import Rational
from sympy import diff
from sympy import Matrix
from sympy import latex
from sympy import Integral
from sympy import I as sympy_I
from sympy.core import Basic
from sympy.core.singleton import S
from sympy.simplify.simplify import nsimplify
from sympy.utilities.lambdify import implemented_function
from sympy.matrices.dense import MutableDenseMatrix
from sympy import Mul, Add
from sympy import Tuple
from sympy import postorder_traversal
from sympy import preorder_traversal
from sympy import Indexed
from sympy import IndexedBase
from sympy import Lambda

from itertools import product

# ...
class glt_function(Function):
    """

    Examples
    ========

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

    @classmethod
    def eval(cls, *_args):
        """."""

        if not _args:
            return

        f = _args[0]
        n = _args[1]
        p = _args[2]

        if isinstance(n, (Tuple, list, tuple)):
            dim = len(n)
        else:
            dim = 1
            n = [n]
            p = [p]
        discretization = {"n_elements": n, "degrees": p}

        F = glt_symbol(f, dim=dim,
                       discretization=discretization,
                       evaluate=True,
                       verbose=False)

        return F
# ...

# ...
class glt_symbol_m(Function):
    """
    A class for the mass symbol
    """
    nargs = 3

    @classmethod
    def eval(cls, n, p, t):
        # ...
        if not 0 <= p:
            raise ValueError("must have 0 <= p")
        if not 0 <= n:
            raise ValueError("must have 0 <= n")
        # ...

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
    def eval(cls, n, p, t):
        # ...
        if not 0 <= p:
            raise ValueError("must have 0 <= p")
        if not 0 <= n:
            raise ValueError("must have 0 <= n")
        # ...

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
    def eval(cls, n, p, t):
        # ...
        if not 0 <= p:
            raise ValueError("must have 0 <= p")
        if not 0 <= n:
            raise ValueError("must have 0 <= n")
        # ...

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

#        # ... make it pure imaginary
#        m *= sympy_I
#        # ...

        return m
# ...

# ...
class glt_symbol_b(Function):
    """
    A class for the bilaplacian symbol
    """
    nargs = 3

    @classmethod
    def eval(cls, n, p, t):
        # ...
        if not 0 <= p:
            raise ValueError("must have 0 <= p")
        if not 0 <= n:
            raise ValueError("must have 0 <= n")
        # ...

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
        m *= n**3
        # ...

        return m
# ...

# ...
def glt_symbol_laplace(discretization, \
                       verbose=False, evaluate=True, \
                       instructions=[], \
                       **settings):
    """
    Returns the Laplace symbol for a given discretization.

    discretization: dict
        a dictionary that contains the used discretization

    verbose: bool
        talk more

    evaluate: bool
        causes the evaluation of the atomic symbols, if true

    instructions: list
        a list to keep track of the applied instructions.

    settings: dict
        dictionary for different settings
    """
    raise NotImplementedError('')
# ...

# ...
def glt_integrate(expr, domain="Omega"):
    """
    Adds the integral to the expression. needed for printing.

    domain: str
        domain over which we integrate the expression.
    """
    return Integral(expr, domain)
# ...

# ...
def glt_formatting(expr, **settings):
    """
    Formatting the glt symbol, prior to calling a printer

    expr: sympy.Expression
        a sympy expression

    settings: dict
        dictionary for different settings
    """

    # ...
    try:
        domain = settings["glt_integrate"]
        domain = sympify(str(domain))

        expr = glt_integrate(expr, domain)
    except:
        pass
    # ...

    # ...
    try:
        fmt = settings["glt_formatting_atoms"]
        if fmt:
            expr = glt_formatting_atoms(expr, **settings)
    except:
        pass
    # ...

    return expr
# ...

# ...
def glt_formatting_atoms(expr, **settings):
    """
    Formatting the glt symbol atoms, prior to calling a printer

    expr: sympy.Expression
        a sympy expression

    settings: dict
        dictionary for different settings
    """
    # TODO do we still need it?

    # ...
    dim    = 3
    prefix = "\mathfrak{"
    suffix = "}"
    # ...

    # ...
    for k in range(0, dim):
        # ...
        t = Symbol('t'+str(k+1))

        for s in ["m", "s", "a", "b"]:
            sk = s + str(k+1)
            s_old = Symbol(sk)
            s_new = Symbol(prefix + sk + suffix)

            expr = expr.subs({s_old: s_new})
        # ...
    # ...

    return expr
# ...



# ...
def latex_title_as_paragraph(title):
    """
    Returns the title as a paragraph.

    title: str
        a string for the paragraph title
    """
    return "\paragraph{" + str(title) + "}"
# ...

# ...
def glt_latex_definitions():
    """
    Returns the definitions of the atomic symbols for the GLT.
    """
    # ...
    t = Symbol('t')
    m = Symbol('m')
    s = Symbol('s')
    a = Symbol('a')
    b = Symbol('b')
    # ...

    # ...
    def formula(symbol):
        """
        returns the latex formula for the mass symbol.
        """
        txt_m = r"\phi_{2p+1}(p+1) + 2 \sum_{k=1}^p \phi_{2p+1}(p+1-k) \cos(k \theta)"
        txt_s = r"- {\phi}''_{2p+1}(p+1) - 2 \sum_{k=1}^p {\phi}''_{2p+1}(p+1-k) \cos(k \theta)"
        txt_a = r"\phi_{2p+1}(p+1) + 2 \sum_{k=1}^p \phi_{2p+1}(p+1-k) \cos(k \theta)"
        txt_b = r"{\phi}''''_{2p+1}(p+1) - 2 \sum_{k=1}^p {\phi}''''_{2p+1}(p+1-k) \cos(k \theta)"

        if str(symbol) == "m":
            return txt_m
        elif str(symbol) == "s":
            return txt_s
        elif str(symbol) == "a":
            return txt_a
        elif str(symbol) == "b":
            return txt_b
        else:
            print ("not yet available.")
    # ...

    # ...
    definitions = {r"m(\theta)": formula(m), \
                   r"s(\theta)": formula(s), \
                   r"a(\theta)": formula(a)}
    # ...

    return definitions
# ...

# ...
def glt_latex_names():
    """
    returns latex names for basis and atoms
    """
    # ...
    dim = 3

    symbol_names = {}
    # ...

    # ... rename basis
    B = "N"
    for i in ["i","j"]:
        Bi = B + i
        symbol_names[Symbol(Bi)] = B + "_" + i
    # ...

    # ... rename basis derivatives in the logical domain
    args_x = ARGS_x[:dim]
    args_u = ARGS_u[:dim]
    B = "N"
    for i in ["i","j"]:
        Bi = B + i
        for u in args_u + args_x:
            Bi_u = Bi + "_" + u
            partial = "\partial_" + u
            symbol_names[Symbol(Bi_u)] = partial + B + "_" + i
    # ...

    # ... rename the tensor basis derivatives
    B = "N"
    for i in ["i","j"]:
        Bi = B + i
        for k in range(0, dim):
            for s in ["", "s", "ss", "sss", "ssss"]:
                Bk  = Bi + str(k+1)
                _Bk = B + "_{" + i + "_" + str(k+1) + "}"

                if len(s) > 0:
                    prime = len(s) * "\prime"

                    Bk += "_" + s
                    _Bk = B + "^{" + prime + "}" \
                            + "_{" + i + "_" + str(k+1) + "}"

                symbol_names[Symbol(Bk)] = _Bk
    # ...

    # ... TODO add flag to choose which kind of printing:
#    for k in range(0, dim):
#        # ...
#        symbol_names[Symbol('m'+str(k+1))] = "\mathfrak{m}_" + str(k+1)
#        symbol_names[Symbol('s'+str(k+1))] = "\mathfrak{s}_" + str(k+1)
#        symbol_names[Symbol('a'+str(k+1))] = "\mathfrak{a}_" + str(k+1)
#        # ...

    degree = "p"
    for k in range(0, dim):
        # ...
        for s in ["m", "s", "a", "b"]:
            symbol_names[Symbol(s+str(k+1))] = r"\mathfrak{" + s + "}_" \
                                             + degree \
                                             + r"(\theta_" \
                                             + str(k+1) + ")"
        # ...
    # ...

    # ...
    for k in range(0, dim):
        symbol_names[Symbol("t"+str(k+1))] = r"\theta_" + str(k+1)
    # ...

    return symbol_names
# ...

# ...
def get_sympy_printer_settings(settings):
    """
    constructs the dictionary for sympy settings needed for the printer.

    settings: dict
        dictionary for different settings
    """
    sets = {}
    for key, value in list(settings.items()):
        if key not in SETTINGS:
            sets[key] = value
    return sets
# ...

# ...
def glt_latex(expr, **settings):
    """
    returns the latex expression of expr.

    expr: sympy.Expression
        a sympy expression

    settings: dict
        dictionary for different settings
    """
    # ...
    if type(expr) == dict:
        d_expr = {}
        try:
            mode = settings["mode"]
        except:
            mode = "plain"

        sets = settings.copy()
        sets["mode"] = "plain"
        for key, txt in list(expr.items()):
            d_expr[key] = glt_latex(txt, **sets)

        return d_expr
    # ...

    # ...
    try:
        from gelato.expression import glt_formatting
        fmt = settings["glt_formatting"]
        if fmt:
            expr = glt_formatting(expr, **settings)
    except:
        pass
    # ...

    # ...
    try:
        smp = settings["glt_simplify"]
        if smp:
            expr = simplify(expr)
    except:
        pass
    # ...

    # ...
    sets = get_sympy_printer_settings(settings)
    # ...

    return latex(expr, symbol_names=glt_latex_names(), **sets)
# ...

# ...
def print_glt_latex(expr, **settings):
    """
    Prints the latex expression of expr.

    settings: dict
        dictionary for different settings
    """
    print((glt_latex(expr, **settings)))
# ...
