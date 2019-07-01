# -*- coding: utf-8 -*-
#

from sympy.core import Symbol
from sympy.printing.latex import LatexPrinter as LatexPrinterSympy

class LatexPrinter(LatexPrinterSympy):

    def __init__(self, settings=None):
        self._enable_fourier_args = settings.pop('enable_fourier_args', True)

        LatexPrinterSympy.__init__(self, settings=settings)

    def _print_BasicGlt(self, name, p, t):

        fourier = ''
        if self._enable_fourier_args:
            t = "{}".format(t)
            fourier = "(" + t + ")"

        p = "{}".format(p)

        name = r"\mathfrak{" + name + "}"
        return name + '_{' + p + "}" + fourier

    def _print_Mass(self, expr, **kwargs):
        return self._print_BasicGlt('m', *expr.args)

    def _print_Stiffness(self, expr, **kwargs):
        return self._print_BasicGlt('s', *expr.args)

    def _print_Advection(self, expr, **kwargs):
        return self._print_BasicGlt('a', *expr.args)

    def _print_Bilaplacian(self, expr, **kwargs):
        return self._print_BasicGlt('b', *expr.args)

def latex(expr, **settings):

    coords = ['x', 'y', 'z']
    names = {}
    for x in coords:
        p = Symbol('p{}'.format(x), integer=True)
        n = Symbol('n{}'.format(x), integer=True)
        t = Symbol('t{}'.format(x))

        names[p] = Symbol('p_{}'.format(x))
        names[t] = Symbol(r'\theta_{}'.format(x))
        names[n] = Symbol('n_{}'.format(x))

    expr = expr.subs(names)

    return LatexPrinter(settings).doprint(expr)
