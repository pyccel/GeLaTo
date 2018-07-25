# -*- coding: utf-8 -*-
#

from sympy.core import Symbol
from sympy import lambdify

from numpy import linspace, pi
from matplotlib import pyplot as plt

from .glt import (glt_symbol_m,
                  glt_symbol_s,
                  glt_symbol_a,
                  glt_symbol_b)


def plot_stiffness_symbols(n, degrees=[2, 3], nx=100):
    """
    Plots the stiffness symbol for different degrees on a given number of elements.

    n: int
        number of elements of the grid

    degrees: list, tuple
        a list/tuple of spline degrees

    nx: int
        number of points used for plot
    """

    colors = ["r", "b", "g", "y", "m", "k", "c"]

    for i,p in enumerate(degrees):

        t = Symbol('t')

        symbol = glt_symbol_s(p, t, n)

        # ... make the symbol a numeric function, that can be evaluated
        f = lambdify(t, symbol, "numpy")

        t1 = linspace(-pi,pi, nx)
        w = f(t1)
        # ...

        plt.plot(t1, w.real, "-"+colors[i], label="$p=" + str(p) + "$")
        plt.legend(loc=9)
    plt.show()
