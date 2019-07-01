# -*- coding: utf-8 -*-
#

from sympy.core import Symbol
from sympy import lambdify

from numpy import linspace, pi
from matplotlib import pyplot as plt

from .glt import Stiffness


def plot_stiffness_symbols(degrees=[2, 3], nx=100):
    """
    Plots the stiffness symbol for different degrees on a given number of elements.

    degrees: list, tuple
        a list/tuple of spline degrees

    nx: int
        number of points used for plot
    """

    colors = ["r", "b", "g", "y", "m", "k", "c"]

    for i,p in enumerate(degrees):

        t = Symbol('t')

        symbol = Stiffness(p, t)

        # ... make the symbol a numeric function, that can be evaluated
        f = lambdify(t, symbol, "numpy")

        t1 = linspace(-pi,pi, nx)
        w = f(t1)
        # ...

        plt.plot(t1, w.real, "-"+colors[i], label="$p=" + str(p) + "$")

    plt.xlabel('fequencies')
    plt.ylabel('$\mathfrak{s}_p$')
    plt.legend(loc=9)
