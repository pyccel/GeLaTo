# -*- coding: utf-8 -*-
#
#
# TODO: - add tabulation for Bilaplacian

"""This module contains different functions to create and treate the GLT symbols."""

from sympy import Symbol
from sympy import Function
from sympy import bspline_basis
from sympy import sympify
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
from sympy import Rational

from itertools import product

# ............................................
# tabular values
# ............................................
P_MAX = 8

# ... phi_{2p+1}
d_phi = {}
d_phi[1] = [Rational(2,3),
            Rational(1,6)]
d_phi[2] = [Rational(11,20),
            Rational(13,60),
            Rational(1,120)]
d_phi[3] = [Rational(151,315),
            Rational(397,1680),
            Rational(1,42),
            Rational(1,5040)]
d_phi[4] = [Rational(15619,36288),
            Rational(44117,181440),
            Rational(913,22680),
            Rational(251,181440),
            Rational(1,362880)]
d_phi[5] = [Rational(655177,1663200),
            Rational(1623019,6652800),
            Rational(1093,19800),
            Rational(50879,13305600),
            Rational(509,9979200),
            Rational(1,39916800)]
d_phi[6] = [Rational(27085381,74131200),
            Rational(125468459,518918400),
            Rational(28218769,415134720),
            Rational(910669,124540416),
            Rational(82207,345945600),
            Rational(1363,1037836800),
            Rational(1,6227020800)]
d_phi[7] = [Rational(2330931341,6810804000),
            Rational(103795866137,435891456000),
            Rational(6423562433,81729648000),
            Rational(15041229521,1307674368000),
            Rational(26502841,40864824000),
            Rational(13824739,1307674368000),
            Rational(2047,81729648000),
            Rational(1,1307674368000)]
d_phi[8] = [Rational(12157712239,37638881280),
            Rational(8313722318537,35568742809600),
            Rational(7763913237097,88921857024000),
            Rational(317627331799,19760412672000),
            Rational(23667665053,17784371404800),
            Rational(297507989,7113748561920),
            Rational(704339,1976041267200),
            Rational(851,2309658624000),
            Rational(1,355687428096000)]
# ...

# ... phi'_{2p+1}
d_phi_r = {}
d_phi_r[1] = [0,
              Rational(1,2)]
d_phi_r[2] = [0,
              Rational(5,12),
              Rational(1,24)]
d_phi_r[3] = [0,
              Rational(49,144),
              Rational(7,90),
              Rational(1,720)]
d_phi_r[4] = [0,
              Rational(809,2880),
              Rational(289,2880),
              Rational(41,6720),
              Rational(1,40320)]
d_phi_r[5] = [0,
              Rational(6787,28800),
              Rational(16973,151200),
              Rational(5203,403200),
              Rational(253,907200),
              Rational(1,3628800)]
d_phi_r[6] = [0,
              Rational(728741,3628800),
              Rational(1700933,14515200),
              Rational(441337,21772800),
              Rational(10777,10886400),
              Rational(2041,239500800),
              Rational(1,479001600)]
d_phi_r[7] = [0,
              Rational(35263201,203212800),
              Rational(4489301,38102400),
              Rational(5532241,203212800),
              Rational(233021,104781600),
              Rational(6323,121927680),
              Rational(31,165110400),
              Rational(1,87178291200)]
d_phi_r[8] = [0,
              Rational(11102502613,73156608000),
              Rational(8480306503,73156608000),
              Rational(7939969,238436352),
              Rational(3146582819,804722688000),
              Rational(353015251,2092278988800),
              Rational(775319,387459072000),
              Rational(32759,10461394944000),
              Rational(1,20922789888000)]
# ...

# ... phi''_{2p+1}
d_phi_rr = {}
d_phi_rr[1] = [-2,
               1]
d_phi_rr[2] = [-1,
               Rational(1,3),
               Rational(1,6)]
d_phi_rr[3] = [Rational(-2,3),
               Rational(1,8),
               Rational(1,5),
               Rational(1,120)]
d_phi_rr[4] = [Rational(-35,72),
               Rational(11,360),
               Rational(17,90),
               Rational(59,2520),
               Rational(1,5040)]
d_phi_rr[5] = [Rational(-809,2160),
               Rational(-1,64),
               Rational(31,189),
               Rational(907,24192),
               Rational(25,18144),
               Rational(1,362880)]
d_phi_rr[6] = [Rational(-4319,14400),
               Rational(-11731,302400),
               Rational(6647,48384),
               Rational(3455,72576),
               Rational(2251,604800),
               Rational(113,2217600),
               Rational(1,39916800)]
d_phi_rr[7] = [Rational(-56057,226800),
               Rational(-104159,2073600),
               Rational(43993,388800),
               Rational(333361,6220800),
               Rational(14623,2138400),
               Rational(16081,68428800),
               Rational(73,55598400),
               Rational(1,6227020800)]
d_phi_rr[8] = [Rational(-35263201,169344000),
               Rational(-253354477,4572288000),
               Rational(30188519,326592000),
               Rational(44897821,798336000),
               Rational(36700199,3592512000),
               Rational(58605299,93405312000),
               Rational(382201,36324288000),
               Rational(131,5230697472),
               Rational(1,1307674368000)]
# ...

## ...
#d_phi[1] = [Rational(,),
#            Rational(,)]
#d_phi[2] = [Rational(,),
#            Rational(,),
#            Rational(,)]
#d_phi[3] = [Rational(,),
#            Rational(,),
#            Rational(,),
#            Rational(,)]
#d_phi[4] = [Rational(,),
#            Rational(,),
#            Rational(,),
#            Rational(,),
#            Rational(,)]
#d_phi[5] = [Rational(,),
#            Rational(,),
#            Rational(,),
#            Rational(,),
#            Rational(,),
#            Rational(,)]
#d_phi[6] = [Rational(,),
#            Rational(,),
#            Rational(,),
#            Rational(,),
#            Rational(,),
#            Rational(,),
#            Rational(,)]
#d_phi[7] = [Rational(,),
#            Rational(,),
#            Rational(,),
#            Rational(,),
#            Rational(,),
#            Rational(,),
#            Rational(,),
#            Rational(,)]
#d_phi[8] = [Rational(,),
#            Rational(,),
#            Rational(,),
#            Rational(,),
#            Rational(,),
#            Rational(,),
#            Rational(,),
#            Rational(,),
#            Rational(,)]
## ...

# ............................................

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

            # ... we use nsimplify to get the rational number
            if p <= P_MAX:
                phi = d_phi[p]

            else:
                # ...
                r  = Symbol('r')

                pp = 2*p + 1
                N = pp + 1
                L = list(range(0, N + pp + 1))

                b0 = bspline_basis(pp, L, 0, r)
                bsp = lambdify(r, b0)
                # ...

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
#            print(phi)
#            print('> mass = ', m)

#            # hack to avoid
#            # sympy.tensor.array.dense_ndim_array.ImmutableDenseNDimArray
#            m = sympify(str(m))

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

            if p <= P_MAX:
                phi = d_phi_rr[p]

            else:
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

#            print('> stiffness = ', m)

#            # hack to avoid
#            # sympy.tensor.array.dense_ndim_array.ImmutableDenseNDimArray
#            m = sympify(str(m))

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

            if p <= P_MAX:
                phi = d_phi_r[p]

            else:
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

#            print('> advection = ', m)

#            # hack to avoid
#            # sympy.tensor.array.dense_ndim_array.ImmutableDenseNDimArray
#            m = sympify(str(m))

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

            if p <= P_MAX:
                raise NotImplementedError('Bilaplacian table not available')
#                phi = d_phi_rrrr[p]

            else:

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

#            print('> bilaplacian = ', m)

#            # hack to avoid
#            # sympy.tensor.array.dense_ndim_array.ImmutableDenseNDimArray
#            m = sympify(str(m))

            return m
# ...
