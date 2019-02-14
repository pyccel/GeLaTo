# coding: utf-8

from sympy import Symbol
from sympy.core.containers import Tuple
from sympy import symbols
from sympy import pi, cos, sin
from sympy import srepr
from sympy import I

from sympde.core import Constant
from sympde.calculus import grad, dot, inner, cross, rot, curl, div
from sympde.calculus import laplace, hessian, bracket, convect
from sympde.topology import (dx, dy, dz)
from sympde.topology import FunctionSpace, VectorFunctionSpace
from sympde.topology import Field, TestFunction
from sympde.topology import Domain
from sympde.topology import Trace, trace_0, trace_1
from sympde.topology import Mapping
from sympde.topology import Square
from sympde.expr.expr import LinearForm, BilinearForm

from gelato import gelatize
from gelato import (Mass,
                    Stiffness,
                    Advection,
                    Bilaplacian)

DIM = 2
domain = Domain('Omega', dim=DIM)

#==============================================================================
def test_gelatize_2d_1():

    V = FunctionSpace('V', domain)

    v = TestFunction(V, name='v')
    u = TestFunction(V, name='u')

    nx, ny = symbols('nx ny', integer=True)
    px, py = symbols('px py', integer=True)
    tx, ty = symbols('tx ty')


    c = Constant('c')

    bx = Constant('bx')
    by = Constant('by')
    b = Tuple(bx, by)

    # ...
    expected = Mass(px,tx)*Mass(py,ty)/(nx*ny)
    assert(gelatize(BilinearForm((u,v), u*v)) == expected)
    # ...

    # ...
    expected = Mass(py,ty)*nx*Stiffness(px,tx)/ny
    assert(gelatize(BilinearForm((u,v), dx(u)*dx(v))) == expected)
    # ...

    # ...
    expected = I*Advection(py,ty)*Mass(px,tx)/nx
    assert(gelatize(BilinearForm((u,v), dy(u) * v)) == expected)
    # ...

    # ...
    expected = I*Advection(px,tx)*Mass(py,ty)/ny
    assert(gelatize(BilinearForm((u,v), dx(u) * v)) == expected)
    # ...

    # ...
    expected = Mass(px,tx)*ny*Stiffness(py,ty)/nx + Mass(py,ty)*nx*Stiffness(px,tx)/ny
    assert(gelatize(BilinearForm((u,v), dot(grad(v), grad(u)))) == expected)
    # ...

    # ...
    expected = (nx*Mass(py,ty)*Stiffness(px,tx)/ny +
                I*Advection(px,tx)*Mass(py,ty)/ny +
                ny*Mass(px,tx)*Stiffness(py,ty)/nx +
                I*Advection(py,ty)*Mass(px,tx)/nx)
    assert(gelatize(BilinearForm((u,v), dot(grad(v), grad(u)) + dx(u)*v + dy(u)*v)) == expected)
    # ...

    # ...
    expected = -bx*I*Advection(px,tx)*Mass(py,ty)/ny - by*I*Advection(py,ty)*Mass(px,tx)/nx
    assert(gelatize(BilinearForm((u,v), dot(b, grad(v)) * u)) == expected)
    # ...

#    # ... TODO
#    expected = bx**2*nx*Mass(py,ty)*Stiffness(px,tx)/ny + 2*bx*by*Advection(px,tx)*Advection(py,ty) + by**2*ny*Mass(px,tx)*Stiffness(py,ty)/nx
#    assert(gelatize(BilinearForm((u,v), dot(b, grad(v)) * dot(b, grad(u)))) == expected)
#    # ...

#    degrees = None
##    degrees = [2, 1]

##    evaluate = True
#    evaluate = False

#    expr = dot(b, grad(v)) * dot(b, grad(u))
#    expr = dot(grad(v), grad(u)) + dx(u)*v + dy(u)*v
#
#    expr = BilinearForm((u,v), expr)
#    print('> input     >>> {0}'.format(expr))
#    print('> gelatized >>> {0}'.format(gelatize(expr, degrees=degrees)))

#==============================================================================
def test_gelatize_2d_3():

    V = VectorFunctionSpace('V', domain)

    v = TestFunction(V, name='v')
    u = TestFunction(V, name='u')

    c = Constant('c')

    expr = c * div(v) * div(u) + rot(v) * rot(u)
    expr = BilinearForm((u,v), expr)
    print('> input     >>> {0}'.format(expr))
    print('> gelatized >>> {0}'.format(gelatize(expr)))

#==============================================================================
def test_gelatize_2d_4():

    V = FunctionSpace('V', domain)

    v = TestFunction(V, name='v')
    u = TestFunction(V, name='u')

    expr = BilinearForm((u,v), dot(grad(v), grad(u)))

    print('> input     >>> {0}'.format(expr))
    print('> gelatized >>> {0}'.format(gelatize(expr)))

#==============================================================================
def test_gelatize_2d_5():

    V = FunctionSpace('V', domain)

    v = TestFunction(V, name='v')
    u = TestFunction(V, name='u')

    c = Constant('c')

    expr = BilinearForm((u,v), dot(grad(v), grad(u)) + c*v*u)

    print('> input     >>> {0}'.format(expr))
    print('> gelatized >>> {0}'.format(gelatize(expr)))

#==============================================================================
def test_gelatize_2d_5_mapping():

    M = Mapping('M', DIM)

    V = FunctionSpace('V', domain)

    v = TestFunction(V, name='v')
    u = TestFunction(V, name='u')

    c = Constant('c')

    expr = BilinearForm((u,v), dot(grad(v), grad(u)) + c*v*u)

    print('> input     >>> {0}'.format(expr))
    print('> gelatized >>> {0}'.format(gelatize(expr, mapping=M, human=True)))

#==============================================================================
def test_gelatize_2d_3_mapping():

    M = Mapping('M', DIM)

    V = VectorFunctionSpace('V', domain)

    v = TestFunction(V, name='v')
    u = TestFunction(V, name='u')

    c = Constant('c')

    expr = c * div(v) * div(u) + rot(v) * rot(u)
    expr = BilinearForm((u,v), expr)
    print('> input     >>> {0}'.format(expr))
    print('> gelatized >>> {0}'.format(gelatize(expr, mapping=M, human=True)))


#==============================================================================
# CLEAN UP SYMPY NAMESPACE
#==============================================================================

def teardown_module():
    from sympy import cache
    cache.clear_cache()

def teardown_function():
    from sympy import cache
    cache.clear_cache()

#test_gelatize_2d_1()
#test_gelatize_2d_3()
#test_gelatize_2d_4()
#test_gelatize_2d_5()
#test_gelatize_2d_5_mapping()
#test_gelatize_2d_3_mapping()
