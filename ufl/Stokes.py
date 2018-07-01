#
# Original file from the UFL project: ufl/demo/Stokes.ufl
#
# TODO: add dx operator

def a( (v, q), (u, p) ):
    return (inner(grad(v), grad(u)) - div(v)*p + q*div(u))

