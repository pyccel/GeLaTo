#
# Original file from the UFL project: ufl/demo/Poisson.ufl
#
# TODO: add dx operator

def a( v, u ):
    return inner(grad(v), grad(u))
