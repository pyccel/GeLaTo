#
# Original file from the UFL project: ufl/demo/Stiffness.ufl
#
# TODO: add dx operator

def a( v, u ):
    return dot(grad(u), grad(v))
