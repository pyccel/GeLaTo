#
# Original file from the UFL project: ufl/demo/NonlinearPoisson.ufl
#
# TODO: add dx operator

def a( v, u ):
    return (1+u0**2)*dot(grad(v), grad(u)) + 2*u0*u*dot(grad(v), grad(u0))
