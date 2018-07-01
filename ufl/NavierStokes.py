#
# Original file from the UFL project: ufl/demo/NavierStokes.ufl
#
# TODO: add dx operator

def a( v, u ):
    Du = grad(u)
    return dot( dot(w, Du), v ) #*dx
