#
# Original file from the UFL project: ufl/demo/ExplicitConvection.ufl
#
# TODO: add dx operator

def a(v, u):
    return dot( dot(w, grad(u)), v ) #* dx
