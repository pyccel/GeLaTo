#
# Original file from the UFL project: ufl/demo/MixedPoisson.ufl
#
# TODO: add dx operator

def a((tau, w), (sigma, u)):
    return (dot(tau, sigma) - div(tau)*u + w*div(sigma)) #*dx

