#
# Original file from the UFL project: ufl/demo/Poisson.ufl
#
# TODO: add dx operator

def a( W, V ):
    w = TestFunction(W, name='w')
    v = TrialFunction(V, name='v')

    return inner(grad(w), grad(v))
