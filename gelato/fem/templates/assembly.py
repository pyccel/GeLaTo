# coding: utf-8

# .............................................
#          ASSEMBLY OF BILINEAR FORM 1D case - scalar
# .............................................
_assembly_bilinear_1d_scalar = """
def {__ASSEMBLY_NAME__}( self, test_space, trial_space ):
    {__DOCSTRING__}
    # Sizes
    [s1] = test_space.vector_space.starts
    [e1] = test_space.vector_space.ends

    [test_p1] = test_space.vector_space.pads
    [trial_p1] = trial_space.vector_space.pads

    # Quadrature data
    k1        = test_space.quad_order
    points_1  = test_space.quad_points
    weights_1 = test_space.quad_weights

    test_spans_1   = test_space.spans
    trial_spans_1  = trial_space.spans

    test_basis_1   = test_space.quad_basis
    trial_basis_1  = trial_space.quad_basis

    # Create global matrices
    from spl.linalg.stencil import StencilMatrix
    M = StencilMatrix( test_space.vector_space, trial_space.vector_space )

    # Create element matrices
    from numpy import zeros
    mat = zeros( (test_p1+1,2*trial_p1+1) )

    # Build global matrices: cycle over elements
    for ie1 in range(s1, e1+1-test_p1):

        # Get spline index, B-splines' values and quadrature weights
        is1 =   test_spans_1[ie1]
        w1  = weights_1[ie1,:]
        u1  = points_1[ie1, :]

        test_bs1  = test_basis_1[ie1,:,:,:]
        trial_bs1 = trial_basis_1[ie1,:,:,:]

        # Compute element matrices
        {__KERNEL_NAME__}( test_p1, trial_p1, k1, test_bs1, trial_bs1, w1, u1, mat )

        # Update global matrices
        M[is1-test_p1:is1+1,:] += mat[:,:]

    return M
"""
# .............................................

# .............................................
#          ASSEMBLY OF BILINEAR FORM 2D case - scalar
# .............................................
# .............................................

# .............................................
#          ASSEMBLY OF BILINEAR FORM 3D case - scalar
# .............................................
# .............................................
