# coding: utf-8

# .............................................
#          ASSEMBLY OF BILINEAR FORM 1D case - scalar
# .............................................
assembly_bilinear_1d_scalar = """
def {__ASSEMBLY_NAME__}( V ):
    {__DOCSTRING__}
    # Sizes
    [s1] = V.vector_space.starts
    [e1] = V.vector_space.ends
    [p1] = V.vector_space.pads

    # Quadrature data
    k1        = V.quad_order
    spans_1   = V.spans
    basis_1   = V.quad_basis
    weights_1 = V.quad_weights

    # Create global matrices
    M = StencilMatrix( V.vector_space, V.vector_space )

    # Create element matrices
    from numpy import zeros
    mat = zeros( (p1+1,2*p1+1) )

    # Build global matrices: cycle over elements
    for ie1 in range(s1, e1+1-p1):

        # Get spline index, B-splines' values and quadrature weights
        is1 =   spans_1[ie1]
        bs1 =   basis_1[ie1,:,:,:]
        w1  = weights_1[ie1,:]

        # Compute element matrices
        {__KERNEL_NAME__}( p1, k1, bs1, w1, mat_m, mat_s )

        # Update global matrices
        M[is1-p1:is1+1,:] += mat[:,:]

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
