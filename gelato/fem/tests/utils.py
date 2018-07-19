# coding: utf-8


def apply_homogeneous_dirichlet_bc_2d(V, matrix=None, rhs=None):
    """ Apply homogeneous dirichlet boundary conditions """

    # assumes a 2D Tensor space
    # add asserts on the space if it is periodic

    V1,V2 = V.spaces

    s1, s2 = V.vector_space.starts
    e1, e2 = V.vector_space.ends

    if not V1.periodic:
        # left  bc at x=0.
        if s1 == 0:
            matrix[0,:,:,:] = 0.
            rhs[0,:] = 0.
        # right bc at x=1.
        if e1 == V1.nbasis-1:
            matrix[e1,:,:,:] = 0.
            rhs [e1,:]     = 0.

    if not V2.periodic:
        # lower bc at y=0.
        if s2 == 0:
            matrix[:,0,:,:] = 0.
            rhs [:,0]     = 0.
        # upper bc at y=1.
        if e2 == V2.nbasis-1:
            matrix[:,e2,:,:] = 0.
            rhs [:,e2] = 0.

    outputs = []
    if matrix:
        outputs.append(matrix)
    if rhs:
        outputs.append(rhs)

    if len(outputs) == 1:
        return outputs[0]

    else:
        return outputs

def apply_homogeneous_dirichlet_bc(V, matrix=None, rhs=None):
    if V.ldim == 1:
        raise NotImplementedError('')

    elif V.ldim == 2:
        return apply_homogeneous_dirichlet_bc_2d(V, matrix=matrix, rhs=rhs)

    elif V.ldim == 3:
        raise NotImplementedError('')
