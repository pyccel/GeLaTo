# coding: utf-8
from numpy import zeros
from collections import OrderedDict

def assemble_matrix(V, kernel, args=None, M=None, fields={}):
    try:
        assemble = eval('assemble_matrix_{}d'.format(V.pdim))
    except:
        raise ValueError('Could not find the corresponding assembly procedure')


    if not( args is None ):
        if isinstance(args, (dict, OrderedDict)):
            d_args = OrderedDict(sorted(list(args.items())))
            args = list(d_args.values())

        if not(isinstance(args, (list, tuple))):
            raise TypeError('Expecting list/tuple for args')

    if fields:
        from spl.fem.splines import Spline

        assert(isinstance(fields, (dict, OrderedDict)))

        for k,f in list(fields.items()):
            # check validity of arguments
            if not isinstance(f, Spline):
                raise TypeError('Field must be a spl.fem.splines.Spline')

    return assemble(V, kernel, args=args, M=M, fields=fields)


def assemble_matrix_1d(V, kernel, args=None, M=None, fields={}):

    from spl.fem.vector  import VectorBasicSobolevSpace

    # ...
    is_block = False
    if isinstance(V, VectorBasicSobolevSpace):
        if not( V.is_block ):
            raise NotImplementedError('Expecting a scalar or vector space.')

        is_block = True
    # ...

    # ... data structure if not given
    if M is None:
        from spl.linalg.stencil import StencilVector, StencilMatrix

        # scalar case
        if not is_block:
            M = StencilMatrix(V.vector_space, V.vector_space)
        else:
            n_components = len(V.spaces)

            V = V.spaces[0]

            M = []
            for i in range(0, n_components):
                line = []
                for j in range(0, n_components):
                    line.append(StencilMatrix(V.vector_space, V.vector_space))
                M.append(line)
    # ...

    # ... sizes
    [s1] = V.vector_space.starts
    [e1] = V.vector_space.ends
    [p1] = V.vector_space.pads

    k1 = V.quad_order
    spans_1 = V.spans
    basis_1 = V.basis
    weights_1 = V.weights
    points_1 = V.points
    # ...

    # ... element matrix
    if not is_block:
        mat = zeros((p1+1,2*p1+1), order='F')
    else:
        mats = []
        for i in range(0, n_components):
            line = []
            for j in range(0, n_components):
                line.append(zeros((p1+1,2*p1+1), order='F'))
            mats.append(line)
    # ...

    # ... build matrices
    # TODO this is only for the parallel case
#    for ie1 in range(s1, e1+1-p1):
    for ie1 in range(0, V.ncells):
        i_span_1 = spans_1[ie1]

        bs = basis_1[:, :, :, ie1]
        w = weights_1[:, ie1]
        u = points_1[:, ie1]

        s1 = i_span_1 - p1 - 1

        field_coeffs = []
        for k,f in list(fields.items()):
            field_coeffs.append(f.coeffs[s1:s1+p1+1])

        if not is_block:
            if args is None:
                kernel(p1, k1, bs, u, w, mat, *field_coeffs)

            else:
                kernel(p1, k1, bs, u, w, mat, *args, *field_coeffs)

        else:
            _mats = []
            for i in range(0, n_components):
                for j in range(0, n_components):
                    _mats.append(mats[i][j])

            if args is None:
                kernel(p1, k1, bs, u, w, *_mats, *field_coeffs)

            else:
                kernel(p1, k1, bs, u, w, *_mats, *args, *field_coeffs)

        if not is_block:
            M._data[s1:s1+p1+1,:] += mat[:,:]
        else:
            for i in range(0, n_components):
                for j in range(0, n_components):
                    _mat = mats[i][j]
                    _M = M[i][j]

            _M._data[s1:s1+p1+1,:] += _mat[:,:]
    # ...

    return M


def assemble_matrix_2d(V, kernel, args=None, M=None, fields={}):

    from spl.fem.vector  import VectorBasicSobolevSpace

    # ...
    is_block = False
    if isinstance(V, VectorBasicSobolevSpace):
        if not( V.is_block ):
            raise NotImplementedError('Expecting a scalar or vector space.')

        is_block = True
    # ...

    # ... data structure if not given
    if M is None:
        from spl.linalg.stencil import StencilVector, StencilMatrix

        # scalar case
        if not is_block:
            M = StencilMatrix(V.vector_space, V.vector_space)
        else:
            n_components = len(V.spaces)

            V = V.spaces[0]

            M = []
            for i in range(0, n_components):
                line = []
                for j in range(0, n_components):
                    line.append(StencilMatrix(V.vector_space, V.vector_space))
                M.append(line)
    # ...

    # ... sizes
    [s1, s2] = V.vector_space.starts
    [e1, e2] = V.vector_space.ends
    [p1, p2] = V.vector_space.pads
    # ...

    # ... seetings
    [k1, k2] = [W.quad_order for W in V.spaces]
    [spans_1, spans_2] = [W.spans for W in V.spaces]
    [basis_1, basis_2] = [W.basis for W in V.spaces]
    [weights_1, weights_2] = [W.weights for W in V.spaces]
    [points_1, points_2] = [W.points for W in V.spaces]
    # ...

    # ... element matrix
    if not is_block:
        mat = zeros((p1+1, p2+1, 2*p1+1, 2*p2+1), order='F')
    else:
        mats = []
        for i in range(0, n_components):
            line = []
            for j in range(0, n_components):
                line.append(zeros((p1+1, p2+1, 2*p1+1, 2*p2+1), order='F'))
            mats.append(line)
    # ...

    # ... build matrices
    # TODO this is only for the parallel case
#    for ie1 in range(s1, e1+1-p1):
#        for ie2 in range(s2, e2+1-p2):
    for ie1 in range(0, V.ncells[0]):
        for ie2 in range(0, V.ncells[1]):
            i_span_1 = spans_1[ie1]
            i_span_2 = spans_2[ie2]

            bs1 = basis_1[:, :, :, ie1]
            bs2 = basis_2[:, :, :, ie2]
            w1 = weights_1[:, ie1]
            w2 = weights_2[:, ie2]
            u1 = points_1[:, ie1]
            u2 = points_2[:, ie2]

            s1 = i_span_1 - p1 - 1
            s2 = i_span_2 - p2 - 1

            field_coeffs = []
            for k,f in list(fields.items()):
                field_coeffs.append(f.coeffs[s1:s1+p1+1,s2:s2+p2+1])

            if not is_block:
                if args is None:
                    kernel(p1, p2,
                           k1, k2,
                           bs1, bs2,
                           u1, u2,
                           w1, w2,
                           mat, *field_coeffs)

                else:
                    kernel(p1, p2,
                           k1, k2,
                           bs1, bs2,
                           u1, u2,
                           w1, w2,
                           mat, *args, *field_coeffs)

            else:
                _mats = []
                for i in range(0, n_components):
                    for j in range(0, n_components):
                        _mats.append(mats[i][j])

                if args is None:
                    kernel(p1, p2,
                           k1, k2,
                           bs1, bs2,
                           u1, u2,
                           w1, w2,
                           *_mats, *field_coeffs)

                else:
                    kernel(p1, p2,
                           k1, k2,
                           bs1, bs2,
                           u1, u2,
                           w1, w2,
                           *_mats, *args, *field_coeffs)

            if not is_block:
                M._data[s1:s1+p1+1,s2:s2+p2+1,:,:] += mat[:,:,:,:]
            else:
                for i in range(0, n_components):
                    for j in range(0, n_components):
                        _mat = mats[i][j]
                        _M = M[i][j]

                        _M._data[s1:s1+p1+1,s2:s2+p2+1,:,:] += _mat[:,:,:,:]
    # ...

    return M


def assemble_matrix_3d(V, kernel, args=None, M=None, fields={}):

    from spl.fem.vector  import VectorBasicSobolevSpace

    # ...
    is_block = False
    if isinstance(V, VectorBasicSobolevSpace):
        if not( V.is_block ):
            raise NotImplementedError('Expecting a scalar or vector space.')

        is_block = True
    # ...

    # ... data structure if not given
    if M is None:
        from spl.linalg.stencil import StencilVector, StencilMatrix

        # scalar case
        if not is_block:
            M = StencilMatrix(V.vector_space, V.vector_space)
        else:
            n_components = len(V.spaces)

            V = V.spaces[0]

            M = []
            for i in range(0, n_components):
                line = []
                for j in range(0, n_components):
                    line.append(StencilMatrix(V.vector_space, V.vector_space))
                M.append(line)
    # ...

    # ... sizes
    [s1, s2, s3] = V.vector_space.starts
    [e1, e2, e3] = V.vector_space.ends
    [p1, p2, p3] = V.vector_space.pads
    # ...

    # ... seetings
    [k1, k2, k3] = [W.quad_order for W in V.spaces]
    [spans_1, spans_2, spans_3] = [W.spans for W in V.spaces]
    [basis_1, basis_2, basis_3] = [W.basis for W in V.spaces]
    [weights_1, weights_2, weights_3] = [W.weights for W in V.spaces]
    [points_1, points_2, points_3] = [W.points for W in V.spaces]
    # ...

    # ... element matrix
    if not is_block:
        mat = zeros((p1+1, p2+1, p3+1, 2*p1+1, 2*p2+1, 2*p3+1), order='F')
    else:
        mats = []
        for i in range(0, n_components):
            line = []
            for j in range(0, n_components):
                line.append(zeros((p1+1, p2+1, p3+1, 2*p1+1, 2*p2+1, 2*p3+1), order='F'))
            mats.append(line)
    # ...

    # ... build matrices
    # TODO this is only for the parallel case
#    for ie1 in range(s1, e1+1-p1):
#        for ie2 in range(s2, e2+1-p2):
#            for ie3 in range(s3, e3+1-p3):
    for ie1 in range(0, V.ncells[0]):
        for ie2 in range(0, V.ncells[1]):
            for ie3 in range(0, V.ncells[2]):
                i_span_1 = spans_1[ie1]
                i_span_2 = spans_2[ie2]
                i_span_3 = spans_3[ie3]

                bs1 = basis_1[:, :, :, ie1]
                bs2 = basis_2[:, :, :, ie2]
                bs3 = basis_3[:, :, :, ie3]
                w1 = weights_1[:, ie1]
                w2 = weights_2[:, ie2]
                w3 = weights_3[:, ie3]
                u1 = points_1[:, ie1]
                u2 = points_2[:, ie2]
                u3 = points_3[:, ie3]

                s1 = i_span_1 - p1 - 1
                s2 = i_span_2 - p2 - 1
                s3 = i_span_3 - p3 - 1

                field_coeffs = []
                for k,f in list(fields.items()):
                    field_coeffs.append(f.coeffs[s1:s1+p1+1,s2:s2+p2+1,s3:s3+p3+1])

                if not is_block:
                    if args is None:
                        kernel(p1, p2, p3,
                               k1, k2, k3,
                               bs1, bs2, bs3,
                               u1, u2, u3,
                               w1, w2, w3,
                               mat, *field_coeffs)

                    else:
                        kernel(p1, p2, p3,
                               k1, k2, k3,
                               bs1, bs2, bs3,
                               u1, u2, u3,
                               w1, w2, w3,
                               mat, *args, *field_coeffs)

                else:
                    _mats = []
                    for i in range(0, n_components):
                        for j in range(0, n_components):
                            _mats.append(mats[i][j])

                    if args is None:
                        kernel(p1, p2, p3,
                               k1, k2, k3,
                               bs1, bs2, bs3,
                               u1, u2, u3,
                               w1, w2, w3,
                               *_mats, *field_coeffs)
                    else:
                        kernel(p1, p2, p3,
                               k1, k2, k3,
                               bs1, bs2, bs3,
                               u1, u2, u3,
                               w1, w2, w3,
                               *_mats, *args, *field_coeffs)

                if not is_block:
                    M._data[s1:s1+p1+1,s2:s2+p2+1,s3:s3+p3+1,:,:,:] += mat[:,:,:,:,:,:]
                else:
                    for i in range(0, n_components):
                        for j in range(0, n_components):
                            _mat = mats[i][j]
                            _M = M[i][j]

                            _M._data[s1:s1+p1+1,s2:s2+p2+1,s3:s3+p3+1,:,:,:] += _mat[:,:,:,:,:,:]
    # ...

    return M
