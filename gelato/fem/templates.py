# coding: utf-8

# TODO: - allow for giving a name for the trial/test basis
#       - Ni/Nj should be Ni_0/Nj_O

# .............................................
#          KERNEL     1D case - scalar
# .............................................
template_1d_scalar = """
def {__KERNEL_NAME__}(test_p1, trial_p1, k1, test_bs1, trial_bs1, u1, w1, mat{__ARGS__}{__FIELD_COEFFS__}):
{__FIELD_EVALUATION__}
    mat[:,:] = 0.
    for il_1 in range(0, test_p1+1):
        for jl_1 in range(0, trial_p1+1):
            v = 0.0
            for g1 in range(0, k1):
{__TEST_FUNCTION__}
{__TRIAL_FUNCTION__}
                x    = u1[g1]
                wvol = w1[g1]
{__FIELD_VALUE__}
                v += ({__WEAK_FORM__}) * wvol
            mat[il_1, p1 + jl_1 - il_1] = v
"""

template_header_1d_scalar = '#$ header procedure {__KERNEL_NAME__}(int, int, int, double [:,:,:], double [:,:,:], double [:], double [:], double [:,:]{__TYPES__}{__FIELD_TYPES__})'
# .............................................

# .............................................
#          KERNEL     2D case - scalar
# .............................................
template_2d_scalar = """
def {__KERNEL_NAME__}(test_p1, test_p2, trial_p1, trial_p2, k1, k2, test_bs1, test_bs2, trial_bs1, trial_bs2, u1, u2, w1, w2, mat{__ARGS__}{__FIELD_COEFFS__}):
{__FIELD_EVALUATION__}
    mat[:,:,:,:] = 0.
    for il_1 in range(0, test_p1+1):
        for jl_1 in range(0, trial_p1+1):
            for il_2 in range(0, test_p2+1):
                for jl_2 in range(0, trial_p2+1):
                    v = 0.0
                    for g1 in range(0, k1):
                        for g2 in range(0, k2):
{__TEST_FUNCTION__}
{__TRIAL_FUNCTION__}
                            x    = u1[g1]
                            y    = u2[g2]
                            wvol = w1[g1] * w2[g2]
{__FIELD_VALUE__}
                            v += ({__WEAK_FORM__}) * wvol
                    mat[il_1, il_2, p1 + jl_1 - il_1, p2 + jl_2 - il_2] = v
"""

template_header_2d_scalar = '#$ header procedure {__KERNEL_NAME__}(int, int, int, int, int, int, double [:,:,:], double [:,:,:], double [:,:,:], double [:,:,:], double [:], double [:], double [:], double [:], double [:,:,:,:]{__TYPES__}{__FIELD_TYPES__})'
# .............................................

# .............................................
#          KERNEL     3D case - scalar
# .............................................


template_3d_scalar = """
def {__KERNEL_NAME__}(test_p1, test_p2, test_p3, trial_p1, trial_p2, trial_p3, k1, k2, k3, test_bs1, test_bs2, test_bs3, trial_bs1, trial_bs2, trial_bs3, u1, u2, u3, w1, w2, w3, mat{__ARGS__}{__FIELD_COEFFS__}):
{__FIELD_EVALUATION__}
    mat[:,:,:,:,:,:] = 0.
    for il_1 in range(0, test_p1+1):
        for jl_1 in range(0, trial_p1+1):
            for il_2 in range(0, test_p2+1):
                for jl_2 in range(0, trial_p2+1):
                    for il_3 in range(0, test_p3+1):
                        for jl_3 in range(0, trial_p3+1):
                            v = 0.0
                            for g1 in range(0, k1):
                                for g2 in range(0, k2):
                                    for g3 in range(0, k3):
{__TEST_FUNCTION__}
{__TRIAL_FUNCTION__}
                                        x    = u1[g1]
                                        y    = u2[g2]
                                        z    = u3[g3]
                                        wvol = w1[g1] * w2[g2] * w3[g3]
{__FIELD_VALUE__}
                                        v += ({__WEAK_FORM__}) * wvol
                            mat[il_1, il_2, il_3, p1 + jl_1 - il_1, p2 + jl_2 - il_2, p3 + jl_3 - il_3] = v
"""

template_header_3d_scalar = '#$ header procedure {__KERNEL_NAME__}(int, int, int, int, int, int, int, int, int, double [:,:,:], double [:,:,:], double [:,:,:], double [:,:,:], double [:,:,:], double [:,:,:], double [:], double [:], double [:], double [:], double [:], double [:], double [:,:,:,:,:,:]{__TYPES__}{__FIELD_TYPES__})'
# .............................................

# .............................................
#          KERNEL     2D case - block
# .............................................
template_1d_block = """
def {__KERNEL_NAME__}(p1, k1, basis, u, w, {__MAT_ARGS__}{__ARGS__}{__FIELD_COEFFS__}):
{__FIELD_EVALUATION__}
{__MAT_INIT__}
    for il_1 in range(0, p1+1):
        for jl_1 in range(0, p1+1):
{__ACCUM_INIT__}
            for g1 in range(0, k1):
{__TEST_FUNCTION__}
{__TRIAL_FUNCTION__}
                x    = u[g1]
                wvol = w[g1]
{__FIELD_VALUE__}
{__ACCUM__}
{__ACCUM_ASSIGN__}
"""

template_header_1d_block = '#$ header procedure {__KERNEL_NAME__}(int, int, double [:,:,:], double [:], double [:], {__MAT_TYPES__}{__TYPES__}{__FIELD_TYPES__})'
# .............................................

# .............................................
#          KERNEL     2D case - block
# .............................................
template_2d_block = """
def {__KERNEL_NAME__}(test_p1, test_p2, trial_p1, trial_p2, k1, k2, test_bs1, test_bs2, trial_bs1, trial_bs2, u1, u2, w1, w2, {__MAT_ARGS__}{__ARGS__}{__FIELD_COEFFS__}):
{__FIELD_EVALUATION__}
{__MAT_INIT__}
    for il_1 in range(0, test_p1+1):
        for jl_1 in range(0, trial_p1+1):
            for il_2 in range(0, test_p2+1):
                for jl_2 in range(0, trial_p2+1):
{__ACCUM_INIT__}
                    for g1 in range(0, k1):
                        for g2 in range(0, k2):
{__TEST_FUNCTION__}
{__TRIAL_FUNCTION__}
                            x    = u1[g1]
                            y    = u2[g2]
                            wvol = w1[g1] * w2[g2]
{__FIELD_VALUE__}
{__ACCUM__}
{__ACCUM_ASSIGN__}
"""

template_header_2d_block = '#$ header procedure {__KERNEL_NAME__}(int, int, int, int, int, int, double [:,:,:], double [:,:,:], double [:,:,:], double [:,:,:], double [:], double [:], double [:], double [:], {__MAT_TYPES__}{__TYPES__}{__FIELD_TYPES__})'
# .............................................

# .............................................
#          KERNEL     3D case - block
# .............................................


template_3d_block = """
def {__KERNEL_NAME__}(test_p1, test_p2, test_p3, trial_p1, trial_p2, trial_p3, k1, k2, k3, test_bs1, test_bs2, test_bs3, trial_bs1, trial_bs2, trial_bs3, u1, u2, u3, w1, w2, w3, {__MAT_ARGS__}{__ARGS__}{__FIELD_COEFFS__}):
{__FIELD_EVALUATION__}
{__MAT_INIT__}
    for il_1 in range(0, test_p1+1):
        for jl_1 in range(0, trial_p1+1):
            for il_2 in range(0, test_p2+1):
                for jl_2 in range(0, trial_p2+1):
                    for il_3 in range(0, test_p3+1):
                        for jl_3 in range(0, trial_p3+1):
{__ACCUM_INIT__}
                            for g1 in range(0, k1):
                                for g2 in range(0, k2):
                                    for g3 in range(0, k3):
{__TEST_FUNCTION__}
{__TRIAL_FUNCTION__}
                                        x    = u1[g1]
                                        y    = u2[g2]
                                        z    = u3[g3]
                                        wvol = w1[g1] * w2[g2] * w3[g3]
{__FIELD_VALUE__}
{__ACCUM__}
{__ACCUM_ASSIGN__}
"""

template_header_3d_block = '#$ header procedure {__KERNEL_NAME__}(int, int, int, int, int, int, int, int, int, double [:,:,:], double [:,:,:], double [:,:,:], double [:,:,:], double [:,:,:], double [:,:,:], double [:], double [:], double [:], double [:], double [:], double [:], {__MAT_TYPES__}{__TYPES__}{__FIELD_TYPES__})'
# .............................................

# .............................................
#          SYMBOL     1D case - scalar
# .............................................
symbol_1d_scalar ="""
def {__SYMBOL_NAME__}(arr_x1, arr_t1, symbol{__ARGS__}{__FIELD_COEFFS__}):
    from numpy import sin
    from numpy import cos

{__FIELD_EVALUATION__}

    n1 = len(arr_x1)
    symbol[:] = 0.
    for i1 in range(0, n1):
        x = arr_x1[i1]
        t1 = arr_t1[i1]

{__FIELD_VALUE__}

        symbol[i1] = {__SYMBOL_EXPR__}
"""

symbol_header_1d_scalar = '#$ header procedure {__SYMBOL_NAME__}(double [:], double [:], double [:]{__TYPES__}{__FIELD_TYPES__})'
# .............................................

# .............................................
#          SYMBOL     2D case - scalar
# .............................................
symbol_2d_scalar ="""
def {__SYMBOL_NAME__}(arr_x1, arr_x2, arr_t1, arr_t2, symbol{__ARGS__}{__FIELD_COEFFS__}):
    from numpy import sin
    from numpy import cos

{__FIELD_EVALUATION__}

    n1 = len(arr_x1)
    n2 = len(arr_x2)
    symbol[:,:] = 0.
    for i1 in range(0, n1):
        for i2 in range(0, n2):
            x = arr_x1[i1]
            y = arr_x2[i2]
            t1 = arr_t1[i1]
            t2 = arr_t2[i2]

{__FIELD_VALUE__}

            symbol[i1, i2] = {__SYMBOL_EXPR__}
"""

symbol_header_2d_scalar = '#$ header procedure {__SYMBOL_NAME__}(double [:], double [:], double [:], double [:], double [:,:]{__TYPES__}{__FIELD_TYPES__})'
# .............................................

# .............................................
#          SYMBOL     3D case - scalar
# .............................................
symbol_3d_scalar ="""
def {__SYMBOL_NAME__}(arr_x1, arr_x2, arr_x3, arr_t1, arr_t2, arr_t3, symbol{__ARGS__}{__FIELD_COEFFS__}):
    from numpy import sin
    from numpy import cos

{__FIELD_EVALUATION__}

    n1 = len(arr_x1)
    n2 = len(arr_x2)
    n3 = len(arr_x3)
    symbol[:,:,:] = 0.
    for i1 in range(0, n1):
        for i2 in range(0, n2):
            for i3 in range(0, n3):
                x = arr_x1[i1]
                y = arr_x2[i2]
                z = arr_x3[i3]
                t1 = arr_t1[i1]
                t2 = arr_t2[i2]
                t3 = arr_t3[i3]

{__FIELD_VALUE__}

                symbol[i1, i2, i3] = {__SYMBOL_EXPR__}
"""

symbol_header_3d_scalar = '#$ header procedure {__SYMBOL_NAME__}(double [:], double [:], double [:], double [:], double [:], double [:], double [:,:,:]{__TYPES__}{__FIELD_TYPES__})'
# .............................................

# .............................................
#          SYMBOL     1D case - block
# .............................................
symbol_1d_block ="""
def {__SYMBOL_NAME__}(arr_x1, arr_t1, symbol{__ARGS__}{__FIELD_COEFFS__}):
    from numpy import sin
    from numpy import cos

{__FIELD_EVALUATION__}

    n1 = len(arr_x1)
    symbol[:,:,:] = 0.
    for i1 in range(0, n1):
        x = arr_x1[i1]
        t1 = arr_t1[i1]

{__FIELD_VALUE__}

{__SYMBOL_EXPR__}
"""

symbol_header_1d_block = '#$ header procedure {__SYMBOL_NAME__}(double [:], double [:], double [:,:,:]{__TYPES__}{__FIELD_TYPES__})'
# .............................................

# .............................................
#          SYMBOL     2D case - block
# .............................................
symbol_2d_block ="""
def {__SYMBOL_NAME__}(arr_x1, arr_x2, arr_t1, arr_t2, symbol{__ARGS__}{__FIELD_COEFFS__}):
    from numpy import sin
    from numpy import cos

{__FIELD_EVALUATION__}

    n1 = len(arr_x1)
    n2 = len(arr_x2)
    symbol[:,:,:,:] = 0.
    for i1 in range(0, n1):
        for i2 in range(0, n2):
            x = arr_x1[i1]
            y = arr_x2[i2]
            t1 = arr_t1[i1]
            t2 = arr_t2[i2]

{__FIELD_VALUE__}

{__SYMBOL_EXPR__}
"""

symbol_header_2d_block = '#$ header procedure {__SYMBOL_NAME__}(double [:], double [:], double [:], double [:], double [:,:,:,:]{__TYPES__}{__FIELD_TYPES__})'
# .............................................

# .............................................
#          SYMBOL     3D case - block
# .............................................
symbol_3d_block ="""
def {__SYMBOL_NAME__}(arr_x1, arr_x2, arr_x3, arr_t1, arr_t2, arr_t3, symbol{__ARGS__}{__FIELD_COEFFS__}):
    from numpy import sin
    from numpy import cos

{__FIELD_EVALUATION__}

    n1 = len(arr_x1)
    n2 = len(arr_x2)
    n3 = len(arr_x3)
    symbol[:,:,:,:,:] = 0.
    for i1 in range(0, n1):
        for i2 in range(0, n2):
            for i3 in range(0, n3):
                x = arr_x1[i1]
                y = arr_x2[i2]
                z = arr_x3[i3]
                t1 = arr_t1[i1]
                t2 = arr_t2[i2]
                t3 = arr_t3[i3]

{__FIELD_VALUE__}

{__SYMBOL_EXPR__}
"""

symbol_header_3d_block = '#$ header procedure {__SYMBOL_NAME__}(double [:], double [:], double [:], double [:], double [:], double [:], double [:,:,:,:,:]{__TYPES__}{__FIELD_TYPES__})'
# .............................................

# .............................................
#          FIELD EVALUATION 1D case - scalar
# .............................................
eval_field_1d_scalar = """
    from numpy import zeros

{__FIELD_INIT__}
    for g1 in range(0, k1):
        for jl_1 in range(0, p1+1):
            Nj = bs1[jl_1, 0, g1]
            Nj_x = bs1[jl_1, 1, g1]

{__FIELD_ACCUM__}
"""
# .............................................

# .............................................
#          FIELD EVALUATION 2D case - scalar
# .............................................
eval_field_2d_scalar = """
    from numpy import zeros

{__FIELD_INIT__}
    for g1 in range(0, k1):
        for g2 in range(0, k2):
            for jl_1 in range(0, p1+1):
                for jl_2 in range(0, p2+1):
                    Nj = bs1[jl_1, 0, g1] * bs2[jl_2, 0, g2]
                    Nj_x = bs1[jl_1, 1, g1] * bs2[jl_2, 0, g2]
                    Nj_y = bs1[jl_1, 0, g1] * bs2[jl_2, 1, g2]

{__FIELD_ACCUM__}
"""
# .............................................

# .............................................
#          FIELD EVALUATION 3D case - scalar
# .............................................
eval_field_3d_scalar = """
    from numpy import zeros

{__FIELD_INIT__}
    for g1 in range(0, k1):
        for g2 in range(0, k2):
            for g3 in range(0, k3):
                for jl_1 in range(0, p1+1):
                    for jl_2 in range(0, p2+1):
                        for jl_3 in range(0, p3+1):
                            Nj = bs1[jl_1, 0, g1] * bs2[jl_2, 0, g2] * bs3[jl_3, 0, g3]
                            Nj_x = bs1[jl_1, 1, g1] * bs2[jl_2, 0, g2] * bs3[jl_3, 0, g3]
                            Nj_y = bs1[jl_1, 0, g1] * bs2[jl_2, 1, g2] * bs3[jl_3, 0, g3]
                            Nj_z = bs1[jl_1, 0, g1] * bs2[jl_2, 0, g2] * bs3[jl_3, 1, g3]

{__FIELD_ACCUM__}
"""
# .............................................
