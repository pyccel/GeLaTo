# coding: utf-8

# .............................................
#          KERNEL     1D case - scalar
# .............................................
_bilinear_form_1d_scalar = """
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

_bilinear_form_header_1d_scalar = '#$ header procedure {__KERNEL_NAME__}(int, int, int, double [:,:,:], double [:,:,:], double [:], double [:], double [:,:]{__TYPES__}{__FIELD_TYPES__})'
# .............................................

# .............................................
#          KERNEL     2D case - scalar
# .............................................
_bilinear_form_2d_scalar = """
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

_bilinear_form_header_2d_scalar = '#$ header procedure {__KERNEL_NAME__}(int, int, int, int, int, int, double [:,:,:], double [:,:,:], double [:,:,:], double [:,:,:], double [:], double [:], double [:], double [:], double [:,:,:,:]{__TYPES__}{__FIELD_TYPES__})'
# .............................................

# .............................................
#          KERNEL     3D case - scalar
# .............................................


_bilinear_form_3d_scalar = """
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

_bilinear_form_header_3d_scalar = '#$ header procedure {__KERNEL_NAME__}(int, int, int, int, int, int, int, int, int, double [:,:,:], double [:,:,:], double [:,:,:], double [:,:,:], double [:,:,:], double [:,:,:], double [:], double [:], double [:], double [:], double [:], double [:], double [:,:,:,:,:,:]{__TYPES__}{__FIELD_TYPES__})'
# .............................................

# .............................................
#          KERNEL     2D case - block
# .............................................
_bilinear_form_1d_block = """
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

_bilinear_form_header_1d_block = '#$ header procedure {__KERNEL_NAME__}(int, int, double [:,:,:], double [:], double [:], {__MAT_TYPES__}{__TYPES__}{__FIELD_TYPES__})'
# .............................................

# .............................................
#          KERNEL     2D case - block
# .............................................
_bilinear_form_2d_block = """
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

_bilinear_form_header_2d_block = '#$ header procedure {__KERNEL_NAME__}(int, int, int, int, int, int, double [:,:,:], double [:,:,:], double [:,:,:], double [:,:,:], double [:], double [:], double [:], double [:], {__MAT_TYPES__}{__TYPES__}{__FIELD_TYPES__})'
# .............................................

# .............................................
#          KERNEL     3D case - block
# .............................................
_bilinear_form_3d_block = """
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

_bilinear_form_header_3d_block = '#$ header procedure {__KERNEL_NAME__}(int, int, int, int, int, int, int, int, int, double [:,:,:], double [:,:,:], double [:,:,:], double [:,:,:], double [:,:,:], double [:,:,:], double [:], double [:], double [:], double [:], double [:], double [:], {__MAT_TYPES__}{__TYPES__}{__FIELD_TYPES__})'
# .............................................
