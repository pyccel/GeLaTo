# coding: utf-8

# .............................................
#          KERNEL     1D case - scalar
# .............................................
_function_form_1d_scalar = """
def {__KERNEL_NAME__}(test_p1, k1, test_bs1, u1, w1{__ARGS__}{__FIELD_COEFFS__}):
{__FIELD_EVALUATION__}
    mat = 0.
    for g1 in range(0, k1):
        x    = u1[g1]
        wvol = w1[g1]
{__FIELD_VALUE__}
        mat += ({__WEAK_FORM__}) * wvol
    return mat
"""

_function_form_header_1d_scalar = '#$ header procedure {__KERNEL_NAME__}(int, int, double [:,:,:], double [:], double [:], double [:,:]{__TYPES__}{__FIELD_TYPES__})'
# .............................................

# .............................................
#          KERNEL     2D case - scalar
# .............................................
_function_form_2d_scalar = """
def {__KERNEL_NAME__}(test_p1, test_p2, k1, k2, test_bs1, test_bs2, u1, u2, w1, w2{__ARGS__}{__FIELD_COEFFS__}):
{__FIELD_EVALUATION__}
    mat = 0.
    for g1 in range(0, k1):
        for g2 in range(0, k2):
            x    = u1[g1]
            y    = u2[g2]
            wvol = w1[g1] * w2[g2]
{__FIELD_VALUE__}
            mat += ({__WEAK_FORM__}) * wvol
    return mat
"""

_function_form_header_2d_scalar = '#$ header procedure {__KERNEL_NAME__}(int, int, int, int, double [:,:,:], double [:,:,:], double [:], double [:], double [:], double [:], double [:,:,:,:]{__TYPES__}{__FIELD_TYPES__})'
# .............................................

# .............................................
#          KERNEL     3D case - scalar
# .............................................
_function_form_3d_scalar = """
def {__KERNEL_NAME__}(test_p1, test_p2, test_p3, k1, k2, k3, test_bs1, test_bs2, test_bs3, u1, u2, u3, w1, w2, w3{__ARGS__}{__FIELD_COEFFS__}):
{__FIELD_EVALUATION__}
    mat = 0.
    for g1 in range(0, k1):
        for g2 in range(0, k2):
            for g3 in range(0, k3):
                x    = u1[g1]
                y    = u2[g2]
                z    = u3[g3]
                wvol = w1[g1] * w2[g2] * w3[g3]
{__FIELD_VALUE__}
                mat += ({__WEAK_FORM__}) * wvol
    return mat
"""

_function_form_header_3d_scalar = '#$ header procedure {__KERNEL_NAME__}(int, int, int, int, int, int, double [:,:,:], double [:,:,:], double [:,:,:], double [:], double [:], double [:], double [:], double [:], double [:], double [:,:,:,:,:,:]{__TYPES__}{__FIELD_TYPES__})'
# .............................................
