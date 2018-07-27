# coding: utf-8

# .............................................
#          SYMBOL     1D case - scalar
# .............................................
_symbol_1d_scalar ="""
def {__SYMBOL_NAME__}(arr_x1, arr_t1, mat{__ARGS__}{__FIELD_COEFFS__}):
    from numpy import sin
    from numpy import cos
{__FIELD_EVALUATION__}
    n1 = len(arr_x1)
    mat[:] = 0.
    for i1 in range(0, n1):
        x = arr_x1[i1]
        tx = arr_t1[i1]
{__FIELD_VALUE__}
        mat[i1] = {__SYMBOL_EXPR__}
"""

_symbol_header_1d_scalar = '#$ header procedure {__SYMBOL_NAME__}(double [:], double [:], double [:]{__TYPES__}{__FIELD_TYPES__})'
# .............................................

# .............................................
#          SYMBOL     2D case - scalar
# .............................................
_symbol_2d_scalar ="""
def {__SYMBOL_NAME__}(arr_x1, arr_x2, arr_t1, arr_t2, mat{__ARGS__}{__FIELD_COEFFS__}):
    from numpy import sin
    from numpy import cos
{__FIELD_EVALUATION__}
    n1 = len(arr_x1)
    n2 = len(arr_x2)
    mat[:,:] = 0.
    for i1 in range(0, n1):
        for i2 in range(0, n2):
            x = arr_x1[i1]
            y = arr_x2[i2]
            tx = arr_t1[i1]
            ty = arr_t2[i2]
{__FIELD_VALUE__}
            mat[i1, i2] = {__SYMBOL_EXPR__}
"""

_symbol_header_2d_scalar = '#$ header procedure {__SYMBOL_NAME__}(double [:], double [:], double [:], double [:], double [:,:]{__TYPES__}{__FIELD_TYPES__})'
# .............................................

# .............................................
#          SYMBOL     3D case - scalar
# .............................................
_symbol_3d_scalar ="""
def {__SYMBOL_NAME__}(arr_x1, arr_x2, arr_x3, arr_t1, arr_t2, arr_t3, mat{__ARGS__}{__FIELD_COEFFS__}):
    from numpy import sin
    from numpy import cos
{__FIELD_EVALUATION__}
    n1 = len(arr_x1)
    n2 = len(arr_x2)
    n3 = len(arr_x3)
    mat[:,:,:] = 0.
    for i1 in range(0, n1):
        for i2 in range(0, n2):
            for i3 in range(0, n3):
                x = arr_x1[i1]
                y = arr_x2[i2]
                z = arr_x3[i3]
                tx = arr_t1[i1]
                ty = arr_t2[i2]
                tz = arr_t3[i3]
{__FIELD_VALUE__}
                mat[i1, i2, i3] = {__SYMBOL_EXPR__}
"""

_symbol_header_3d_scalar = '#$ header procedure {__SYMBOL_NAME__}(double [:], double [:], double [:], double [:], double [:], double [:], double [:,:,:]{__TYPES__}{__FIELD_TYPES__})'
# .............................................

# .............................................
#          SYMBOL     1D case - block
# .............................................
_symbol_1d_block ="""
def {__SYMBOL_NAME__}(arr_x1, arr_t1, mat{__ARGS__}{__FIELD_COEFFS__}):
    from numpy import sin
    from numpy import cos
{__FIELD_EVALUATION__}
    n1 = len(arr_x1)
    mat[:,:,:] = 0.
    for i1 in range(0, n1):
        x = arr_x1[i1]
        tx = arr_t1[i1]
{__FIELD_VALUE__}
{__SYMBOL_EXPR__}
"""

_symbol_header_1d_block = '#$ header procedure {__SYMBOL_NAME__}(double [:], double [:], double [:,:,:]{__TYPES__}{__FIELD_TYPES__})'
# .............................................

# .............................................
#          SYMBOL     2D case - block
# .............................................
_symbol_2d_block ="""
def {__SYMBOL_NAME__}(arr_x1, arr_x2, arr_t1, arr_t2, mat{__ARGS__}{__FIELD_COEFFS__}):
    from numpy import sin
    from numpy import cos
{__FIELD_EVALUATION__}
    n1 = len(arr_x1)
    n2 = len(arr_x2)
    mat[:,:,:,:] = 0.
    for i1 in range(0, n1):
        for i2 in range(0, n2):
            x = arr_x1[i1]
            y = arr_x2[i2]
            tx = arr_t1[i1]
            ty = arr_t2[i2]
{__FIELD_VALUE__}
{__SYMBOL_EXPR__}
"""

_symbol_header_2d_block = '#$ header procedure {__SYMBOL_NAME__}(double [:], double [:], double [:], double [:], double [:,:,:,:]{__TYPES__}{__FIELD_TYPES__})'
# .............................................

# .............................................
#          SYMBOL     3D case - block
# .............................................
_symbol_3d_block ="""
def {__SYMBOL_NAME__}(arr_x1, arr_x2, arr_x3, arr_t1, arr_t2, arr_t3, mat{__ARGS__}{__FIELD_COEFFS__}):
    from numpy import sin
    from numpy import cos
{__FIELD_EVALUATION__}
    n1 = len(arr_x1)
    n2 = len(arr_x2)
    n3 = len(arr_x3)
    mat[:,:,:,:,:] = 0.
    for i1 in range(0, n1):
        for i2 in range(0, n2):
            for i3 in range(0, n3):
                x = arr_x1[i1]
                y = arr_x2[i2]
                z = arr_x3[i3]
                tx = arr_t1[i1]
                ty = arr_t2[i2]
                tz = arr_t3[i3]
{__FIELD_VALUE__}
{__SYMBOL_EXPR__}
"""

_symbol_header_3d_block = '#$ header procedure {__SYMBOL_NAME__}(double [:], double [:], double [:], double [:], double [:], double [:], double [:,:,:,:,:]{__TYPES__}{__FIELD_TYPES__})'
# .............................................
