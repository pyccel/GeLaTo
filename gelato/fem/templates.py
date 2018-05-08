# coding: utf-8

# TODO: - allow for giving a name for the trial/test basis
#       - Ni/Nj should be Ni_0/Nj_O

# .............................................
#          KERNEL     1D case - scalar
# .............................................
template_1d_scalar = """
def {__KERNEL_NAME__}(p1, k1, basis, u, w, mat {__ARGS__}):
    mat[:,:] = 0.
    for il_1 in range(0, p1+1):
        for jl_1 in range(0, p1+1):

            v = 0.0
            for g1 in range(0, k1):
                Ni = basis[il_1, 0, g1]
                Ni_x = basis[il_1, 1, g1]

                Nj = basis[jl_1, 0, g1]
                Nj_x = basis[jl_1, 1, g1]

                x    = u[g1]
                wvol = w[g1]

                v += ({__WEAK_FORM__}) * wvol
            mat[il_1, p1 + jl_1 - il_1] = v
"""

template_header_1d_scalar = '#$ header procedure {__KERNEL_NAME__}(int, int, double [:,:,:], double [:], double [:], double [:,:] {__TYPES__})'
# .............................................

# .............................................
#          KERNEL     2D case - scalar
# .............................................
template_2d_scalar = """
def {__KERNEL_NAME__}(p1, p2, k1, k2, bs1, bs2, u1, u2, w1, w2, mat {__ARGS__}):
    mat[:,:,:,:] = 0.
    for il_1 in range(0, p1+1):
        for jl_1 in range(0, p1+1):
            for il_2 in range(0, p2+1):
                for jl_2 in range(0, p2+1):

                    v = 0.0
                    for g1 in range(0, k1):
                        for g2 in range(0, k2):
                            Ni = bs1[il_1, 0, g1] * bs2[il_2, 0, g2]
                            Ni_x = bs1[il_1, 1, g1] * bs2[il_2, 0, g2]
                            Ni_y = bs1[il_1, 0, g1] * bs2[il_2, 1, g2]

                            Nj = bs1[jl_1, 0, g1] * bs2[jl_2, 0, g2]
                            Nj_x = bs1[jl_1, 1, g1] * bs2[jl_2, 0, g2]
                            Nj_y = bs1[jl_1, 0, g1] * bs2[jl_2, 1, g2]

                            x    = u1[g1]
                            y    = u2[g2]
                            wvol = w1[g1] * w2[g2]

                            v += ({__WEAK_FORM__}) * wvol
                    mat[il_1, il_2, p1 + jl_1 - il_1, p2 + jl_2 - il_2] = v
"""

template_header_2d_scalar = '#$ header procedure {__KERNEL_NAME__}(int, int, int, int, double [:,:,:], double [:,:,:], double [:], double [:], double [:], double [:], double [:,:,:,:] {__TYPES__})'
# .............................................

# .............................................
#          KERNEL     3D case - scalar
# .............................................
template_3d_scalar = """
def {__KERNEL_NAME__}(p1, p2, p3, k1, k2, k3, bs1, bs2, bs3, u1, u2, u3, w1, w2, w3, mat {__ARGS__}):
    mat[:,:,:,:,:,:] = 0.
    for il_1 in range(0, p1+1):
        for jl_1 in range(0, p1+1):
            for il_2 in range(0, p2+1):
                for jl_2 in range(0, p2+1):
                    for il_3 in range(0, p3+1):
                        for jl_3 in range(0, p3+1):

                            v = 0.0
                            for g1 in range(0, k1):
                                for g2 in range(0, k2):
                                    for g3 in range(0, k3):
                                        Ni = bs1[il_1, 0, g1] * bs2[il_2, 0, g2] * bs3[il_3, 0, g3]
                                        Ni_x = bs1[il_1, 1, g1] * bs2[il_2, 0, g2] * bs3[il_3, 0, g3]
                                        Ni_y = bs1[il_1, 0, g1] * bs2[il_2, 1, g2] * bs3[il_3, 0, g3]
                                        Ni_z = bs1[il_1, 0, g1] * bs2[il_2, 0, g2] * bs3[il_3, 1, g3]

                                        Nj = bs1[jl_1, 0, g1] * bs2[jl_2, 0, g2] * bs3[jl_3, 0, g3]
                                        Nj_x = bs1[jl_1, 1, g1] * bs2[jl_2, 0, g2] * bs3[jl_3, 0, g3]
                                        Nj_y = bs1[jl_1, 0, g1] * bs2[jl_2, 1, g2] * bs3[jl_3, 0, g3]
                                        Nj_z = bs1[jl_1, 0, g1] * bs2[jl_2, 0, g2] * bs3[jl_3, 1, g3]

                                        x    = u1[g1]
                                        y    = u2[g2]
                                        z    = u3[g3]
                                        wvol = w1[g1] * w2[g2] * w3[g3]

                                        v += ({__WEAK_FORM__}) * wvol
                            mat[il_1, il_2, il_3, p1 + jl_1 - il_1, p2 + jl_2 - il_2, p3 + jl_3 - il_3] = v
"""

template_header_3d_scalar = '#$ header procedure {__KERNEL_NAME__}(int, int, int, int, int, int, double [:,:,:], double [:,:,:], double [:,:,:], double [:], double [:], double [:], double [:], double [:], double [:], double [:,:,:,:,:,:] {__TYPES__})'
# .............................................

# .............................................
#          KERNEL     2D case - block
# .............................................
template_1d_block = """
def {__KERNEL_NAME__}(p1, k1, basis, u, w, {__MAT_ARGS__} {__ARGS__}):
{__MAT_INIT__}
    for il_1 in range(0, p1+1):
        for jl_1 in range(0, p1+1):

{__ACCUM_INIT__}
            for g1 in range(0, k1):
                Ni = basis[il_1, 0, g1]
                Ni_x = basis[il_1, 1, g1]

                Nj = basis[jl_1, 0, g1]
                Nj_x = basis[jl_1, 1, g1]

                x    = u[g1]
                wvol = w[g1]

{__ACCUM__}
{__ACCUM_ASSIGN__}
"""

template_header_1d_block = '#$ header procedure {__KERNEL_NAME__}(int, int, double [:,:,:], double [:], double [:], {__MAT_TYPES__} {__TYPES__})'
# .............................................

# .............................................
#          KERNEL     2D case - block
# .............................................
template_2d_block = """
def {__KERNEL_NAME__}(p1, p2, k1, k2, bs1, bs2, u1, u2, w1, w2, {__MAT_ARGS__} {__ARGS__}):
{__MAT_INIT__}
    for il_1 in range(0, p1+1):
        for jl_1 in range(0, p1+1):
            for il_2 in range(0, p2+1):
                for jl_2 in range(0, p2+1):

{__ACCUM_INIT__}
                    for g1 in range(0, k1):
                        for g2 in range(0, k2):
                            Ni = bs1[il_1, 0, g1] * bs2[il_2, 0, g2]
                            Ni_x = bs1[il_1, 1, g1] * bs2[il_2, 0, g2]
                            Ni_y = bs1[il_1, 0, g1] * bs2[il_2, 1, g2]

                            Nj = bs1[jl_1, 0, g1] * bs2[jl_2, 0, g2]
                            Nj_x = bs1[jl_1, 1, g1] * bs2[jl_2, 0, g2]
                            Nj_y = bs1[jl_1, 0, g1] * bs2[jl_2, 1, g2]

                            x    = u1[g1]
                            y    = u2[g2]
                            wvol = w1[g1] * w2[g2]

{__ACCUM__}
{__ACCUM_ASSIGN__}
"""

template_header_2d_block = '#$ header procedure {__KERNEL_NAME__}(int, int, int, int, double [:,:,:], double [:,:,:], double [:], double [:], double [:], double [:], {__MAT_TYPES__} {__TYPES__})'
# .............................................

# .............................................
#          KERNEL     3D case - block
# .............................................
template_3d_block = """
def {__KERNEL_NAME__}(p1, p2, p3, k1, k2, k3, bs1, bs2, bs3, u1, u2, u3, w1, w2, w3, {__MAT_ARGS__} {__ARGS__}):
{__MAT_INIT__}
    for il_1 in range(0, p1+1):
        for jl_1 in range(0, p1+1):
            for il_2 in range(0, p2+1):
                for jl_2 in range(0, p2+1):
                    for il_3 in range(0, p3+1):
                        for jl_3 in range(0, p3+1):

{__ACCUM_INIT__}
                            for g1 in range(0, k1):
                                for g2 in range(0, k2):
                                    for g3 in range(0, k3):
                                        Ni = bs1[il_1, 0, g1] * bs2[il_2, 0, g2] * bs3[il_3, 0, g3]
                                        Ni_x = bs1[il_1, 1, g1] * bs2[il_2, 0, g2] * bs3[il_3, 0, g3]
                                        Ni_y = bs1[il_1, 0, g1] * bs2[il_2, 1, g2] * bs3[il_3, 0, g3]
                                        Ni_z = bs1[il_1, 0, g1] * bs2[il_2, 0, g2] * bs3[il_3, 1, g3]

                                        Nj = bs1[jl_1, 0, g1] * bs2[jl_2, 0, g2] * bs3[jl_3, 0, g3]
                                        Nj_x = bs1[jl_1, 1, g1] * bs2[jl_2, 0, g2] * bs3[jl_3, 0, g3]
                                        Nj_y = bs1[jl_1, 0, g1] * bs2[jl_2, 1, g2] * bs3[jl_3, 0, g3]
                                        Nj_z = bs1[jl_1, 0, g1] * bs2[jl_2, 0, g2] * bs3[jl_3, 1, g3]

                                        x    = u1[g1]
                                        y    = u2[g2]
                                        z    = u3[g3]
                                        wvol = w1[g1] * w2[g2] * w3[g3]

{__ACCUM__}
{__ACCUM_ASSIGN__}
"""

template_header_3d_block = '#$ header procedure {__KERNEL_NAME__}(int, int, int, int, int, int, double [:,:,:], double [:,:,:], double [:,:,:], double [:], double [:], double [:], double [:], double [:], double [:], {__MAT_TYPES__} {__TYPES__})'
# .............................................

# .............................................
#          SYMBOL     1D case - scalar
# .............................................
symbol_1d_scalar ="""
def {__SYMBOL_NAME__}(arr_x, arr_t1, symbol {__ARGS__}):
    from numpy import sin
    from numpy import cos

    n1 = len(arr_x)
    symbol[:] = 0.
    for i1 in range(0, n1):
        x = arr_x[i1]
        t1 = arr_t1[i1]

        symbol[i1] = {__SYMBOL_EXPR__}
"""

symbol_header_1d_scalar = '#$ header procedure {__SYMBOL_NAME__}(double [:], double [:], double [:] {__TYPES__})'
# .............................................

# .............................................
#          SYMBOL     2D case - scalar
# .............................................
symbol_2d_scalar ="""
def {__SYMBOL_NAME__}({__SYMBOL_ARGS__}, expr {__ARGS__}):
    from numpy import sin
    from numpy import cos

    expr = {__SYMBOL_EXPR__}
    return expr
"""

symbol_header_2d_scalar = '#$ header procedure {__SYMBOL_NAME__}(double, double, double, double, double {__TYPES__})'
# .............................................

# .............................................
#          SYMBOL     3D case - scalar
# .............................................
symbol_3d_scalar ="""
def {__SYMBOL_NAME__}({__SYMBOL_ARGS__}, expr {__ARGS__}):
    from numpy import sin
    from numpy import cos

    expr = {__SYMBOL_EXPR__}
    return expr
"""

symbol_header_3d_scalar = '#$ header procedure {__SYMBOL_NAME__}(double, double, double, double, double {__TYPES__})'
# .............................................
