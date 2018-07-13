# coding: utf-8

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
