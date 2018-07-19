# coding: utf-8

from .utils import _convert_int_to_float

def _element_matrix_name(i,j):
    return 'mat_{i}{j}'.format(i=i,j=j)

def _global_matrix_name(i,j):
    return 'M_{i}{j}'.format(i=i,j=j)

def construct_element_matrix_names(n_rows, n_cols):
    mat_args = []
    for i in range(0, n_rows):
        ls = []
        for j in range(0, n_cols):
            mat = _element_matrix_name(i,j)
            ls.append(mat)
        mat_args.append(ls)

    return mat_args

def print_element_matrix_args(n_rows, n_cols, mat_args):
    ls = []
    for i in range(0, n_rows):
        for j in range(0, n_cols):
            ls.append(mat_args[i][j])
    mat_args_str = ', '.join(mat for mat in ls)
    return mat_args_str

def print_element_matrix_init(n_rows, n_cols, mat_args, size, tab):
    slices = ','.join(':' for i in range(0, size))

    lines = []
    for i in range(0, n_rows):
        for j in range(0, n_cols):
            mat = mat_args[i][j]

            line = '{mat}[{slices}] = 0.0'.format(mat=mat,slices=slices)
            line = tab + line

            lines.append(line)

    mat_init_str = '\n'.join(line for line in lines)
    return mat_init_str

def print_accumulation_var_init(n_rows, n_cols, tab):
    lines = []
    for i in range(0, n_rows):
        for j in range(0, n_cols):
            line = 'v_{i}{j} = 0.0'.format(i=i,j=j)
            line = tab + line

            lines.append(line)

    accum_init_str = '\n'.join(line for line in lines)
    return accum_init_str

def print_accumulation_var(n_rows, n_cols, expr, tab):
    lines = []
    for i in range(0, n_rows):
        for j in range(0, n_cols):
            line = 'v_{i}{j} += ({__WEAK_FORM__}) * wvol'
            e = _convert_int_to_float(expr[i,j].evalf())
            # we call evalf to avoid having fortran doing the evaluation of rational
            # division
            line = line.format(i=i, j=j, __WEAK_FORM__=e)
            line = tab + line

            lines.append(line)

    accum_str = '\n'.join(line for line in lines)
    return accum_str

def print_bilinear_accumulation_assign(n_rows, n_cols, dim, tab):
    if dim == 1:
        e_pattern = 'mat_{i}{j}[il_1, test_p1 + jl_1 - il_1] = v_{i}{j}'

    elif dim == 2:
        e_pattern = 'mat_{i}{j}[il_1, il_2, test_p1 + jl_1 - il_1, test_p2 + jl_2 - il_2] = v_{i}{j}'

    elif dim ==3:
        e_pattern = 'mat_{i}{j}[il_1, il_2, il_3, test_p1 + jl_1 - il_1, test_p2 + jl_2 - il_2, test_p3 + jl_3 - il_3] = v_{i}{j}'

    else:
        raise NotImplementedError('only 1d, 2d and 3d are available')

    lines = []
    for i in range(0, n_rows):
        for j in range(0, n_cols):
            line = e_pattern.format(i=i,j=j)
            line = tab + line

            lines.append(line)

    accum_assign_str = '\n'.join(line for line in lines)
    return accum_assign_str

def print_linear_accumulation_assign(n_rows, n_cols, dim, tab):
    if dim == 1:
        e_pattern = 'mat_{i}{j}[il_1] = v_{i}{j}'

    elif dim == 2:
        e_pattern = 'mat_{i}{j}[il_1, il_2] = v_{i}{j}'

    elif dim ==3:
        e_pattern = 'mat_{i}{j}[il_1, il_2, il_3] = v_{i}{j}'

    else:
        raise NotImplementedError('only 1d, 2d and 3d are available')

    lines = []
    for i in range(0, n_rows):
        for j in range(0, n_cols):
            line = e_pattern.format(i=i,j=j)
            line = tab + line

            lines.append(line)

    accum_assign_str = '\n'.join(line for line in lines)
    return accum_assign_str

# TODO use space degrees
def print_element_matrix_decs(n_rows, n_cols, dim, mat_args, tab):
    if dim == 1:
        pattern = 'zeros( (test_p1+1, 2*trial_p1+1) )'

    elif dim == 2:
        pattern = 'zeros( (test_p1+1, test_p2+1, 2*trial_p1+1, 2*trial_p2+1) )'

    elif dim == 3:
        pattern = 'zeros( (test_p1+1, test_p2+1, test_p3+1, 2*trial_p1+1, 2*trial_p2+1, 2*trial_p3+1) )'

    else:
        raise NotImplementedError('only 1d, 2d and 3d are available')

    lines = []
    for i in range(0, n_rows):
        for j in range(0, n_cols):
            mat = mat_args[i][j]
            line = '{mat} = {pattern}'.format(mat=mat, pattern=pattern)
            line = tab + line
            lines.append(line)

    decs = '\n'.join(i for i in lines)
    return decs

def construct_global_matrix_names(n_rows, n_cols):
    mat_args = []
    for i in range(0, n_rows):
        ls = []
        for j in range(0, n_cols):
            mat = _global_matrix_name(i,j)
            ls.append(mat)
        mat_args.append(ls)

    return mat_args

def print_global_matrix_args(n_rows, n_cols, mat_args):
    ls = []
    for i in range(0, n_rows):
        for j in range(0, n_cols):
            ls.append(mat_args[i][j])
    mat_args_str = ', '.join(mat for mat in ls)
    return mat_args_str

# TODO use spaces
def print_global_matrix_decs(n_rows, n_cols, mat_args):
    pattern = 'StencilMatrix( test_space.vector_space, trial_space.vector_space )'

    lines = []
    for i in range(0, n_rows):
        for j in range(0, n_cols):
            mat = mat_args[i][j]
            line = '{mat} = {pattern}'.format(mat=mat, pattern=pattern)
            lines.append(line)

    decs = '\n'.join(i for i in lines)
    return decs

# TODO use spaces
def print_global_matrix_update(n_rows, n_cols, dim,
                               element_mat_args,
                               global_mat_args, tab):

    lslices = ','.join(':' for i in range(0, 2*dim))
    suffix = ','.join(':' for i in range(0, dim))

    lines = []
    for i in range(0, n_rows):
        for j in range(0, n_cols):
            lmat = element_mat_args[i][j]
            gmat = global_mat_args[i][j]

            # ... every matrix should have access to its corresponding space and
            # degrees
            if dim == 1:
                gslices = 'is1-test_p1:is1+1'

            elif dim == 2:
                gslices = 'is1-test_p1:is1+1, is2-test_p2:is2+1'

            elif dim == 3:
                gslices = 'is1-test_p1:is1+1, is2-test_p2:is2+1, is3-test_p3:is3+1'

            gslices = '{gslices}, {suffix}'.format(gslices=gslices,
                                                   suffix=suffix)
            # ...

            line = '{gmat}[{gslices}] += {lmat}[{lslices}]'.format(lmat=lmat,
                                                                   gmat=gmat,
                                                                   lslices=lslices,
                                                                   gslices=gslices)
            line = tab + line
            lines.append(line)

    decs = '\n'.join(i for i in lines)
    return decs

def construct_argument_matrix_name(n_rows, n_cols):
    if (n_rows == 1) and (n_cols == 1):
        return _global_matrix_name(0,0)
    else:
        return 'd_matrix'

def print_argument_matrix_kwargs(argument_mat):
    return ', {}=None'.format(argument_mat)

def print_import_stencil_matrix():
    return 'from spl.linalg.stencil import StencilMatrix'

_template_define_global_matrix = """
if {__ARGUMENT_MAT__} is None:
    {__IMPORT_STENCIL__}
    {__DECS__}
{__ELSE__}
"""
# TODO add comment to the generated code
def print_define_global_matrix(n_rows, n_cols, global_mat_args, argument_mat, tab):
    global_mat_decs_str = print_global_matrix_decs(n_rows, n_cols, global_mat_args)

    pattern = _template_define_global_matrix

    import_str = print_import_stencil_matrix()

    # TODO in the case of block space
    else_str = ''

    code = pattern.format(__ARGUMENT_MAT__=argument_mat,
                          __IMPORT_STENCIL__=import_str,
                          __DECS__=global_mat_decs_str,
                          __ELSE__=else_str)

    lines = []
    for line in code.split('\n'):
        line = tab + line
        lines.append(line)

    code = '\n'.join(line for line in lines)

    return code
