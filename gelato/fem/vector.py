# coding: utf-8

from .utils import _convert_int_to_float


def construct_element_vector_names(n_rows):
    vec_args = []
    for i in range(0, n_rows):
        vec = 'vec_{i}'.format(i=i)
        vec_args.append(vec)

    return vec_args

def print_element_vector_args(n_rows, vec_args):
    vec_args_str = ', '.join(vec for vec in vec_args)
    return vec_args_str

def print_element_vector_init(n_rows, vec_args, size, tab):
    slices = ','.join(':' for i in range(0, size))

    lines = []
    for i in range(0, n_rows):
        vec = vec_args[i]

        line = '{vec}[{slices}] = 0.0'.format(vec=vec,slices=slices)
        line = tab + line

        lines.append(line)

    vec_init_str = '\n'.join(line for line in lines)
    return vec_init_str

def print_accumulation_var_init(n_rows, tab):
    lines = []
    for i in range(0, n_rows):
        line = 'v_{i} = 0.0'.format(i=i,)
        line = tab + line

        lines.append(line)

    accum_init_str = '\n'.join(line for line in lines)
    return accum_init_str

def print_accumulation_var(n_rows, expr, tab):
    lines = []
    for i in range(0, n_rows):
        line = 'v_{i} += ({__WEAK_FORM__}) * wvol'
        e = _convert_int_to_float(expr[i].evalf())
        # we call evalf to avoid having fortran doing the evaluation of rational
        # division
        line = line.format(i=i, __WEAK_FORM__=e)
        line = tab + line

        lines.append(line)

    accum_str = '\n'.join(line for line in lines)
    return accum_str

def print_bilinear_accumulation_assign(n_rows, dim, tab):
    if dim == 1:
        e_pattern = 'vec_{i}[il_1] = v_{i}'

    elif dim == 2:
        e_pattern = 'vec_{i}[il_1, il_2] = v_{i}'

    elif dim ==3:
        e_pattern = 'vec_{i}[il_1, il_2, il_3] = v_{i}'

    else:
        raise NotImplementedError('only 1d, 2d and 3d are available')

    lines = []
    for i in range(0, n_rows):
        line = e_pattern.format(i=i)
        line = tab + line

        lines.append(line)

    accum_assign_str = '\n'.join(line for line in lines)
    return accum_assign_str

def print_linear_accumulation_assign(n_rows, dim, tab):
    if dim == 1:
        e_pattern = 'vec_{i}[il_1] = v_{i}'

    elif dim == 2:
        e_pattern = 'vec_{i}[il_1, il_2] = v_{i}'

    elif dim ==3:
        e_pattern = 'vec_{i}[il_1, il_2, il_3] = v_{i}'

    else:
        raise NotImplementedError('only 1d, 2d and 3d are available')

    lines = []
    for i in range(0, n_rows):
        line = e_pattern.format(i=i)
        line = tab + line

        lines.append(line)

    accum_assign_str = '\n'.join(line for line in lines)
    return accum_assign_str

# TODO use space degrees
def print_element_vector_decs(n_rows, dim, vec_args, tab):
    if dim == 1:
        pattern = 'zeros( test_p1+1 )'

    elif dim == 2:
        pattern = 'zeros( (test_p1+1, test_p2+1) )'

    elif dim == 3:
        pattern = 'zeros( (test_p1+1, test_p2+1, test_p3+1) )'

    else:
        raise NotImplementedError('only 1d, 2d and 3d are available')

    lines = []
    for i in range(0, n_rows):
        vec = vec_args[i]
        line = '{vec} = {pattern}'.format(vec=vec, pattern=pattern)
        line = tab + line
        lines.append(line)

    decs = '\n'.join(i for i in lines)
    return decs

def construct_global_vector_names(n_rows):
    vec_args = []
    for i in range(0, n_rows):
        vec = 'rhs_{i}'.format(i=i)
        vec_args.append(vec)

    return vec_args

def print_global_vector_args(n_rows, vec_args):
    vec_args_str = ', '.join(vec for vec in vec_args)
    return vec_args_str

# TODO use spaces
def print_global_vector_decs(n_rows, vec_args, tab):
    pattern = 'StencilVector( test_space.vector_space )'

    lines = []
    for i in range(0, n_rows):
        vec = vec_args[i]
        line = '{vec} = {pattern}'.format(vec=vec, pattern=pattern)
        line = tab + line
        lines.append(line)

    decs = '\n'.join(i for i in lines)
    return decs

# TODO use spaces
def print_global_vector_update(n_rows, dim,
                               element_vec_args,
                               global_vec_args, tab):

    lslices = ','.join(':' for i in range(0, dim))

    lines = []
    for i in range(0, n_rows):
        lvec = element_vec_args[i]
        gvec = global_vec_args[i]

        # ... every vector should have access to its corresponding space and
        # degrees
        if dim == 1:
            gslices = 'is1-test_p1:is1+1'

        elif dim == 2:
            gslices = 'is1-test_p1:is1+1, is2-test_p2:is2+1'

        elif dim == 3:
            gslices = 'is1-test_p1:is1+1, is2-test_p2:is2+1, is3-test_p3:is3+1'
        # ...

        line = '{gvec}[{gslices}] += {lvec}[{lslices}]'.format(lvec=lvec,
                                                               gvec=gvec,
                                                               lslices=lslices,
                                                               gslices=gslices)
        line = tab + line
        lines.append(line)

    decs = '\n'.join(i for i in lines)
    return decs

