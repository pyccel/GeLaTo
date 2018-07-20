# coding: utf-8

from .utils import _convert_int_to_float

def _element_array_name(i):
    return 'norm_{i}'.format(i=i)

def _global_array_name(i):
    return 'rhs_{i}'.format(i=i)

def construct_element_array_names(n_rows):
    norm_args = []
    for i in range(0, n_rows):
        arr = _element_array_name(i)
        norm_args.append(arr)

    return norm_args

def print_element_array_args(n_rows, norm_args):
    norm_args_str = ', '.join(arr for arr in norm_args)
    return norm_args_str

def print_element_array_init(n_rows, norm_args, size, tab):
    slices = ','.join(':' for i in range(0, size))

    lines = []
    for i in range(0, n_rows):
        arr = norm_args[i]

        line = '{arr}[{slices}] = 0.0'.format(arr=arr,slices=slices)
        line = tab + line

        lines.append(line)

    norm_init_str = '\n'.join(line for line in lines)
    return norm_init_str

def print_accumulation_var_init(n_rows, tab):
    lines = []
    for i in range(0, n_rows):
        line = 'v_{i} = 0.0'.format(i=i)
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
        e_pattern = 'norm_{i}[il_1] = v_{i}'

    elif dim == 2:
        e_pattern = 'norm_{i}[il_1, il_2] = v_{i}'

    elif dim ==3:
        e_pattern = 'norm_{i}[il_1, il_2, il_3] = v_{i}'

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
        e_pattern = 'norm_{i}[il_1] = v_{i}'

    elif dim == 2:
        e_pattern = 'norm_{i}[il_1, il_2] = v_{i}'

    elif dim ==3:
        e_pattern = 'norm_{i}[il_1, il_2, il_3] = v_{i}'

    else:
        raise NotImplementedError('only 1d, 2d and 3d are available')

    lines = []
    for i in range(0, n_rows):
        line = e_pattern.format(i=i)
        line = tab + line

        lines.append(line)

    accum_assign_str = '\n'.join(line for line in lines)
    return accum_assign_str

# TODO use size from function
def print_element_array_decs(n_rows, dim, norm_args, tab):
    if dim == 1:
        pattern = '0.'

    elif dim == 2:
        pattern = '0.'

    elif dim == 3:
        pattern = '0.'

    else:
        raise NotImplementedError('only 1d, 2d and 3d are available')

    lines = []
    for i in range(0, n_rows):
        arr = norm_args[i]
        line = '{arr} = {pattern}'.format(arr=arr, pattern=pattern)
        line = tab + line
        lines.append(line)

    decs = '\n'.join(i for i in lines)
    return decs

def construct_global_array_names(n_rows):
    norm_args = []
    for i in range(0, n_rows):
        arr = _global_array_name(i)
        norm_args.append(arr)

    return norm_args

def print_global_array_args(n_rows, norm_args):
    norm_args_str = ', '.join(arr for arr in norm_args)
    return norm_args_str

# TODO use spaces
def print_global_array_decs(n_rows, dim, norm_args):
    if dim == 1:
        pattern = 'zeros( n_cells1 )'

    elif dim == 2:
        pattern = 'zeros( (n_cells1, n_cells2) )'

    elif dim == 3:
        pattern = 'zeros( (n_cells1, n_cells2, n_cells3) )'

    lines = []
    for i in range(0, n_rows):
        arr = norm_args[i]
        line = '{arr} = {pattern}'.format(arr=arr, pattern=pattern)
        lines.append(line)

    decs = '\n'.join(i for i in lines)
    return decs

# TODO use spaces
def print_global_array_update(n_rows, dim,
                               element_norm_args,
                               global_norm_args, tab):

    lines = []
    for i in range(0, n_rows):
        larr = element_norm_args[i]
        garr = global_norm_args[i]

        # ... every array should have access to its corresponding space and
        # degrees
        if dim == 1:
            gslices = 'ie1'

        elif dim == 2:
            gslices = 'ie1, ie2'

        elif dim == 3:
            gslices = 'ie1, ie2, ie3'
        # ...

        line = '{garr}[{gslices}] += {larr}'.format(larr=larr,
                                                    garr=garr,
                                                    gslices=gslices)
        line = tab + line
        lines.append(line)

    decs = '\n'.join(i for i in lines)
    return decs

def construct_argument_array_name(n_rows):
    if (n_rows == 1):
        return _global_array_name(0)
    else:
        return 'd_array'

def print_argument_array_kwargs(argument_mat):
    return ', {}=None'.format(argument_mat)

# TODO
def print_import_array():
    return 'from numpy import zeros'

def print_set_dict_array(n_rows, argument_arr, norm_args):
    if n_rows == 1:
        return ''

    lines = [argument_arr + ' = {}']
    for i in range(0, n_rows):
        M = _global_array_name(i)
        line = '{d}[{i}] = {M}'.format(d=argument_arr, i=i, M=M)
        lines.append(line)

    return '\n'.join(i for i in lines)

def print_get_dict_array(n_rows, argument_arr, norm_args):
    if n_rows == 1:
        return ''

    lines = []
    for i in range(0, n_rows):
        M = _global_array_name(i)
        line = '{M} = {d}[{i}]'.format(d=argument_arr, i=i, M=M)
        lines.append(line)

    return '\n'.join(i for i in lines)

_template_define_global_array = """
if {__ARGUMENT_VEC__} is None:
{__IMPORT_STENCIL__}
{__DECS__}
{__SET_DICT__}
{__ELSE__}
{__GET_DICT__}
"""
# TODO add comment to the generated code
def print_define_global_array(n_rows, dim, global_norm_args, argument_arr, tab):
    # ...
    def _indent_block(txt):
        indent = ' '*4

        lines = []
        for line in txt.split('\n'):
            line = indent + line
            lines.append(line)

        return '\n'.join(line for line in lines)
    # ...

    # ...
    global_norm_decs_str = print_global_array_decs(n_rows, dim, global_norm_args)
    global_norm_decs_str = _indent_block( global_norm_decs_str )
    # ...

    # ...
    import_str = print_import_array()
    import_str = _indent_block( import_str )
    # ...

    # ...
    set_dict_str = print_set_dict_array(n_rows, argument_arr, global_norm_args)
    set_dict_str = _indent_block( set_dict_str )
    # ...

    # ...
    get_dict_str = print_get_dict_array(n_rows, argument_arr, global_norm_args)
    get_dict_str = _indent_block( get_dict_str )
    # ...

    # ...
    if n_rows == 1:
        else_str = ''
    else:
        else_str = 'else:'
    # ...

    pattern = _template_define_global_array
    code = pattern.format(__ARGUMENT_VEC__=argument_arr,
                          __IMPORT_STENCIL__=import_str,
                          __DECS__=global_norm_decs_str,
                          __SET_DICT__=set_dict_str,
                          __GET_DICT__=get_dict_str,
                          __ELSE__=else_str)

    lines = []
    for line in code.split('\n'):
        line = tab + line
        lines.append(line)

    code = '\n'.join(line for line in lines)

    return code
