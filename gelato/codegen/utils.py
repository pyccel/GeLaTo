# -*- coding: utf-8 -*-

def _matrix_name(i,j):
    return 'M_{i}{j}'.format(i=i,j=j)

def construct_x_args_names(dim):
    pattern = 'arr_x{}'
    ls = []
    for i in range(1, dim+1):
        x = pattern.format(i)
        ls.append(x)
    return ls

def construct_t_args_names(dim):
    pattern = 'arr_t{}'
    ls = []
    for i in range(1, dim+1):
        x = pattern.format(i)
        ls.append(x)
    return ls

def print_position_args(x_args):
    ls = ', '.join(i for i in x_args)
    return ls

def print_fourier_args(t_args):
    ls = ', '.join(i for i in t_args)
    return ', {}'.format(ls)

# TODO improve in the case of block
def print_mat_args():
    return ', mat'

_docstring_header = """
Parameters
----------
"""

_pattern_docstring_argument = """
{__ARG__} : {__TYPE__}
   {__LABEL__}
"""

def docstring_arguments(constants, d_args):
    if len(constants) == 0:
        return ''

    pattern = _pattern_docstring_argument

    lines = []

    # ... constants
    for c in constants:
        dtype = d_args[c.name]
        arg = pattern.format(__ARG__=c.name, __TYPE__=dtype, __LABEL__=c.label)
        lines += [arg]
    # ...

    code = '\n'.join(line for line in lines)

    txt = '{header}{arguments}'.format(header=_docstring_header,
                                       arguments=code)
    return txt




def print_import_zeros():
    return 'from numpy import zeros'

def print_matrix_decs(n_rows, n_cols, x_args, mat_args):
    dim = len(x_args)
    sizes = ['len({})'.format(i) for i in x_args]

    pattern = 'zeros( {lpar}{size}{rpar} )'
    if dim == 1:
        pattern = pattern.format(size=sizes[0], lpar='', rpar='')
    else:
        sizes = ', '.join(str(i) for i in sizes)
        pattern = pattern.format(size=sizes, lpar='(', rpar=')')

    lines = []
    for i in range(0, n_rows):
        for j in range(0, n_cols):
            mat = mat_args[i][j]
            line = '{mat} = {pattern}'.format(mat=mat, pattern=pattern)
            lines.append(line)

    decs = '\n'.join(i for i in lines)
    return decs

def print_set_dict_matrix(n_rows, n_cols, argument_mat, mat_args):
    if (n_rows == 1) and (n_cols == 1):
        return ''

    lines = [argument_mat + ' = {}']
    for i in range(0, n_rows):
        for j in range(0, n_cols):
            M = _matrix_name(i,j)
            line = '{d}[{i},{j}] = {M}'.format(d=argument_mat, i=i, j=j, M=M)
            lines.append(line)

    return '\n'.join(i for i in lines)

def print_get_dict_matrix(n_rows, n_cols, argument_mat, mat_args):
    if (n_rows == 1) and (n_cols == 1):
        return ''

    lines = []
    for i in range(0, n_rows):
        for j in range(0, n_cols):
            M = _matrix_name(i,j)
            line = '{M} = {d}[{i},{j}]'.format(d=argument_mat, i=i, j=j, M=M)
            lines.append(line)

    return '\n'.join(i for i in lines)


def construct_argument_matrix_name(n_rows, n_cols):
    if (n_rows == 1) and (n_cols == 1):
        return _matrix_name(0,0)
    else:
        return 'd_matrix'

def construct_argument_matrix_name(n_rows, n_cols):
    if (n_rows == 1) and (n_cols == 1):
        return _matrix_name(0,0)
    else:
        return 'd_matrix'

def print_argument_matrix_kwargs(argument_mat):
    return ', {}=None'.format(argument_mat)

def construct_matrix_names(n_rows, n_cols):
    mat_args = []
    for i in range(0, n_rows):
        ls = []
        for j in range(0, n_cols):
            mat = _matrix_name(i,j)
            ls.append(mat)
        mat_args.append(ls)

    return mat_args

def print_matrix_args(n_rows, n_cols, mat_args):
    ls = []
    for i in range(0, n_rows):
        for j in range(0, n_cols):
            ls.append(mat_args[i][j])
    mat_args_str = ', '.join(mat for mat in ls)
    return mat_args_str
#    return ', {}'.format(mat_args_str)


_template_define_matrix = """
if {__ARGUMENT_MAT__} is None:
{__IMPORT_ZEROS__}
{__DECS__}
{__SET_DICT__}
{__ELSE__}
{__GET_DICT__}
"""
# TODO add comment to the generated code
def print_define_matrix(n_rows, n_cols, x_args, mat_args, argument_mat, tab):
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
    mat_decs_str = print_matrix_decs(n_rows, n_cols, x_args, mat_args)
    mat_decs_str = _indent_block( mat_decs_str )
    # ...

    # ...
    import_str = print_import_zeros()
    import_str = _indent_block( import_str )
    # ...

    # ...
    set_dict_str = print_set_dict_matrix(n_rows, n_cols, argument_mat, mat_args)
    set_dict_str = _indent_block( set_dict_str )
    # ...

    # ...
    get_dict_str = print_get_dict_matrix(n_rows, n_cols, argument_mat, mat_args)
    get_dict_str = _indent_block( get_dict_str )
    # ...

    # ...
    if (n_rows == 1) and (n_cols == 1):
        else_str = ''
    else:
        else_str = 'else:'
    # ...

    pattern = _template_define_matrix
    code = pattern.format(__ARGUMENT_MAT__=argument_mat,
                          __IMPORT_ZEROS__=import_str,
                          __DECS__=mat_decs_str,
                          __SET_DICT__=set_dict_str,
                          __GET_DICT__=get_dict_str,
                          __ELSE__=else_str)

    lines = []
    for line in code.split('\n'):
        line = tab + line
        lines.append(line)

    code = '\n'.join(line for line in lines)

    return code
