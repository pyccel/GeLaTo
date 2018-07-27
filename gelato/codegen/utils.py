# -*- coding: utf-8 -*-

def print_position_args(dim):
    pattern = 'arr_x{}'
    ls = []
    for i in range(1, dim+1):
        x = pattern.format(i)
        ls.append(x)
    ls = ', '.join(i for i in ls)
    return ls

def print_fourier_args(dim):
    pattern = 'arr_t{}'
    ls = []
    for i in range(1, dim+1):
        x = pattern.format(i)
        ls.append(x)
    ls = ', '.join(i for i in ls)
    return ', {}'.format(ls)

# TODO improve in the case of block
def print_mat_args():
    return ', mat'

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
