# coding: utf-8

# TODO: - eval field as a function
#       - gather fields by space and then print eval on each of them
#       - eval for block field
#       -

from collections import OrderedDict

from gelato.core import get_index_derivatives_atom

def construct_field_names(fields):
    fields_str = ', '.join(i.name for i in fields)
    fields_str = ', {}'.format(fields_str)

    return fields_str

def construct_field_coeffs_names(fields):
    field_coeffs = OrderedDict()
    for f in fields:
        coeffs = 'coeff_{}'.format(f.name)
        field_coeffs[str(f.name)] = coeffs
    return field_coeffs

def print_field_coeffs(field_coeffs):
    ls = [v for v in list(field_coeffs.values())]
    field_coeffs_str = ', '.join(i for i in ls)

    # add ',' for kernel signature
    field_coeffs_str = ', {}'.format(field_coeffs_str)
    return field_coeffs_str

def print_field_coeffs_slices(fields, dim):
    # ...
    if dim == 1:
        slices = 'is1-test_p1:is1+1'

    elif dim == 2:
        slices = 'is1-test_p1:is1+1,is2-test_p2:is2+1'

    elif dim == 3:
        slices = 'is1-test_p1:is1+1,is2-test_p2:is2+1,is3-test_p3:is3+1'

    coeffs = []
    for F in fields:
        coeff_str = '{field}.coeffs[{slices}]'.format(field=F.name,
                                                      slices=slices)
        coeffs.append(coeff_str)
    # ...

    fields_coeffs_str = ', '.join(c for c in coeffs)
    fields_coeffs_str = ', {}'.format(fields_coeffs_str)

    return fields_coeffs_str

def print_field_coeffs_types(field_coeffs, dim):
    field_types = []
    slices = ','.join(':' for i in range(0, dim))
    for v in list(field_coeffs.values()):
        field_types.append('double [{slices}]'.format(slices=slices))

    field_types_str = ', '.join(i for i in field_types)
    field_types_str = ', {}'.format(field_types_str)
    return field_types_str

def construct_field_values_names(expr, fields, verbose=False):
    field_values = OrderedDict()
    for F in fields:
        indices = get_index_derivatives_atom(expr, F)

        ls = []

        # TODO remove this
        #      right now we need to add it manually
        ls.append(F.name)

        # ... change code to a string like '_x'
        for index in indices:
            code = ''
            for k,n in list(index.items()):
                code += k*n
            name = '{name}_{code}'.format(name=F.name, code=code)
            ls.append(name)
        # ...

        field_values[F.name] = ls

    # ...
    if verbose:
        print('>>> field_values = ', field_values)
    # ...

    return field_values

def print_eval_field(expr, dim, fields, field_coeffs, field_values,
                     verbose=False):
    """."""
    # ... identation (def function body)
    tab = ' '*4
    # ...

    # ... field values init
    sizes = ','.join('k{}'.format(i) for i in range(1, dim+1))
    if dim > 1:
        sizes = '({})'.format(sizes)

    lines = []
    for k, fs in list(field_values.items()):
        for f in fs:
            line = '{field}_values = zeros({sizes})'.format(field=f, sizes=sizes)
            line = tab + line

            lines.append(line)

    field_init_str = '\n'.join(line for line in lines)
    # ...

    # ... update identation to be inside the loop
    for i in range(0, 2*dim):
        tab += ' '*4

    tab_base = tab
    # ...

    # ...
    if dim == 1:
        e_pattern = '{field}{deriv}_values[g1] += {coeff}[jl_1]*Nj{deriv}'
    elif dim == 2:
        e_pattern = '{field}{deriv}_values[g1,g2] += {coeff}[jl_1,jl_2]*Nj{deriv}'
    elif dim ==3:
        e_pattern = '{field}{deriv}_values[g1,g2,g3] += {coeff}[jl_1,jl_2,jl_3]*Nj{deriv}'
    else:
        raise NotImplementedError('only 1d, 2d and 3d are available')

    lines = []
    for k, fs in list(field_values.items()):
        coeff = field_coeffs[k]
        for f in fs:
            ls = f.split('_')
            if len(ls) == 1:
                deriv = ''
            else:
                deriv = '_{}'.format(ls[-1])
            line = e_pattern.format(field=k, deriv=deriv, coeff=coeff)
            line = tab + line

            lines.append(line)

    field_accum_str = '\n'.join(line for line in lines)
    # ...

    # ...
    # TODO
    pattern = 'scalar'
    # ...

    # ...
    template_str = 'eval_field_{dim}d_{pattern}'.format(dim=dim, pattern=pattern)
    try:
        from .templates import eval_field_1d_scalar
        from .templates import eval_field_2d_scalar
        from .templates import eval_field_3d_scalar

        template = eval(template_str)
    except:
        raise ValueError('Could not find the corresponding template {}'.format(template_str))
    # ...

    # ...
    field_coeffs_str = ', '.join('{}'.format(c) for c in list(field_coeffs.values()))

    code = template.format(__FIELD_INIT__=field_init_str,
                           __FIELD_ACCUM__=field_accum_str)
    # ...

    return code


def print_assign_field(expr, dim, fields, field_coeffs, field_values, tab,
                       verbose=False):
    """."""
    # ...
    if verbose:
        print('> input     := {0}'.format(expr))
    # ...

    # ...
    if dim == 1:
        e_pattern = '{field}{deriv} = {field}{deriv}_values[g1]'
    elif dim == 2:
        e_pattern = '{field}{deriv} = {field}{deriv}_values[g1,g2]'
    elif dim ==3:
        e_pattern = '{field}{deriv} = {field}{deriv}_values[g1,g2,g3]'
    else:
        raise NotImplementedError('only 1d, 2d and 3d are available')

    lines = []
    for k, fs in list(field_values.items()):
        coeff = field_coeffs[k]
        for f in fs:
            ls = f.split('_')
            if len(ls) == 1:
                deriv = ''
            else:
                deriv = '_{}'.format(ls[-1])
            line = e_pattern.format(field=k, deriv=deriv)
            line = tab + line

            lines.append(line)

    field_value_str = '\n'.join(line for line in lines)
    # ...

    return field_value_str
