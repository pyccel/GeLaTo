# -*- coding: utf-8 -*-

import importlib

from symfe.core import BilinearForm
from symfe.codegen import arguments_datatypes_as_dict
from symfe.codegen import arguments_datatypes_split
from symfe.codegen.utils import _convert_int_to_float
from symfe.codegen.utils import write_code

from gelato.core import gelatize

from .utils import (print_position_args, print_fourier_args, print_mat_args)

def compile_symbol(name, a,
                   degrees,
                   n_elements=None,
                   verbose=False,
                   namespace=globals(),
                   context=None,
                   backend='python',
                   export_pyfile=True):
    """."""
    if not isinstance(a, BilinearForm):
           raise TypeError('Expecting a BilinearForm')

    # TODO: nderiv must be computed from the weak form
    nderiv = 1

    # ... weak form attributs
    dim = a.ldim
    fields = a.fields

    # TODO improve
    is_block = False
    is_vector = False

    if verbose:
        print('> dim    = ', dim)
        print('> Fields = ', fields)
    # ...

    # ... contants
    d_args = arguments_datatypes_as_dict(a.constants)
    args, dtypes = arguments_datatypes_split(d_args)
    # ...

    # TODO check what are the free_symbols of expr,
    #      to make sure the final code will compile
    #      the remaining free symbols must be the trial/test basis functions,
    #      and the coordinates

    # ...
    if is_vector and not( is_block ):
        raise NotImplementedError('We only treat the case of a block space, for '
                                  'which all components have are identical.')
    # ...

    # ...
    if is_block:
        pattern = 'block'
    elif is_vector:
        raise NotImplementedError('TODO.')
    else:
        pattern = 'scalar'
    # ...

    # ... get name of the template to be used
    template_str = '_symbol_{dim}d_{pattern}'.format(dim=dim, pattern=pattern)
    # ...

    # ... import the variable from the templates module
    #     NOTE: THE PATH IS HARD CODED HERE!
    try:
        package = importlib.import_module("gelato.codegen.templates.symbol")

    except:
        raise ImportError('could not import {0}'.format(name))

    template = getattr(package, template_str)
    # ...

    # ...
    expr = gelatize(a, degrees=degrees, n_elements=n_elements)
    # ...

    # ... TODO add an attribute called symbol_expr to the BilinearForm?
#    setattr(a, 'symbol_expr', expr)
    # ...

    # ... identation (def function body)
    tab_base = ' '*4
    tab = tab_base
    # ...

    # ... append n_elements as argument of the generated symbol function
    if n_elements:
        n_elements_str = ''
        n_elements_types_str = ''

    else:
        ns = ['nx', 'ny', 'nz'][:dim]
        n_elements_str = ', '.join(n for n in ns)
        n_elements_str = ', {}'.format(n_elements_str)

        n_elements_types_str = ', '.join('int' for n in ns)
        n_elements_types_str = ', {}'.format(n_elements_types_str)
    # ...

    # ... field coeffs
    if fields:
        raise NotImplementedError()

        field_coeffs = construct_field_coeffs_names(fields)
        field_coeffs_str = print_field_coeffs(field_coeffs)
        field_types_str = print_field_coeffs_types(field_coeffs, dim)

        field_values = construct_field_values_names(expr, fields)

        eval_field_str = print_eval_field(expr, dim,
                                          fields,
                                          field_coeffs,
                                          field_values,
                                          verbose=verbose)


        # ... update identation to be inside the loop
        if is_bilinear_form:
            for i in range(0, 3*dim):
                tab += ' '*4

        elif is_linear_form:
            for i in range(0, 2*dim):
                tab += ' '*4

        elif is_function_form:
            for i in range(0, dim):
                tab += ' '*4

        field_value_str = print_assign_field(expr, dim,
                                             fields,
                                             field_coeffs,
                                             field_values,
                                             tab, verbose=verbose)

        tab = tab_base
        # ...

    else:
        field_coeffs_str = ''
        eval_field_str   = ''
        field_value_str  = ''
        field_types_str  = ''
    # ...

    # ...
    if fields:
        raise NotImplementedError()

        # we call normalize a second time, and activate the modification of the
        # fields, this will turn terms like dx(F) into F_x if F is a field
        expr = normalize(expr, enable_fields=True)
    # ...

    # ...
    x_args_str = print_position_args(dim)
    t_args_str = print_fourier_args(dim)
    # ...

    # ... TODO be careful of conflict when adding block case
    mat_args_str = print_mat_args()
    # ...

    # ... compute indentation
    # TODO
#    tab += ' '*4
    # ...

    # ...
    tab = tab_base
    # ...

    # ...
    if is_block:
        raise NotImplementedError('')

        # ... - initializing element matrices
        #     - define arguments
        # test functions and trial functions
        if is_bilinear_form:
            size = 2*dim

        # test functions
        elif is_linear_form:
            size = dim

        elif is_function_form:
            raise NotImplementedError('')

        n_rows = test_n_components
        n_cols = trial_n_components
        mat_args = construct_element_matrix_names(n_rows, n_cols)
        mat_args_str = print_element_matrix_args(n_rows, n_cols, mat_args)
        mat_init_str = print_element_matrix_init(n_rows, n_cols, mat_args, size, tab)

        # ... update identation to be inside the loop
        for i in range(0, size):
            tab += ' '*4

        tab_base = tab
        # ...

        # ... initializing accumulation variables
        accum_init_str = print_accumulation_var_init(n_rows, n_cols, tab)
        # ...

        # .. update indentation
        for i in range(0, dim):
            tab += ' '*4
        # ...

        # ... accumulation contributions
        accum_str = print_accumulation_var(n_rows, n_cols, expr, tab)
        # ...

        # ... assign accumulated values to element matrix
        if is_bilinear_form:
            accum_assign_str = print_bilinear_accumulation_assign(n_rows,
                                                                  n_cols, dim,
                                                                  tab_base)

        elif is_linear_form:
            accum_assign_str = print_linear_accumulation_assign(n_rows,
                                                                n_cols, dim,
                                                                tab_base)
        # ...

        code = template.format(__KERNEL_NAME__=name,
                               __X_ARGS__=x_args_str,
                               __T_ARGS__=t_args_str,
                               __MAT_ARGS__=mat_args_str,
                               __N_ELEMENTS__=n_elements_str,
                               __FIELD_COEFFS__=field_coeffs_str,
                               __FIELD_EVALUATION__=eval_field_str,
                               __MAT_INIT__=mat_init_str,
                               __ACCUM_INIT__=accum_init_str,
                               __FIELD_VALUE__=field_value_str,
                               __TEST_FUNCTION__=test_function_str,
                               __TRIAL_FUNCTION__=trial_function_str,
                               __ACCUM__=accum_str,
                               __ACCUM_ASSIGN__=accum_assign_str,
                               __ARGS__=args)

    else:
        # we call evalf to avoid having fortran doing the evaluation of rational
        # division

        e = _convert_int_to_float(expr.evalf())
        code = template.format(__SYMBOL_NAME__=name,
                               __SYMBOL_EXPR__=e.evalf(),
                               __X_ARGS__=x_args_str,
                               __T_ARGS__=t_args_str,
                               __MAT_ARGS__=mat_args_str,
                               __N_ELEMENTS__=n_elements_str,
                               __FIELD_COEFFS__=field_coeffs_str,
                               __FIELD_EVALUATION__=eval_field_str,
                               __FIELD_VALUE__=field_value_str,
                               __ARGS__=args)
    # ...

#    print('--------------')
#    print(code)
#    print('--------------')
#    import sys; sys.exit(0)

    # ...
    if context:
        from pyccel.epyccel import ContextPyccel

        if isinstance(context, ContextPyccel):
            context = [context]
        elif isinstance(context, (list, tuple)):
            for i in context:
                assert(isinstance(i, ContextPyccel))
        else:
            raise TypeError('Expecting a ContextPyccel or list/tuple of ContextPyccel')

        # append functions to the namespace
        for c in context:
            for k,v in list(c.functions.items()):
                namespace[k] = v[0]
    # ...

    # ...
    exec(code, namespace)
    kernel = namespace[name]
    # ...

    # ... export the python code of the module
    if export_pyfile:
        write_code(name, code, ext='py', folder='.pyccel')
    # ...

    return kernel
