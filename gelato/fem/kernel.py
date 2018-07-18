# coding: utf-8

# TODO: - allow for giving a name for the trial/test basis
#       - Ni/Nj should be Ni_0/Nj_0
#       - define templates as proper python functions
#       - use redbaron to modify the template
#       - check what are the free_symbols of expr,

#     NOTE: THE PATH OF TEMPLATES IS HARD CODED!


from sympy.core.containers import Tuple
from sympy import Matrix
from sympy import Integer, Float

from numbers import Number
from collections import OrderedDict

from numpy import unique
import os
import importlib

from gelato.core import gelatize, normalize
from gelato.core import BilinearForm, LinearForm
from gelato.core import Constant
from gelato.core import Field

from .utils import _is_base_function
from .utils import _convert_int_to_float
from .utils import _count_letter
from .utils import construct_test_functions
from .utils import construct_trial_functions
from .utils import mkdir_p
from .utils import write_code
from .utils import arguments_datatypes_as_dict
from .utils import arguments_datatypes_split

from .field import construct_field_coeffs_names
from .field import print_field_coeffs
from .field import print_field_coeffs_types
from .field import construct_field_values_names
from .field import print_eval_field
from .field import print_assign_field

from .matrix import construct_element_matrix_names
from .matrix import print_element_matrix_args
from .matrix import print_element_matrix_init
from .matrix import print_accumulation_var_init
from .matrix import print_accumulation_var
from .matrix import print_bilinear_accumulation_assign
from .matrix import print_linear_accumulation_assign


def compile_kernel(name, a,
                   verbose=False,
                   namespace=globals(),
                   context=None,
                   backend='python',
                   export_pyfile=True):
    """."""
    if not isinstance(a, (BilinearForm, LinearForm)):
           raise TypeError('Expecting BilinearForm, LinearForm')

    # TODO: nderiv must be computed from the weak form
    nderiv = 1

    # ... weak form attributs
    dim = a.ldim
    fields = a.fields
    is_bilinear_form = isinstance(a, BilinearForm)
    is_linear_form = isinstance(a, LinearForm)

    if verbose:
        print('> dim    = ', dim)
        print('> Fields = ', fields)
    # ...

    # ... TODO improve
    TEST_BASIS = 'Ni'
    TRIAL_BASIS = 'Nj'

    is_test_function = lambda a: _is_base_function(a, TEST_BASIS)
    is_trial_function = lambda a: _is_base_function(a, TRIAL_BASIS)

    d_basis = None
    if is_bilinear_form:
        U = a.test_spaces[0]
        V = a.trial_spaces[0]
        d_basis = {V: TRIAL_BASIS, U: TEST_BASIS}

    elif is_linear_form:
        U = a.test_spaces[0]
        d_basis = {U: TEST_BASIS}

    else:
        raise NotImplementedError('Only Bilinear and Linear forms are available')

    expr = gelatize(a, basis=d_basis)

    if verbose:
        print('> gelatized  >>> {0}'.format(expr))

    is_block = isinstance(expr, Matrix)
    test_n_components = None
    trial_n_components = None
    if is_block:
        if is_bilinear_form:
            test_n_components = expr.shape[0]
            trial_n_components = expr.shape[1]

        elif is_linear_form:
            test_n_components = expr.shape[0]
            trial_n_components = 1 # since we are using Matrix from sympy

        else:
            raise NotImplementedError('Only Bilinear and Linear forms are available')

    # TODO add once we handle feec
    is_vector = False

    # we add attributs of Bilinear/Linear form
    setattr(a, 'is_block', is_block)
    setattr(a, 'n_rows', test_n_components)
    setattr(a, 'n_cols', trial_n_components)
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
    if is_bilinear_form:
        template_str = '_bilinear_form_{dim}d_{pattern}'.format(dim=dim, pattern=pattern)
    elif is_linear_form:
        template_str = '_linear_form_{dim}d_{pattern}'.format(dim=dim, pattern=pattern)
    else:
        raise NotImplementedError('Only Bilinear and Linear forms are available')
    # ...

    # ... import the variable from the templates module
    #     NOTE: THE PATH IS HARD CODED HERE!
    try:
        if is_bilinear_form:
            package = importlib.import_module("gelato.fem.templates.bilinear_form")
        elif is_linear_form:
            package = importlib.import_module("gelato.fem.templates.linear_form")
        else:
            raise ValueError('only linear and bilinear form are available')

    except:
        raise ImportError('could not import {0}'.format(name))

    template = getattr(package, template_str)
    # ...

    # ... identation (def function body)
    tab = ' '*4
    # ...

    # ... field coeffs
    if fields:
        field_coeffs = construct_field_coeffs_names(fields)
        field_coeffs_str = print_field_coeffs(field_coeffs)
        field_types_str = print_field_coeffs_types(field_coeffs, dim)

        field_values = construct_field_values_names(expr, fields)

        eval_field_str = print_eval_field(expr, dim,
                                          fields,
                                          field_coeffs,
                                          field_values,
                                          verbose=verbose)
        field_value_str = print_assign_field(expr, dim,
                                             fields,
                                             field_coeffs,
                                             field_values,
                                             tab, verbose=verbose)

    else:
        field_coeffs_str = ''
        eval_field_str   = ''
        field_value_str  = ''
        field_types_str  = ''
    # ...

    # ...
    if fields:
        # we call normalize a second time, and activate the modification of the
        # fields, this will turn terms like dx(F) into F_x if F is a field
        expr = normalize(expr, enable_fields=True)
    # ...

    # ... compute indentation
    tab_base = tab
    if is_bilinear_form:
        # dim loops for test functions
        # and dim loops for trial functions
        for i in range(0, 2*dim):
            tab += ' '*4

    elif is_linear_form:
        # dim loops for test functions
        for i in range(0, dim):
            tab += ' '*4

    # dim loops for quadrature points
    for i in range(0, dim):
        tab += ' '*4
    # ...

    # ... TODO improve following the same idea as for Fields
    test_function_str = ''
    if is_bilinear_form or is_linear_form:
        # ... print test functions
        d_test_basis = construct_test_functions(nderiv, dim)
        test_names = [i.name for i in expr.free_symbols if is_test_function(i)]
        test_names.sort()

        lines = []
        for arg in test_names:
            if arg == TEST_BASIS:
                basis = ' * '.join(d_test_basis[k,0] for k in range(1, dim+1))
                line = '{test} = {basis}'.format(test=TEST_BASIS, basis=basis)
            else:
                deriv = arg.split('_')[-1]
                nx = _count_letter(deriv, 'x')
                ny = _count_letter(deriv, 'y')
                nz = _count_letter(deriv, 'z')
                basis = ' * '.join(d_test_basis[k,d] for k,d in zip(range(1, dim+1), [nx,ny,nz]))
                line = '{test}_{deriv} = {basis}'.format(test=TEST_BASIS, deriv=deriv, basis=basis)
            lines.append(tab+line)
        test_function_str = '\n'.join(l for l in lines)
        # ...
    # ...

    # ... TODO improve following the same idea as for Fields
    trial_function_str = ''
    if is_bilinear_form:
        # ... print trial functions
        d_trial_basis = construct_trial_functions(nderiv, dim)
        trial_names = [i.name for i in expr.free_symbols if is_trial_function(i)]
        trial_names.sort()

        lines = []
        for arg in trial_names:
            if arg == TRIAL_BASIS:
                basis = ' * '.join(d_trial_basis[k,0] for k in range(1, dim+1))
                line = '{trial} = {basis}'.format(trial=TRIAL_BASIS, basis=basis)
            else:
                deriv = arg.split('_')[-1]
                nx = _count_letter(deriv, 'x')
                ny = _count_letter(deriv, 'y')
                nz = _count_letter(deriv, 'z')
                basis = ' * '.join(d_trial_basis[k,d] for k,d in zip(range(1, dim+1), [nx,ny,nz]))
                line = '{trial}_{deriv} = {basis}'.format(trial=TRIAL_BASIS, deriv=deriv, basis=basis)
            lines.append(tab+line)
        trial_function_str = '\n'.join(l for l in lines)
    # ...

    # ...
    tab = tab_base
    # ...

    # ...
    if is_block:
        # ... - initializing element matrices
        #     - define arguments
        # test functions and trial functions
        if is_bilinear_form:
            size = 2*dim
        # test functions
        elif is_linear_form:
            size = dim

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
                               __MAT_ARGS__=mat_args_str,
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
        e = _convert_int_to_float(expr.evalf())
        # we call evalf to avoid having fortran doing the evaluation of rational
        # division
        code = template.format(__KERNEL_NAME__=name,
                               __FIELD_COEFFS__=field_coeffs_str,
                               __FIELD_EVALUATION__=eval_field_str,
                               __FIELD_VALUE__=field_value_str,
                               __TEST_FUNCTION__=test_function_str,
                               __TRIAL_FUNCTION__=trial_function_str,
                               __WEAK_FORM__=e,
                               __ARGS__=args)
    # ...

#    print('--------------')
#    print(code)
#    print('--------------')

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
