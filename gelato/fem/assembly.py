# coding: utf-8

# TODO: - allow for giving a name for the trial/test basis
#       - Ni/Nj should be Ni_0/Nj_0
#       - define templates as proper python functions
#       - use redbaron to modify the template

#     NOTE: THE PATH OF TEMPLATES IS HARD CODED!


from sympy.core.containers import Tuple
from sympy import Matrix
from sympy import Integer, Float

from numbers import Number
from collections import OrderedDict

from numpy import unique
import os
import importlib

from gelato.core import gelatize
from gelato.core import BilinearForm, LinearForm, FunctionForm
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

from .matrix import construct_element_matrix_names
from .matrix import print_element_matrix_args
from .matrix import print_element_matrix_decs
from .matrix import construct_global_matrix_names
from .matrix import print_global_matrix_args
from .matrix import print_global_matrix_update
from .matrix import construct_argument_matrix_name
from .matrix import print_argument_matrix_kwargs
from .matrix import print_define_global_matrix

from .vector import construct_element_vector_names
from .vector import print_element_vector_args
from .vector import print_element_vector_decs
from .vector import construct_global_vector_names
from .vector import print_global_vector_args
from .vector import print_global_vector_decs
from .vector import print_global_vector_update
from .vector import construct_argument_vector_name
from .vector import print_argument_vector_kwargs
from .vector import print_define_global_vector

from .array import construct_element_array_names
from .array import print_element_array_args
from .array import print_element_array_decs
from .array import construct_global_array_names
from .array import print_global_array_args
from .array import print_global_array_decs
from .array import print_global_array_update
from .array import construct_argument_array_name
from .array import print_argument_array_kwargs
from .array import print_define_global_array

# NOTE: test_n_components, trial_n_components  will be provided after calling compile_kernel
def compile_assembly(name, a, kernel_name=None,
                     verbose=False,
                     namespace=globals(),
                     context=None,
                     is_vector=False,
                     test_n_components=1,
                     trial_n_components=1,
                     backend='python',
                     export_pyfile=True):
    """."""
    # ...
    assembly_name = name
    docstring     = ''
    if kernel_name is None:
        kernel_name = 'kernel_{}'.format(name)
    # ...

    # ... weak form attributs
    dim = a.ldim
    fields = a.fields
    is_bilinear_form = isinstance(a, BilinearForm)
    is_linear_form = isinstance(a, LinearForm)
    is_function_form = isinstance(a, FunctionForm)

    if is_bilinear_form:
        form = 'bilinear'

    elif is_linear_form:
        form = 'linear'

    elif is_function_form:
        form = 'function'
    # ...

    # ... contants
    d_args = arguments_datatypes_as_dict(a.constants)
    args, dtypes = arguments_datatypes_split(d_args)
    # ...

    # ... fields
    if fields:
        fields_str = ', '.join(i.name for i in fields)
        fields_str = ', {}'.format(fields_str)

        # ...
        # TODO must use span, to use local index in kernel
        slices = None
        if dim == 1:
            slices = 's1:s1+test_p1+1'

        elif dim == 2:
            slices = 's1:s1+test_p1+1,s2:s2+test_p2+1'

        elif dim == 3:
            slices = 's1:s1+test_p1+1,s2:s2+test_p2+1,s3:s3+test_p3+1'

        coeffs = []
        for F in fields:
            coeff_str = '{field}.coeffs[{slices}]'.format(field=F.name,
                                                          slices=slices)
            coeffs.append(coeff_str)
        # ...

        fields_coeffs_str = ', '.join(c for c in coeffs)
        fields_coeffs_str = ', {}'.format(fields_coeffs_str)

    else:
        fields_str = ''
        fields_coeffs_str = ''
    # ...

    # ... get name of the template to be used
    template_str = '_assembly_{form}_{dim}d'.format(dim=dim,
                                                    form=form)
    # ...

    # ... import the variable from the templates module
    #     NOTE: THE PATH IS HARD CODED HERE!
    try:
        package = importlib.import_module("gelato.fem.templates.assembly")

    except:
        raise ImportError('could not import {0}'.format(name))

    template = getattr(package, template_str)
    # ...

    # ... identation (def function body)
    tab = ' '*4
    # ...

    # ...
    n_rows = test_n_components
    n_cols = trial_n_components
    # ...

    # ...
    argument_mat = ''
    argument_mat_kwargs = ''

    global_mat_args = ''
    global_mat_args_str = ''
    global_mat_decs_str = ''
    global_mat_update_str = ''

    element_mat_args = ''
    element_mat_args_str = ''
    element_mat_decs_str = ''

    argument_vec = ''
    argument_vec_kwargs = ''

    global_vec_args = ''
    global_vec_args_str = ''
    global_vec_decs_str = ''
    global_vec_update_str = ''

    element_vec_args = ''
    element_vec_args_str = ''
    element_vec_decs_str = ''

    argument_arr = ''
    argument_arr_kwargs = ''

    global_arr_args = ''
    global_arr_args_str = ''
    global_arr_decs_str = ''
    global_arr_update_str = ''

    element_arr_args = ''
    element_arr_args_str = ''
    element_arr_decs_str = ''

    element_wise_str = ''

    if is_bilinear_form:
        argument_mat = construct_argument_matrix_name(n_rows, n_cols)
        argument_mat_kwargs = print_argument_matrix_kwargs(argument_mat)

        global_mat_args = construct_global_matrix_names(n_rows, n_cols)
        global_mat_args_str = print_global_matrix_args(n_rows, n_cols, global_mat_args)
        global_mat_decs_str = print_define_global_matrix(n_rows, n_cols, global_mat_args, argument_mat, tab)

        element_mat_args = construct_element_matrix_names(n_rows, n_cols)
        element_mat_args_str = print_element_matrix_args(n_rows, n_cols, element_mat_args)
        element_mat_decs_str = print_element_matrix_decs(n_rows, n_cols, dim, element_mat_args, tab)

        # ...
        for i in range(0, dim):
            tab += ' '*4

        global_mat_update_str = print_global_matrix_update(n_rows, n_cols, dim,
                                                           element_mat_args,
                                                           global_mat_args,
                                                           tab)
        # ...

    elif is_linear_form:
        argument_vec = construct_argument_vector_name(n_rows)
        argument_vec_kwargs = print_argument_vector_kwargs(argument_vec)

        global_vec_args = construct_global_vector_names(n_rows)
        global_vec_args_str = print_global_vector_args(n_rows, global_vec_args)
        global_vec_decs_str = print_define_global_vector(n_rows, global_vec_args, argument_vec, tab)

        element_vec_args = construct_element_vector_names(n_rows)
        element_vec_args_str = print_element_vector_args(n_rows, element_vec_args)
        element_vec_decs_str = print_element_vector_decs(n_rows, dim, element_vec_args, tab)

        # ...
        for i in range(0, dim):
            tab += ' '*4

        global_vec_update_str = print_global_vector_update(n_rows, dim,
                                                           element_vec_args,
                                                           global_vec_args,
                                                           tab)
        # ...

    elif is_function_form:
        argument_arr = construct_argument_array_name(n_rows)
        argument_arr_kwargs = print_argument_array_kwargs(argument_arr)

        global_arr_args = construct_global_array_names(n_rows)
        global_arr_args_str = print_global_array_args(n_rows, global_arr_args)
        global_arr_decs_str = print_define_global_array(n_rows, dim, global_arr_args, argument_arr, tab)

        element_arr_args = construct_element_array_names(n_rows)
        element_arr_args_str = print_element_array_args(n_rows, element_arr_args)
        element_arr_decs_str = print_element_array_decs(n_rows, dim, element_arr_args, tab)

        # ...
        for i in range(0, dim):
            tab += ' '*4

        global_arr_update_str = print_global_array_update(n_rows, dim,
                                                           element_arr_args,
                                                           global_arr_args,
                                                           tab)
        # ...

        element_wise_str = ', element_wise=False'
        # TODO improve, not sure this can be pyccelized!!
        argument_arr = '{arg} if element_wise else sum({arg})'.format(arg=argument_arr)
    # ...

    # ...
    code = template.format(__ASSEMBLY_NAME__=assembly_name,
                           __ARGS__=args,
                           __DOCSTRING__=docstring,
                           __FIELDS__=fields_str,
                           __FIELDS_COEFFS__=fields_coeffs_str,

                           __ARGUMENT_MAT_KWARGS__=argument_mat_kwargs,
                           __GLOBAL_MAT_DEC__=global_mat_decs_str,
                           __GLOBAL_MAT_ARGS__=argument_mat,
                           __ELEMENT_MAT_DEC__=element_mat_decs_str,
                           __ELEMENT_MAT_ARGS__=element_mat_args_str,
                           __GLOBAL_MAT_UPDATE__=global_mat_update_str,

                           __ARGUMENT_VEC_KWARGS__=argument_vec_kwargs,
                           __GLOBAL_VEC_DEC__=global_vec_decs_str,
                           __GLOBAL_VEC_ARGS__=argument_vec,
                           __ELEMENT_VEC_DEC__=element_vec_decs_str,
                           __ELEMENT_VEC_ARGS__=element_vec_args_str,
                           __GLOBAL_VEC_UPDATE__=global_vec_update_str,

                           __ARGUMENT_ARR_KWARGS__=argument_arr_kwargs,
                           __GLOBAL_ARR_DEC__=global_arr_decs_str,
                           __GLOBAL_ARR_ARGS__=argument_arr,
                           __ELEMENT_ARR_DEC__=element_arr_decs_str,
                           __ELEMENT_ARR_ARGS__=element_arr_args_str,
                           __GLOBAL_ARR_UPDATE__=global_arr_update_str,

                           __ELEMENT_WISE__=element_wise_str,

                           __KERNEL_NAME__=kernel_name)
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
    assembly = namespace[name]
    # ...

    # ... export the python code of the module
    if export_pyfile:
        write_code(name, code, ext='py', folder='.pyccel')
    # ...

    return assembly
