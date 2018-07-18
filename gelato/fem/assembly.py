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

from .matrix import construct_element_matrix_names
from .matrix import print_element_matrix_args
from .matrix import print_element_matrix_decs
from .matrix import construct_global_matrix_names
from .matrix import print_global_matrix_args
from .matrix import print_global_matrix_decs
from .matrix import print_global_matrix_update

from .vector import construct_element_vector_names
from .vector import print_element_vector_args
from .vector import print_element_vector_decs
from .vector import construct_global_vector_names
from .vector import print_global_vector_args
from .vector import print_global_vector_decs
from .vector import print_global_vector_update

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

    if is_bilinear_form:
        form = 'bilinear'
    elif is_linear_form:
        form = 'linear'
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
    element_mat_args = ''
    element_mat_args_str = ''
    element_mat_decs_str = ''

    global_mat_args = ''
    global_mat_args_str = ''
    global_mat_decs_str = ''
    global_mat_update_str = ''

    element_vec_args = ''
    element_vec_args_str = ''
    element_vec_decs_str = ''

    global_vec_args = ''
    global_vec_args_str = ''
    global_vec_decs_str = ''
    global_vec_update_str = ''

    if is_bilinear_form:
        element_mat_args = construct_element_matrix_names(n_rows, n_cols)
        element_mat_args_str = print_element_matrix_args(n_rows, n_cols, element_mat_args)
        element_mat_decs_str = print_element_matrix_decs(n_rows, n_cols, dim, element_mat_args, tab)

        global_mat_args = construct_global_matrix_names(n_rows, n_cols)
        global_mat_args_str = print_global_matrix_args(n_rows, n_cols, global_mat_args)
        global_mat_decs_str = print_global_matrix_decs(n_rows, n_cols, global_mat_args, tab)

        # ...
        for i in range(0, dim):
            tab += ' '*4

        global_mat_update_str = print_global_matrix_update(n_rows, n_cols, dim,
                                                           element_mat_args,
                                                           global_mat_args,
                                                           tab)
        # ...

    elif is_linear_form:
        element_vec_args = construct_element_vector_names(n_rows)
        element_vec_args_str = print_element_vector_args(n_rows, element_vec_args)
        element_vec_decs_str = print_element_vector_decs(n_rows, dim, element_vec_args, tab)

        global_vec_args = construct_global_vector_names(n_rows)
        global_vec_args_str = print_global_vector_args(n_rows, global_vec_args)
        global_vec_decs_str = print_global_vector_decs(n_rows, global_vec_args, tab)

        # ...
        for i in range(0, dim):
            tab += ' '*4

        global_vec_update_str = print_global_vector_update(n_rows, dim,
                                                           element_vec_args,
                                                           global_vec_args,
                                                           tab)
        # ...
    # ...

    # ...
    code = template.format(__ASSEMBLY_NAME__=assembly_name,
                           __ARGS__=args,
                           __DOCSTRING__=docstring,
                           __FIELDS__=fields_str,
                           __FIELDS_COEFFS__=fields_coeffs_str,
                           __GLOBAL_MAT_DEC__=global_mat_decs_str,
                           __GLOBAL_MAT_ARGS__=global_mat_args_str,
                           __ELEMENT_MAT_DEC__=element_mat_decs_str,
                           __ELEMENT_MAT_ARGS__=element_mat_args_str,
                           __GLOBAL_MAT_UPDATE__=global_mat_update_str,
                           __GLOBAL_VEC_DEC__=global_vec_decs_str,
                           __GLOBAL_VEC_ARGS__=global_vec_args_str,
                           __ELEMENT_VEC_DEC__=element_vec_decs_str,
                           __ELEMENT_VEC_ARGS__=element_vec_args_str,
                           __GLOBAL_VEC_UPDATE__=global_vec_update_str,
                           __KERNEL_NAME__=kernel_name)
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
    assembly = namespace[name]
    # ...

    # ... export the python code of the module
    if export_pyfile:
        write_code(name, code, ext='py', folder='.pyccel')
    # ...

    return assembly
