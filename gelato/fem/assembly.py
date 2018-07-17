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

from .kernel import compile_kernel


def compile_assembly(name, a, kernel_name=None,
                     verbose=False,
                     namespace=globals(),
                     context=None,
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

    # ... TODO is_block must be set inside compile_kernel?
#    if is_block:
#        pattern = 'block'
#    elif is_vector:
#        raise NotImplementedError('TODO.')
#    else:
#        pattern = 'scalar'
    pattern = 'scalar'
    # ...

    # ... get name of the template to be used
    template_str = '_assembly_{form}_{dim}d_{pattern}'.format(dim=dim,
                                                              pattern=pattern,
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
    code = template.format(__ASSEMBLY_NAME__=assembly_name,
                           __ARGS__=args,
                           __DOCSTRING__=docstring,
                           __FIELDS__=fields_str,
                           __FIELDS_COEFFS__=fields_coeffs_str,
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
