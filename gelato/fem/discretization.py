# coding: utf-8

# ... TODO: - add args


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

from .kernel import compile_kernel
from .assembly import compile_assembly

import types

# TODO change target to meta var, and update spaces_str
_template ="""
def {__NAME__}( target{__ARGS__}{__FIELDS__}{__ELEMENT_WISE__} ):
    {__DOCSTRING__}
    return {__ASSEMBLY_NAME__}( target, {__SPACE_ARGS__}{__ARGS__}{__FIELDS__}{__ELEMENT_WISE_KWARG__} )
"""

_template_docstring = """
\"\"\"
Assembly method for {__KIND_FORM__}.

This method is calling two functions that have been automatically generated:
    - {__ASSEMBLY_NAME__} low-level assembly function
    - {__KERNEL_NAME__}   low-level kernel function
{__PARAMETERS__}
\"\"\"
"""

_docstring_header = """
Parameters
----------
"""

_pattern_docstring_argument = """
{__ARG__} : {__TYPE__}
   {__LABEL__}
"""

# TODO add element_wise to docstring
def docstring_arguments(constants, d_args, is_function_form=False):
    if len(constants) == 0 and not is_function_form:
        return ''

    pattern = _pattern_docstring_argument

    lines = []

    # ... constants
    for c in constants:
        dtype = d_args[c.name]
        arg = pattern.format(__ARG__=c.name, __TYPE__=dtype, __LABEL__=c.label)
        lines += [arg]
    # ...

    # ... element wise for FunctionForm
    if is_function_form:
        label = 'Assemble FunctionForm on every element if True. [Default: False]'

        arg = pattern.format(__ARG__='element_wise',
                             __TYPE__='bool',
                             __LABEL__=label)
        lines += [arg]
    # ...

    code = '\n'.join(line for line in lines)

    txt = '{header}{arguments}'.format(header=_docstring_header,
                                       arguments=code)
    return txt

# TODO add check on spaces
# TODO pass root to compile kernel and assembly
def discretize(a, spaces,
               verbose=False,
               namespace=globals(),
               context=None,
               name=None,
               root=None,
               backend='python',
               export_pyfile=True):
    """."""
    # ...
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

    # ...
    if root is None:
        root = '.pyccel'
    # ...

    # ...
    if name is None:
        name = '{form}_{hash}'.format(form=form, hash=abs(hash(a)))

    assembly_name = 'assembly_{}'.format(name)
    kernel_name   = 'kernel_{}'.format(name)
    # ...

    # ...
    kernel = compile_kernel(kernel_name, a,
                            verbose=verbose,
                            namespace=namespace,
                            context=context,
                            backend=backend,
                            export_pyfile=export_pyfile)
    # ...

    # ...
    assembly = compile_assembly(assembly_name, a, kernel_name,
                                test_n_components=a.n_rows,
                                trial_n_components=a.n_cols,
                                verbose=verbose,
                                namespace=namespace,
                                context=context,
                                backend=backend,
                                export_pyfile=export_pyfile)
    # ...

    # ...
    if is_bilinear_form:
        spaces_str = 'target.discrete_spaces[0], target.discrete_spaces[1]'

    elif is_linear_form:
        spaces_str = 'target.discrete_spaces'

    elif is_function_form:
        spaces_str = 'target.discrete_spaces'
    # ...

    # ... contants
    d_args = arguments_datatypes_as_dict(a.constants)
    args, dtypes = arguments_datatypes_split(d_args)
    # ...

    # ... fields
    fields_str = ''
    if fields:
        fields_str = ', '.join(i.name for i in fields)
        fields_str = ', {}'.format(fields_str)
    # ...

    # ...
    args_docstring = docstring_arguments(a.constants, d_args,
                                         is_function_form=is_function_form)
    # ...

    # ...
    element_wise_str = ''
    element_wise_kwarg_str = ''
    if is_function_form:
        element_wise_str = ', element_wise=False'
        element_wise_kwarg_str = ', element_wise=element_wise'
    # ...

    # ...
    docstring = _template_docstring.format(__KIND_FORM__=form,
                                           __ASSEMBLY_NAME__=assembly_name,
                                           __KERNEL_NAME__=kernel_name,
                                           __PARAMETERS__=args_docstring)
    # identation (def function body)
    tab = ' '*4

    lines = []
    for line in docstring.split('\n'):
        lines.append(tab + line)
    docstring = '\n'.join(i for i in lines)
    # ...

    # ...
    code = _template.format(__NAME__=name,
                            __ASSEMBLY_NAME__=assembly_name,
                            __SPACE_ARGS__=spaces_str,
                            __ARGS__=args,
                            __FIELDS__=fields_str,
                            __ELEMENT_WISE__=element_wise_str,
                            __ELEMENT_WISE_KWARG__=element_wise_kwarg_str,
                            __DOCSTRING__=docstring)
    # ...

    # ...
    exec(code, namespace)
    _assemble = namespace[name]
    # ...

    # ...
    setattr(a, 'discrete_spaces', spaces)
    a.assemble = types.MethodType(_assemble, a)
    # ...
