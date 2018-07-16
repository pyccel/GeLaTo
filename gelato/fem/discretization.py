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

from .kernel import compile_kernel
from .assembly import compile_assembly

import types

# TODO change target to meta var, and update spaces_str
_template ="""
def {__NAME__}(target):
    {__DOCSTRING__}
    return {__ASSEMBLY_NAME__}(target, {__SPACE_ARGS__})
"""

_template_docstring = """
\"\"\"
Assembly method for {__KIND_FORM__}.

This method is calling two functions that have been automatically generated:
    {__ASSEMBLY_NAME__} low-level assembly function
    {__KERNEL_NAME__}   low-level kernel function
\"\"\"
"""


# TODO add check on spaces
# TODO pass root to compile kernel and assembly
def discretize(a, spaces,
               d_constants={},
               d_args={},
               verbose=False,
               namespace=globals(),
               context=None,
               name=None,
               root=None,
               backend='python',
               export_pyfile=True):
    """."""
    # ...
    is_bilinear_form = isinstance(a, BilinearForm)
    is_linear_form = isinstance(a, LinearForm)

    if is_bilinear_form:
        form = 'bilinear'
    elif is_linear_form:
        form = 'linear'
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
                            spaces=spaces,
                            d_constants=d_constants,
                            d_args=d_args,
                            verbose=verbose,
                            namespace=namespace,
                            context=context,
                            backend=backend,
                            export_pyfile=export_pyfile)
    # ...

    # ...
    assembly = compile_assembly(assembly_name, a, kernel_name,
                                spaces=spaces,
                                d_constants=d_constants,
                                d_args=d_args,
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
    # ...

    # ...
    docstring = _template_docstring.format(__KIND_FORM__=form,
                                           __ASSEMBLY_NAME__=assembly_name,
                                           __KERNEL_NAME__=kernel_name)
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
