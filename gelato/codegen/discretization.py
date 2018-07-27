# -*- coding: utf-8 -*-



from sympy.core.containers import Tuple
from sympy import Matrix
from sympy import Integer, Float

from numbers import Number
from collections import OrderedDict

from numpy import unique
import os
import importlib

from symfe.core import gelatize
from symfe.core import BilinearForm, LinearForm, FunctionForm
from symfe.core import Constant
from symfe.core import Field
from symfe.codegen.utils import arguments_datatypes_as_dict
from symfe.codegen.utils import arguments_datatypes_split

from spl.fem.splines import SplineSpace
from spl.fem.tensor  import TensorFemSpace

from .symbol import compile_symbol

import types

# ... TODO move this in another file
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
# ...

# TODO change target to meta var, and update spaces_str
_template ="""
def {__NAME__}( target, {__X_ARGS__}{__T_ARGS__}{__ARGS__}{__FIELDS__}{__MAT_ARGS__} ):
    {__DOCSTRING__}
    return {__SYMBOL_NAME__}( {__X_ARGS__}{__T_ARGS__}{__ARGS__}{__FIELDS__}{__MAT_ARGS__} )
"""

_template_docstring = """
\"\"\"
Evaluates the GLT symbol.

This method is calling a function that has been automatically generated:
    - {__SYMBOL_NAME__} low-level symbol evaluation function
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
def discretize_symbol(a, spaces,
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
    dim = a.ldim
    # ...

    # ...
    if root is None:
        root = '.pyccel'
    # ...

    # ...
    if name is None:
        name = 'glt_symbol'

    symbol_name = '{name}_{hash}'.format(name=name,
                                         hash=abs(hash(a)))
    # ...

    # ... TODO improve + checks
    V = spaces[0]
    if isinstance(V, SplineSpace):
        degrees = [V.degree]
        n_elements = [V.ncells]

    elif isinstance(V, TensorFemSpace):
        degrees = V.degree
        n_elements = V.ncells

    else:
        raise NotImplementedError('')
    # ...

    # ...
    symbol = compile_symbol(symbol_name, a, degrees,
                            n_elements=n_elements,
                            verbose=verbose,
                            namespace=namespace,
                            context=context,
                            backend=backend,
                            export_pyfile=export_pyfile)
    # ...

    # ...
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
    args_docstring = docstring_arguments(a.constants, d_args)
    # ...

    # ...
    docstring = _template_docstring.format(__SYMBOL_NAME__=symbol_name,
                                           __PARAMETERS__=args_docstring)
    # identation (def function body)
    tab = ' '*4

    lines = []
    for line in docstring.split('\n'):
        lines.append(tab + line)
    docstring = '\n'.join(i for i in lines)
    # ...

    # ...
    x_args_str = print_position_args(dim)
    t_args_str = print_fourier_args(dim)
    # ...

    # ...
    mat_args_str = print_mat_args()
    # ...

    # ...
    code = _template.format(__NAME__=name,
                            __SYMBOL_NAME__=symbol_name,
                            __X_ARGS__=x_args_str,
                            __T_ARGS__=t_args_str,
                            __MAT_ARGS__=mat_args_str,
                            __ARGS__=args,
                            __FIELDS__=fields_str,
                            __DOCSTRING__=docstring)
    # ...

#    print('--------------')
#    print(code)
#    print('--------------')
#    import sys; sys.exit(0)

    # ...
    exec(code, namespace)
    _glt_symbol = namespace[name]
    # ...

    # ...
    setattr(a, 'discrete_spaces', spaces)
    a.symbol = types.MethodType(_glt_symbol, a)
    # ...
