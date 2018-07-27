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

from .symbol import compile_symbol

import types

# TODO change target to meta var, and update spaces_str
_template ="""
def {__NAME__}( target{__ARGS__}{__FIELDS__} ):
    {__DOCSTRING__}
    return {__SYMBOL_NAME__}( target, {__SPACE_ARGS__}{__ARGS__}{__FIELDS__} )
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

    # ...
    symbol = compile_symbol(symbol_name, a,
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
    code = _template.format(__NAME__=name,
                            __SYMBOL_NAME__=symbol_name,
                            __SPACE_ARGS__=spaces_str,
                            __ARGS__=args,
                            __FIELDS__=fields_str,
                            __DOCSTRING__=docstring)
    # ...

    # ...
    exec(code, namespace)
    _glt_symbol = namespace[name]
    # ...

    # ...
    setattr(a, 'discrete_spaces', spaces)
    a.symbol = types.MethodType(_glt_symbol, a)
    # ...
