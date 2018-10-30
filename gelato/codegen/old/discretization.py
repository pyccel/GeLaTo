# -*- coding: utf-8 -*-

# TODO: - improve docstrings with M arguments
#       - compute n_rows and n_cols from BilinearForm
#       - add check on spaces
#       - pass root to compile kernel and assembly



from sympy.core.containers import Tuple
from sympy import Matrix
from sympy import Integer, Float

from numbers import Number
from collections import OrderedDict

from numpy import unique
import os
import importlib

from sympde.core import gelatize
from sympde.core import BilinearForm, LinearForm, FunctionForm
from sympde.core import Constant
from sympde.core import Field
from sympde.codegen.utils import arguments_datatypes_as_dict
from sympde.codegen.utils import arguments_datatypes_split

from spl.fem.splines import SplineSpace
from spl.fem.tensor  import TensorFemSpace

from .symbol import compile_symbol
from .utils import (print_position_args, print_fourier_args,
                    construct_x_args_names,
                    construct_t_args_names,
                    print_mat_args,
                    construct_argument_matrix_name,
                    print_argument_matrix_kwargs,
                    construct_matrix_names,
                    print_matrix_args,
                    print_define_matrix,
                    docstring_arguments)

import types

_template ="""
def {__NAME__}( target, {__X_ARGS__}{__T_ARGS__}{__ARGS__}{__FIELDS__}{__MAT_KWARGS__} ):
    {__DOCSTRING__}
    {__MAT_DEC__}
    {__SYMBOL_NAME__}( {__X_ARGS__}{__T_ARGS__}, {__MAT_ARGS__}{__ARGS__}{__FIELDS__} )
    return {__MAT_ARGS__}
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

    n_rows = 1 ; n_cols = 1
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
    x_args = construct_x_args_names(dim)
    t_args = construct_t_args_names(dim)

    x_args_str = print_position_args(x_args)
    t_args_str = print_fourier_args(t_args)
    # ...

    # ...
    mat_args_str = print_mat_args()
    # ...

    # ...
    argument_mat = construct_argument_matrix_name(n_rows, n_cols)
    mat_args = construct_matrix_names(n_rows, n_cols)

    mat_kwargs = print_argument_matrix_kwargs(argument_mat)
    mat_args_str = print_matrix_args(n_rows, n_cols, mat_args)
    mat_decs_str = print_define_matrix(n_rows, n_cols, x_args, mat_args, argument_mat, tab)
    # ...

    # ...
    code = _template.format(__NAME__=name,
                            __SYMBOL_NAME__=symbol_name,
                            __X_ARGS__=x_args_str,
                            __T_ARGS__=t_args_str,
                            __MAT_KWARGS__=mat_kwargs,
                            __MAT_DEC__=mat_decs_str,
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
