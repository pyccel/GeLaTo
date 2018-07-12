# coding: utf-8

# TODO: - allow for giving a name for the trial/test basis
#       - Ni/Nj should be Ni_0/Nj_O
#       - define templates as proper python functions
#       - use redbaron to modify the template

from sympy.core.containers import Tuple

from gelato.expression import _is_base_function
from gelato.expression import BASIS_PREFIX

from gelato.fem.expr import gelatize
from gelato.fem.expr import BilinearForm, LinearForm
from gelato.calculus   import Constant
from gelato.calculus   import Field
from gelato.fem.templates import template_1d_scalar, template_header_1d_scalar
from gelato.fem.templates import template_2d_scalar, template_header_2d_scalar
from gelato.fem.templates import template_3d_scalar, template_header_3d_scalar

from gelato.fem.templates import template_1d_block, template_header_1d_block
from gelato.fem.templates import template_2d_block, template_header_2d_block
from gelato.fem.templates import template_3d_block, template_header_3d_block

from gelato.fem.templates import symbol_1d_scalar, symbol_header_1d_scalar
from gelato.fem.templates import symbol_2d_scalar, symbol_header_2d_scalar
from gelato.fem.templates import symbol_3d_scalar, symbol_header_3d_scalar

from gelato.fem.templates import symbol_1d_block, symbol_header_1d_block
from gelato.fem.templates import symbol_2d_block, symbol_header_2d_block
from gelato.fem.templates import symbol_3d_block, symbol_header_3d_block

from gelato.fem.templates import eval_field_1d_scalar
from gelato.fem.templates import eval_field_2d_scalar
from gelato.fem.templates import eval_field_3d_scalar


from numbers import Number
from collections import OrderedDict
from sympy import Integer, Float
import os

def _convert_int_to_float(expr):
    sub = zip( expr.atoms(Integer), map(Float, expr.atoms(Integer)) )
    expr = expr.subs(sub)
    return expr

def _count_letter(word, char):
    count = 0
    for c in word:
        if c == char:
            count += 1
    return count

def construct_test_functions(nderiv, dim):
    """constructs test functions and their derivatives for every direction k.
    on return, we get a list of statements, that we need to indent later
    """
    d_basis = OrderedDict()
    for k in range(1, dim+1):
        d_basis[k,0] = 'test_bs{k}[il_{k}, 0, g{k}]'.format(k=k)
        for d in range(1, nderiv+1):
            d_basis[k,d] = 'test_bs{k}[il_{k}, {d}, g{k}]'.format(d=d, k=k)

    return d_basis

def construct_trial_functions(nderiv, dim):
    """constructs trial functions and their derivatives for every direction k.
    on return, we get a list of statements, that we need to indent later
    """
    d_basis = OrderedDict()
    for k in range(1, dim+1):
        d_basis[k,0] = 'trial_bs{k}[jl_{k}, 0, g{k}]'.format(k=k)
        for d in range(1, nderiv+1):
            d_basis[k,d] = 'trial_bs{k}[jl_{k}, {d}, g{k}]'.format(d=d, k=k)

    return d_basis


def mkdir_p(dir):
    if os.path.isdir(dir):
        return
    os.makedirs(dir)

def write_code(name, code, ext='py', folder='.pyccel'):
    filename = '{name}.{ext}'.format(name=name, ext=ext)
    if folder:
        mkdir_p(folder)
        filename = os.path.join(folder, filename)

    f = open(filename, 'w')
    for line in code:
        f.write(line)
    f.close()


def compile_kernel(name, a, spaces,
                   d_constants={},
                   d_args={},
                   verbose=False,
                   namespace=globals(),
                   context=None,
                   backend='python',
                   export_pyfile=True):
    """."""
    from spl.fem.vector  import VectorFemSpace
    from spl.fem.splines import SplineSpace
    from spl.fem.tensor  import TensorFemSpace

    if not isinstance(a, (BilinearForm, LinearForm)):
           raise TypeError('Expecting BilinearForm, LinearForm')

    if not isinstance(spaces, (tuple, list, Tuple)):
        raise TypeError('Expecting tuple, list, Tuple')

    test_space = spaces[0]
    trial_space = spaces[1]

    if not(test_space is trial_space):
        raise NotImplementedError('TODO')

    # TODO: nderiv must be computed from the weak form
    nderiv = 1

    # dimension
    dim = a.ldim

    # ...
    fields = a.fields
    if verbose:
        print('> Fields = ', fields)
    # ...

    # ... TODO improve
    V = a.trial_spaces[0]
    U = a.test_spaces[0]
    TEST_BASIS = 'Ni'
    TRIAL_BASIS = 'Nj'
    expr = gelatize(a, basis={V: TRIAL_BASIS, U: TEST_BASIS})

    if verbose:
        print('> gelatized  >>> {0}'.format(expr))

    is_test_function = lambda a: _is_base_function(a, TEST_BASIS)
    is_trial_function = lambda a: _is_base_function(a, TRIAL_BASIS)
    # ...

    # ... contants
    #     for each argument, we compute its datatype (needed for Pyccel)
    #     case of Numeric Native Python types
    #     this means that a has a given value (1, 1.0 etc)
    if d_constants:
        for k, arg in list(d_constants.items()):
            if not isinstance(arg, Number):
                raise TypeError('Expecting a Python Numeric object')

        # update the weak formulation using the given arguments
        _d = {}
        for k,v in list(d_constants.items()):
            if isinstance(k, str):
                _d[Constant(k)] = v
            else:
                _d[k] = v

        expr = expr.subs(_d)

    args = ''
    dtypes = ''
    if d_args:
        # ... additional arguments
        #     for each argument, we compute its datatype (needed for Pyccel)
        for k, arg in list(d_args.items()):
            # otherwise it can be a string, that specifies its type
            if not isinstance(arg, str):
                raise TypeError('Expecting a string')

            if not arg in ['int', 'double', 'complex']:
                raise TypeError('Wrong type for {} :: {}'.format(k, arg))

        # we convert the dictionaries to OrderedDict, to avoid wrong ordering
        d_args = OrderedDict(sorted(list(d_args.items())))

        names = []
        dtypes = []
        for n,d in list(d_args.items()):
            names.append(n)
            dtypes.append(d)

        args = ', '.join('{}'.format(arg) for arg in names)
        dtypes = ', '.join('{}'.format(arg) for arg in dtypes)

        args = ', {}'.format(args)
        dtypes = ', {}'.format(dtypes)

        # TODO check what are the free_symbols of expr,
        #      to make sure the final code will compile
        #      the remaining free symbols must be the trial/test basis functions,
        #      and the coordinates
    # ...

    # ...
    if isinstance(V, VectorFemSpace) and not( V.is_block ):
        raise NotImplementedError('We only treat the case of a block space, for '
                                  'which all components have are identical.')
    # ...

    # ...
    pattern = 'scalar'
    if isinstance(V, VectorFemSpace):
        if V.is_block:
            pattern = 'block'

        else:
            raise NotImplementedError('We only treat the case of a block space, for '
                                      'which all components have are identical.')
    # ...

    # ...
    template_str = 'template_{dim}d_{pattern}'.format(dim=dim, pattern=pattern)
    try:
        template = eval(template_str)
    except:
        raise ValueError('Could not find the corresponding template {}'.format(template_str))
    # ...

    # ... identation (def function body)
    tab = ' '*4
    # ...

    # ... field coeffs
    if fields:
        raise NotImplementedError('')

    else:
        field_coeffs_str = ''
        eval_field_str   = ''
        field_value_str  = ''
        field_types_str  = ''
    # ...

    # ... compute indentation
    tab_base = tab
    for i in range(0, 3*dim):
        tab += ' '*4
    # ...

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
    if isinstance(V, VectorFemSpace):
        raise NotImplementedError('TODO')
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

    print('--------------')
    print(code)
    print('--------------')

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




