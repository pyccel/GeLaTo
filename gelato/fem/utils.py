# coding: utf-8

# TODO: - allow for giving a name for the trial/test basis
#       - Ni/Nj should be Ni_0/Nj_O
#       - define templates as proper python functions
#       - use redbaron to modify the template

from gelato.expression import construct_weak_form
from gelato.calculus   import (Dot, Cross, Grad, Curl, Rot, Div)
from gelato.fem.templates import template_1d_scalar, template_header_1d_scalar
from gelato.fem.templates import template_2d_scalar, template_header_2d_scalar
from gelato.fem.templates import template_3d_scalar, template_header_3d_scalar

from gelato.fem.templates import template_1d_block, template_header_1d_block
from gelato.fem.templates import template_2d_block, template_header_2d_block
from gelato.fem.templates import template_3d_block, template_header_3d_block

from numbers import Number
from collections import OrderedDict

def compile_kernel(name, expr, V,
                   namespace=globals(),
                   verbose=False,
                   d_constants={},
                   d_args={},
                   backend='python'):
    """returns a kernel from a Lambda expression on a Finite Elements space."""

    from spl.fem.vector  import VectorFemSpace

    # ... parametric dimension
    dim = V.pdim
    # ...

    # ...
    if verbose:
        print('> input     := {0}'.format(expr))
    # ...

    # ...
    expr = construct_weak_form(expr, dim=dim,
                               is_block=isinstance(V, VectorFemSpace))
    if verbose:
        print('> weak form := {0}'.format(expr))
    # ...

    # ... contants
    #     for each argument, we compute its datatype (needed for Pyccel)
    #     case of Numeric Native Python types
    #     this means that a has a given value (1, 1.0 etc)
    if d_constants:
        for k, a in list(d_constants.items()):
            if not isinstance(a, Number):
                raise TypeError('Expecting a Python Numeric object')

        # update the weak formulation using the given arguments
        expr = expr.subs(d_constants)

    args = ''
    dtypes = ''
    if d_args:
        # ... additional arguments
        #     for each argument, we compute its datatype (needed for Pyccel)
        for k, a in list(d_args.items()):
            # otherwise it can be a string, that specifies its type
            if not isinstance(a, str):
                raise TypeError('Expecting a string')

            if not a in ['int', 'double', 'complex']:
                raise TypeError('Wrong type for {} :: {}'.format(k, a))

        # we convert the dictionaries to OrderedDict, to avoid wrong ordering
        d_args = OrderedDict(d_args)

        names = []
        dtypes = []
        for n,d in list(d_args.items()):
            names.append(n)
            dtypes.append(d)

        args = ', '.join('{}'.format(a) for a in names)
        dtypes = ', '.join('{}'.format(a) for a in dtypes)

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

    # ...
    if isinstance(V, VectorFemSpace):
        if V.is_block:
            n_components = len(V.spaces)

            # ... identation
            tab = ' '*4
            for i in range(0, 2*dim):
                tab += ' '*4

            tab_base = tab
            # ...

            # ... initializing accumulation variables
            lines = []
            for i in range(0, n_components):
                for j in range(0, n_components):
                    line = 'v_{i}{j} = 0.0'.format(i=i,j=j)
                    line = tab + line

                    lines.append(line)

            accum_init_str = '\n'.join(line for line in lines)
            # ...

            # .. update indentation
            for i in range(0, dim):
                tab += ' '*4
            # ...

            # ... accumulation contributions
            lines = []
            for i in range(0, n_components):
                for j in range(0, n_components):
                    line = 'v_{i}{j} += ({__WEAK_FORM__}) * wvol'
                    line = line.format(i=i, j=j, __WEAK_FORM__=expr[i,j])
                    line = tab + line

                    lines.append(line)

            accum_str = '\n'.join(line for line in lines)
            # ...

            # ... assign accumulated values to element matrix
            e_pattern = 'mat[{i}, {j}, il_1, il_2, p1 + jl_1 - il_1, p2 + jl_2 - il_2] = v_{i}{j}'
            tab = tab_base
            lines = []
            for i in range(0, n_components):
                for j in range(0, n_components):
                    line = e_pattern.format(i=i,j=j)
                    line = tab + line

                    lines.append(line)

            accum_assign_str = '\n'.join(line for line in lines)
            # ...

            code = template.format(__KERNEL_NAME__=name,
                                   __ACCUM_INIT__=accum_init_str,
                                   __ACCUM__=accum_str,
                                   __ACCUM_ASSIGN__=accum_assign_str,
                                   __ARGS__=args)

            print(code)
            import sys; sys.exit(0)

        else:
            raise NotImplementedError('We only treat the case of a block space, for '
                                      'which all components have are identical.')

    else:
        code = template.format(__KERNEL_NAME__=name,
                               __WEAK_FORM__=expr,
                               __ARGS__=args)

    # ...

    # ... always execute it to check if it is python valid
    exec(code, namespace)
    kernel = namespace[name]
    # ...

    # ...
    if backend == 'fortran':
        try:
            # import epyccel function
            from pyccel.epyccel import epyccel

            #  ... define a header to specify the arguments types for kernel
            try:
                template = eval('template_header_{dim}d_{pattern}'.format(dim=dim,
                                                                          pattern=pattern))
            except:
                raise ValueError('Could not find the corresponding template')

            header = template.format(__KERNEL_NAME__=name,
                                     __TYPES__=dtypes)
            # ...

            # compile the kernel
            kernel = epyccel(code, header, name=name)
        except:
            print('> COULD NOT CONVERT KERNEL TO FORTRAN')
            print('  THE PYTHON BACKEND WILL BE USED')
    # ...

    return kernel
