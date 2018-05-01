# coding: utf-8

# TODO: - allow for giving a name for the trial/test basis
#       - Ni/Nj should be Ni_0/Nj_O

from gelato.expression import construct_weak_form
from gelato.calculus   import (Dot, Cross, Grad, Curl, Rot, Div)
from gelato.fem.templates import template_1d, template_header_1d
from gelato.fem.templates import template_2d, template_header_2d
from gelato.fem.templates import template_3d, template_header_3d

from numbers import Number
from collections import OrderedDict

def compile_kernel(name, expr, V,
                   namespace=globals(),
                   verbose=False,
                   d_constants={},
                   d_args={},
                   backend='python'):
    """returns a kernel from a Lambda expression on a Finite Elements space."""

    # ... parametric dimension
    dim = V.pdim
    # ...

    # ...
    if verbose:
        print('> input     := {0}'.format(expr))
    # ...

    # ...
    expr = construct_weak_form(expr, dim=dim)
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
    try:
        template = eval('template_{}d'.format(dim))
    except:
        raise ValueError('Could not find the corresponding template')

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
                template = eval('template_header_{}d'.format(dim))
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
