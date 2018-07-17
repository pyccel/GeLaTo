# coding: utf-8

# TODO: - allow for giving a name for the trial/test basis
#       - Ni/Nj should be Ni_0/Nj_0
#       - define templates as proper python functions
#       - use redbaron to modify the template

#     NOTE: THE PATH OF TEMPLATES IS HARD CODED!


from sympy import Integer, Float

from collections import OrderedDict

from numpy import unique
import os
import importlib

# ...
def _is_base_function(a, base):
    """Returns True if the atom is a test/trial function or a partial derivative of a
    test/trial function.

    base: str
        must be `Ni` or `Nj`
    """
    name = a.name
    if name == base:
        return True

    if not name.startswith('{}_'.format(base)):
        return False

    # name starts with `Ni_` or `Nj_`
    names = name[3:]

    # letters
    ls = [i for i in names]
    ls = unique(ls)

    # we should have only x,y or z
    for i in ls:
        if not(i in ['x', 'y', 'z']):
            return False

    return True
# ...

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

def arguments_datatypes_as_dict(arguments):
    d_args = {}
    for c in arguments:
        if c.is_integer:
            dtype = 'int'

        elif c.is_real:
            dtype = 'double'

        elif c.is_complex:
            dtype = 'complex'

        else:
            dtype = 'double'
            print('> found a constant without prescribed type. Default type (double) will be used!')

        d_args[c.name] = dtype
    # ...

    return d_args

def arguments_datatypes_split(d_args):

    # ... contants
    #     for each argument, we compute its datatype (needed for Pyccel)
    #     case of Numeric Native Python types
    #     this means that a has a given value (1, 1.0 etc)
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
    # ...

    return args, dtypes
