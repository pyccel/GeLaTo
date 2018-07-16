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

from .kernel import compile_kernel
from .assembly import compile_assembly

import types


# TODO add check on spaces
def discretize(name, a, spaces,
               d_constants={},
               d_args={},
               verbose=False,
               namespace=globals(),
               context=None,
               backend='python',
               export_pyfile=True):
    """."""
    # ...
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

    # ... add the assemble method to the bilinear/linear form
    if isinstance(a, BilinearForm):
        def _assemble(target):
            return assembly(target, spaces[0], spaces[1])

    elif isinstance(a, LinearForm):
        def _assemble(target):
            return assembly(target, spaces)

    a.assemble = types.MethodType(_assemble, a)
    # ...
