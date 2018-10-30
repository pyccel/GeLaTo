# coding: utf-8

from collections import OrderedDict
from collections import namedtuple

from pyccel.ast import Nil

from sympde.core import BilinearForm

from spl.api.codegen.utils import write_code

from gelato.codegen.ast import Kernel
from gelato.codegen.ast import Interface
from gelato.printing.pycode import pycode

import os
import importlib
import string
import random
import math

SPL_DEFAULT_FOLDER = '__pycache__/gelato'

def random_string( n ):
    chars    = string.ascii_uppercase + string.ascii_lowercase + string.digits
    selector = random.SystemRandom()
    return ''.join( selector.choice( chars ) for _ in range( n ) )


class DiscreteSymbol(object):

    def __init__(self, a, Vh, namespace=globals(), to_compile=True, module_name=None):

        if not isinstance(a, BilinearForm):
            raise TypeError('> Expecting a symbolic BilinearForm')

        # ...
        kernel = Kernel(a, Vh)
        interface = Interface(kernel)
        # ...

        # ...
        self._expr = a
        self._discrete_space = Vh
        self._tag = kernel.tag
        self._mapping = None
        self._interface = interface
        self._dependencies = self.interface.dependencies
        # ...

        # generate python code as strings for dependencies
        self._dependencies_code = self._generate_code()

        self._dependencies_fname = None
        self._interface_code = None
        self._func = None
        if to_compile:
            # save dependencies code
            self._save_code(module_name=module_name)

            # generate code for Python interface
            self._generate_interface_code()

            # compile code
            self._compile(namespace)

    @property
    def expr(self):
        return self._expr

    @property
    def discrete_space(self):
        return self._discrete_space

    @property
    def tag(self):
        return self._tag

    @property
    def mapping(self):
        return self._mapping

    @property
    def interface(self):
        return self._interface

    @property
    def dependencies(self):
        return self._dependencies

    @property
    def interface_code(self):
        return self._interface_code

    @property
    def dependencies_code(self):
        return self._dependencies_code

    @property
    def dependencies_fname(self):
        return self._dependencies_fname

    @property
    def dependencies_modname(self):
        module_name = os.path.splitext(self.dependencies_fname)[0]
        module_name = module_name.replace('/', '.')
        return module_name

    @property
    def func(self):
        return self._func

    def _generate_code(self):
        # ... generate code that can be pyccelized
        code = ''
        for dep in self.dependencies:
            code = '{code}\n{dep}'.format(code=code, dep=pycode(dep))
        # ...
        return code

    def _save_code(self, module_name=None):
        folder = SPL_DEFAULT_FOLDER

        code = self.dependencies_code
        if module_name is None:
            module_name = 'dependencies_{}'.format(self.tag)
        self._dependencies_fname = write_code(module_name, code, ext='py', folder=folder)

    def _generate_interface_code(self, module_name=None):
        imports = []

        # ... generate imports from dependencies module
        pattern = 'from {module} import {dep}'

        if module_name is None:
            module_name = self.dependencies_modname

        for dep in self.dependencies:
            txt = pattern.format(module=module_name, dep=dep.name)
            imports.append(txt)
        # ...

        # ...
        imports = '\n'.join(imports)
        # ...

        code = pycode(self.interface)

        self._interface_code = '{imports}\n{code}'.format(imports=imports, code=code)

    def _compile(self, namespace, module_name=None):

        if module_name is None:
            module_name = self.dependencies_modname

        # ...
        dependencies_module = importlib.import_module(module_name)
        # ...

        # ...
        code = self.interface_code
        name = self.interface.name

        exec(code, namespace)
        interface = namespace[name]
        # ...

        self._func = interface

    def _check_arguments(self, **kwargs):

        # TODO do we need a method from Interface to map the dictionary of arguments
        # that are passed for the call (in the same spirit of build_arguments)
        # the idea is to be sure of their order, since they can be passed to
        # Fortran

        _kwargs = {}

        # ... mandatory arguments
        sym_args = self.interface.in_arguments
        keys = [str(a) for a in sym_args]
        for key in keys:
            try:
                _kwargs[key] = kwargs[key]
            except:
                raise KeyError('Unconsistent argument with interface')
        # ...

        # ... optional (inout) arguments
        sym_args = self.interface.inout_arguments
        keys = [str(a) for a in sym_args]
        for key in keys:
            try:
                _kwargs[key] = kwargs[key]
            except:
                pass
        # ...

        return _kwargs

    @property
    def spaces(self):
        return self._spaces

    def evaluate(self, *args, **kwargs):
#        newargs = tuple(self.discrete_spaces)
#
#        if self.mapping:
#            newargs = newargs + (self.mapping,)

        kwargs = self._check_arguments(**kwargs)

#        return self.func(*newargs, **kwargs)
        return self.func(*args, **kwargs)
