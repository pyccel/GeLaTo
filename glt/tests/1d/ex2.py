# coding: utf-8
#
# Usage:
#  python test_glt_1d.py
#

from vale import construct_model
import numpy as np
import matplotlib.pyplot       as plt

# ...
def run(filename):
    # ...
    from caid.cad_geometry import line
    geometry = line()

    from clapp.spl.mapping import Mapping
    mapping = Mapping(geometry=geometry)
    # ...

    # ... creates discretization parameters
    from clapp.disco.parameters.bspline import BSpline

    bspline_params = BSpline([256], [2], \
                             bc_min=[0], \
                             bc_max=[0])
    # ...

    # ... create a context from discretization
    from clapp.fema.context        import Context

    context = Context(dirname="input", \
                      discretization_params=bspline_params)
    # ...

    # ...
    pde = construct_model(filename, backend="clapp", \
                          context=context, mapping=mapping)
    # ...

    # ... accessing the pde declarations
    V      = pde["V"]
    form_a = pde["a"]
    # ...

    # ...
    assembler_a = form_a.assembler
    matrix      = form_a.matrix
    # ...

    # ... define the constants
    constants = {"r": 0.5}
    # ...

    # ...
    assembler_a.set_constants(constants, verbose=False)
    assembler_a.assemble()
    # ...

    # ... compute and plot the glt symbol
    from glt.expression import glt_symbol_from_weak_formulation
    from glt.expression import glt_plot_eigenvalues

    discretization = {"n_elements": bspline_params["N_ELEMENTS"], \
                      "degrees": bspline_params["DEGREES"]}

    expr = glt_symbol_from_weak_formulation(form_a, \
                                            verbose=True, evaluate=True, \
                                            discretization=discretization, \
                                            user_constants=constants)
    print "GLT symbol : ", expr
    # ...

    # ...
    glt_plot_eigenvalues(expr, discretization, mapping=mapping, matrix=matrix)
    # ...

    # ...
    filename_out = "1d_"+filename.split('/')[-1].split('.')[0] + ".png"
    plt.legend(loc=2)
    plt.savefig(filename_out)
    # ...

    # ...
    cmd = "rm -rf input"
    os.system(cmd)
    # ...

    # ...
    plt.clf()
    # ...

    print ("> run using ", filename, " passed.")
# ...

import clapp.common.utils      as clapp_utils

# ... initializing Clapp
clapp_utils.initialize()
# ...

import os

cmd = "rm -rf input"
os.system(cmd)

run(filename="inputs/laplace.vl")

cmd = "rm -rf input"
os.system(cmd)

# ... Finalizing Clapp
clapp_utils.finalize()
# ...
