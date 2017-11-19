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
    constants = {"alpha": 1., "beta": 0.5}
    # ...

    # ...
    discretization = {"n_elements": bspline_params["N_ELEMENTS"], \
                      "degrees": bspline_params["DEGREES"]}
    # ...

    # ...
    def export_plot():
        # ...
        assembler_a.set_constants(constants, verbose=False)
        assembler_a.assemble()
        # ...

        # ... compute and plot the glt symbol
        from gelato.expression import glt_symbol_from_weak_formulation
        from gelato.expression import glt_plot_eigenvalues

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
    # ...

    # ...
    def export_proof():
        # ...
        instructions = []
        # ...

        # ... compute and plot the glt symbol
        from gelato.expression import glt_symbol_from_weak_formulation

        expr = glt_symbol_from_weak_formulation(form_a, \
                                                verbose=True, evaluate=False, \
                                                discretization=discretization, \
                                                user_constants=constants, \
                                                instructions=instructions, \
                                                mode="equation", \
                                                glt_formatting=True, \
                                                glt_simplify=True, \
                                                glt_integrate="Omega")

        print "GLT symbol : ", expr
        # ...

        # ...
        from gelato.printing.latex import glt_latex_definitions
        from gelato.codegen.latex import Latex

        proof = Latex(definitions=glt_latex_definitions())
        # ...

        # ... append with definitions
        if proof.definitions is not None:
            proof.append("prelude", proof.definitions)
        # ...

        # ... append with notations
        if proof.notations is not None:
            proof.append("prelude", proof.notations)
        # ...

        # ... append instructions to the body
        proof.append("body", instructions)
        # ...

        # ...
        actual_packages = ("amsmath", "amsfonts")
        actual_packages += ("euler",)
        package_includes = "\n" + "\n".join(["\\usepackage{%s}" % p
                                             for p in actual_packages])

        preamble = r"""\documentclass[a4paper,9pt]{article}
%s


\usepackage{geometry}
 \geometry{
 a4paper,
 total={170mm,257mm},
 left=20mm,
 top=20mm,
 }

\begin{document}
""" % (package_includes)
        # ...

        # ...
        from sympy.printing.preview import preview
        latex_txt = str(proof)
        preview(latex_txt, outputTexFile='symbol.tex', preamble=preamble)
        # ...
    # ...

    # ...
#    export_plot()
    export_proof()
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
