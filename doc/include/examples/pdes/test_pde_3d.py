# coding: utf-8
#
# Usage:
#  python test_pde_3d.py
#

from clapp.vale.vale           import Vale
from clapp.fema.context        import Context
from clapp.fema.assembler      import Assembler
from clapp.plaf.parameters.ddm import Ddm as DDM_parameters
from clapp.spl.mapping         import Mapping
from clapp.plaf.linear_solver  import Linear_solver
from clapp.plaf.vector         import Vector
import clapp.common.utils      as clapp_utils
import matplotlib.pyplot       as plt

# ...
def run(filename):
    # ...
    from caid.cad_geometry import cube
    geometry = cube()
    mapping = Mapping(geometry=geometry)
    # ...

    # ...
    import os

    cmd = "mkdir -p output"
    os.system(cmd)
    # ...

    # ... creates discretization parameters
    from clapp.disco.parameters.bspline import BSpline

    bspline_params = BSpline([8,8,8], [2,2,2], \
                             bc_min=[0,0,0], \
                             bc_max=[0,0,0])
    # ...

    # ... create a context from discretization
    context = Context(dirname="input", \
                      discretization_params=bspline_params)
    # ...

    # ...
    vale = Vale(filename=filename)
    pde  = vale.parse(context, mapping)
    # ...

    # ... accessing the pde declarations
    V           = pde["V"]
    u           = pde["u"]
    form_a      = pde["a"]
    form_b      = pde["b"]

    assembler_a = form_a.assembler
    matrix      = form_a.matrix
    assembler_b = form_b.assembler
    rhs         = form_b.vector

    assembler_a.assemble()
    assembler_b.assemble()
    # ...

    # ...
    from clapp.plaf.parameters.linear_solver import LAPACK_LU
    from clapp.plaf.parameters.linear_solver import DRIVER

    params = DRIVER(solver=LAPACK_LU())
    linsol = Linear_solver(matrix=matrix, dirname="input", parameters=params)
    # ...

    # ...
    y = linsol.solve(rhs)
    # ...

    # ... exports the field
    u.set(y)
    # ...

    # ... exports phi to vtk file. Can be used in Visit and Paraview
    filename_out = "uh_3d_"+filename.split('/')[-1].split('.')[0] + ".vtk"
    u.to_vtk(filename_out, mapping=mapping, n_pts=20)
    # ...

    # ... compute L2 error
    x = u.compute_l2_error(mapping=mapping, \
                           function_name="u_analytic")
    print (">> norms : ", x)
    # ...

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

run("pdes/3d/poisson.vl"  )

cmd = "rm -rf input"
os.system(cmd)

# ... Finalizing Clapp
clapp_utils.finalize()
# ...
