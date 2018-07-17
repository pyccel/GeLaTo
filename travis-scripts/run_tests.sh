# this file can be used by developers to run all tests localy before travis
# it must be called from the root directory

# TODO not working
#python3 gelato/core/tests/test_calculus.py

python3 gelato/core/tests/test_derivatives.py
python3 gelato/core/tests/test_space.py
python3 gelato/core/tests/test_expr_1d.py
python3 gelato/core/tests/test_expr_2d.py
python3 gelato/core/tests/test_expr_3d.py

python3 gelato/fem/tests/test_kernel_1d.py
python3 gelato/fem/tests/test_kernel_2d.py
python3 gelato/fem/tests/test_kernel_3d.py
python3 gelato/fem/tests/test_assembly_1d.py
python3 gelato/fem/tests/test_assembly_2d.py
python3 gelato/fem/tests/test_assembly_3d.py
python3 gelato/fem/tests/test_pde_1d.py
python3 gelato/fem/tests/test_pde_2d.py
python3 gelato/fem/tests/test_pde_3d.py
