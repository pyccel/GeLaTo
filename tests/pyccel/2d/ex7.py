# coding: utf-8

# TODO: - improve how we pass arguments to weak_form

# ... import symbolic tools
weak_formulation = load('pyccel.symbolic.gelato', 'weak_formulation', True, 2)
dx               = load('pyccel.symbolic.gelato', 'dx', False, 1)
dy               = load('pyccel.symbolic.gelato', 'dy', False, 1)
# ...

# ... weak formulation
eps = 0.1
bracket = lambda u,v: dy(u)*dx(v) - dx(u)*dy(v)
a       = lambda x,y,u,v: (1.0 + eps*x)*dx(u) * dx(v) + bracket(u,v)

wf        = weak_formulation(a, 2)
weak_form = lambdify(wf)
# ...

# ... to be improved
Ni_x = 0.2
Ni_y = 0.2
Nj_x = -0.2
Nj_y = -0.1

y = weak_form(0.1, 0.1, Ni_x, Ni_y, Nj_x, Nj_y)
# ...

# ... a Lambda expression can be printed
print((' a             := ', a))
print(('weak-form      := ', wf))
print(('eval weak-form := ', y))
print('')
# ...
