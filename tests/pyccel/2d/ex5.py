# coding: utf-8

# TODO: - improve how we pass arguments to weak_form

# ... import symbolic tools
weak_formulation = load('gelato.expression', 'weak_formulation', True, 2)
Grad             = load('gelato.calculus', 'Grad', False, 1)
Dot              = load('gelato.calculus', 'Dot', False, 1)
# ...

# ... weak formulation
a  = lambda x,y,u,v: Dot(Grad(u), Grad(v)) + u * v

wf        = weak_formulation(a, 2)
weak_form = lambdify(wf)
# ...

# ... to be improved
Ni   = 0.5
Ni_x = 0.2
Ni_y = 0.2
Nj   = 0.4
Nj_x = -0.2
Nj_y = -0.1

y = weak_form(0.1, 0.1, Ni, Ni_x, Ni_y, Nj, Nj_x, Nj_y)
# ...

# ... a Lambda expression can be printed
print(' a             := ', a)
print('weak-form      := ', wf)
print('eval weak-form := ', y)
print('')
# ...
