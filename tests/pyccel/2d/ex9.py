# coding: utf-8

# TODO - improve lambdifying wf: use Matrix instead of Tuple?

# ... import symbolic tools
weak_formulation = load('pyccel.symbolic.gelato', 'weak_formulation', True, 2)
dx               = load('pyccel.symbolic.gelato', 'dx', False, 1)
dy               = load('pyccel.symbolic.gelato', 'dy', False, 1)
# ...

# ... weak formulation
eps = 0.1

a11 = lambda x,y,v,u: dx(u) * dx(v) + dy(u) * dy(v)
a21 = lambda x,y,v,u: 2 * dx(u) * dx(v)
a12 = lambda x,y,v,u: dy(u) * dy(v)
a22 = lambda x,y,v,u: eps * u * v

a = lambda x,y,v1,v2,u1,u2: a11(x,y,u1,v1) + a12(x,y,u1,v2) + a21(x,y,u2,v1) + a22(x,y,u2,v2)
# ...

# ...
wf = weak_formulation(a, 2)

#weak_form = lambdify(wf)
# ...

# ... a Lambda expression can be printed
print(' a             := ', a)
print('weak-form      := ', wf)
#print('eval weak-form := ', y)
print('')
# ...
