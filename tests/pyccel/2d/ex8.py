# coding: utf-8

# TODO - Matrix object are not printed => check fcode
#        or convert it to a Tuple (not recommended)

# ... import symbolic tools
glt_function = load('pyccel.symbolic.gelato', 'glt_function', True, 3)
dx           = load('pyccel.symbolic.gelato', 'dx', False, 1)
dy           = load('pyccel.symbolic.gelato', 'dy', False, 1)
# ...

# ... weak formulation
a11 = lambda x,y,v,u: dx(u) * dx(v) + dy(u) * dy(v)
a21 = lambda x,y,v,u: 2 * dx(u) * dx(v)
a12 = lambda x,y,v,u: dy(u) * dy(v)

a   = lambda x,y,v1,v2,u1,u2: a11(x,y,u1,v1) - a12(x,y,u1,v2) + a21(x,y,u2,v1)
# ...

# ...
ga = glt_function(a, [4, 4], [2, 2])

#g = lambdify(ga)
#y = g(0.5, 0.5, 0.1, 0.3)
# ...

# ... a Lambda expression can be printed
print(' a          := ', a)
print(' glt symbol := ', ga)
print('')
#print(' symbol (0.5, 0.5, 0.1, 0.3) = ', y)
# ...
