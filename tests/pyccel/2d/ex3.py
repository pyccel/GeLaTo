# coding: utf-8

# ... import symbolic tools
glt_function = load('gelato.glt', 'glt_function', True, 3)
dx           = load('gelato.calculus', 'dx', False, 1)
dy           = load('gelato.calculus', 'dy', False, 1)
# ...

# ... weak formulation
#     here the weak formulation is based on a user-defined operator 'bracket'
#     we can omit the coordinates 'x,y' when defining 'bracket' as long as we
#     call it in 'b' in a consistent way.
bracket = lambda u,v: dy(u)*dx(v) - dx(u)*dy(v)
b       = lambda x,y,u,v: dx(u) * dx(v) + bracket(u,v)
# ...

# ... computing the glt symbol
gb = glt_function(b, [4, 4], [2, 2])
# ...

# ...
print(' b          := ', b)
print(' glt symbol := ', gb)
# ...
