from .algebra      import *
from .basic        import *
from .derivatives  import *
from .generic      import *
from .expr import *

from sympy import Integer, Float

# ... core registeries
_coeffs_registery = (Integer, Float, Constant)

_generic_ops  = (Dot, Cross,
                 Grad, Curl, Rot, Div)


# ... alias for ufl compatibility
cross = Cross
dot = Dot

Inner = Dot # TODO do we need to add the Inner class Function?
inner = Inner

grad = Grad
curl = Curl
rot = Rot
div = Div

_calculus_operators = (Grad, Dot, Inner, Cross, Rot, Curl, Div)
# ...
