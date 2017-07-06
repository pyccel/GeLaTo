# coding: utf-8
#
#
from sympy import Symbol
from sympy import latex
from sympy.simplify.simplify import simplify

# TODO find a better solution.
#      this code is duplicated in gelato.expression
ARGS_x       = ["x", "y", "z"]
ARGS_u       = ["u", "v", "w"]
ARGS_s       = ["s", "ss"]
BASIS_TEST   = "Ni"
BASIS_TRIAL  = "Nj"
BASIS_PREFIX = ["x", "y", "z", "xx", "yy", "zz", "xy", "yz", "xz"]
TOLERANCE    = 1.e-10
SETTINGS     = ["glt_integrate", \
                "glt_simplify", \
                "glt_formatting", \
                "glt_formatting_atoms"]

# ...
def latex_title_as_paragraph(title):
    """
    Returns the title as a paragraph.

    title: str
        a string for the paragraph title
    """
    return "\paragraph{" + str(title) + "}"
# ...

# ...
def glt_latex_definitions():
    """
    Returns the definitions of the atomic symbols for the GLT.
    """
    # ...
    t = Symbol('t')
    m = Symbol('m')
    s = Symbol('s')
    a = Symbol('a')
    # ...

    # ...
    def formula(symbol):
        """
        returns the latex formula for the mass symbol.
        """
        txt_m = r"\phi_{2p+1}(p+1) + 2 \sum_{k=1}^p \phi_{2p+1}(p+1-k) \cos(k \theta)"
        txt_s = r"- {\phi}''_{2p+1}(p+1) - 2 \sum_{k=1}^p {\phi}''_{2p+1}(p+1-k) \cos(k \theta)"
        txt_a = r"\phi_{2p+1}(p+1) + 2 \sum_{k=1}^p \phi_{2p+1}(p+1-k) \cos(k \theta)"

        if str(symbol) == "m":
            return txt_m
        elif str(symbol) == "s":
            return txt_s
        elif str(symbol) == "a":
            return txt_a
        else:
            print ("not yet available.")
    # ...

    # ...
    definitions = {r"m(\theta)": formula(m), \
                   r"s(\theta)": formula(s), \
                   r"a(\theta)": formula(a)}
    # ...

    return definitions
# ...

# ...
def glt_latex_names():
    """
    returns latex names for basis and atoms
    """
    # ...
    dim = 3

    symbol_names = {}
    # ...

    # ... rename basis
    B = "N"
    for i in ["i","j"]:
        Bi = B + i
        symbol_names[Symbol(Bi)] = B + "_" + i
    # ...

    # ... rename basis derivatives in the logical domain
    args_x = ARGS_x[:dim]
    args_u = ARGS_u[:dim]
    B = "N"
    for i in ["i","j"]:
        Bi = B + i
        for u in args_u + args_x:
            Bi_u = Bi + "_" + u
            partial = "\partial_" + u
            symbol_names[Symbol(Bi_u)] = partial + B + "_" + i
    # ...

    # ... rename the tensor basis derivatives
    B = "N"
    for i in ["i","j"]:
        Bi = B + i
        for k in range(0, dim):
            for s in ["", "s", "ss"]:
                Bk  = Bi + str(k+1)
                _Bk = B + "_{" + i + "_" + str(k+1) + "}"

                if len(s) > 0:
                    prime = len(s) * "\prime"

                    Bk += "_" + s
                    _Bk = B + "^{" + prime + "}" \
                            + "_{" + i + "_" + str(k+1) + "}"

                symbol_names[Symbol(Bk)] = _Bk
    # ...

    # ... TODO add flag to choose which kind of printing:
#    for k in range(0, dim):
#        # ...
#        symbol_names[Symbol('m'+str(k+1))] = "\mathfrak{m}_" + str(k+1)
#        symbol_names[Symbol('s'+str(k+1))] = "\mathfrak{s}_" + str(k+1)
#        symbol_names[Symbol('a'+str(k+1))] = "\mathfrak{a}_" + str(k+1)
#        # ...

    degree = "p"
    for k in range(0, dim):
        # ...
        for s in ["m", "s", "a"]:
            symbol_names[Symbol(s+str(k+1))] = r"\mathfrak{" + s + "}_" \
                                             + degree \
                                             + r"(\theta_" \
                                             + str(k+1) + ")"
        # ...
    # ...

    # ...
    for k in range(0, dim):
        symbol_names[Symbol("t"+str(k+1))] = r"\theta_" + str(k+1)
    # ...

    return symbol_names
# ...

# ...
def get_sympy_printer_settings(settings):
    """
    constructs the dictionary for sympy settings needed for the printer.

    settings: dict
        dictionary for different settings
    """
    sets = {}
    for key, value in settings.items():
        if key not in SETTINGS:
            sets[key] = value
    return sets
# ...

# ...
def glt_latex(expr, **settings):
    """
    returns the latex expression of expr.

    expr: sympy.Expression
        a sympy expression

    settings: dict
        dictionary for different settings
    """
    # ...
    if type(expr) == dict:
        d_expr = {}
        try:
            mode = settings["mode"]
        except:
            mode = "plain"

        sets = settings.copy()
        sets["mode"] = "plain"
        for key, txt in expr.items():
            d_expr[key] = glt_latex(txt, **sets)

        return d_expr
    # ...

    # ...
    try:
        from gelato.expression import glt_formatting
        fmt = settings["glt_formatting"]
        if fmt:
            expr = glt_formatting(expr, **settings)
    except:
        pass
    # ...

    # ...
    try:
        smp = settings["glt_simplify"]
        if smp:
            expr = simplify(expr)
    except:
        pass
    # ...

    # ...
    sets = get_sympy_printer_settings(settings)
    # ...

    return latex(expr, symbol_names=glt_latex_names(), **sets)
# ...

# ...
def print_glt_latex(expr, **settings):
    """
    Prints the latex expression of expr.

    settings: dict
        dictionary for different settings
    """
    print(glt_latex(expr, **settings))
# ...

