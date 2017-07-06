# coding: utf-8
#
#

from sympy import Symbol

def assign(a, b):
    """
    prints the instruction that assigns b to a.

    a: str
        a string expression

    b: str
        a string expression
    """
    txt = str(a) + " = " + str(b)
    return txt

class Latex(object):
    """
    A Latex generator class.
    """
    def __init__(self, notations=None, definitions=None, lemmas=None, rules=None):
        """
        Initialization of the latex.
        """
        self._definitions = definitions
        self._notations   = notations
        self._lemmas      = lemmas
        self._rules       = rules
        self._txt         = ""
        self._body        = []
        self._prelude     = []

    @property
    def definitions(self):
        """
        Returns the definitions section
        """
        return self._definitions

    @property
    def notations(self):
        """
        Returns the notations section
        """
        return self._notations

    @property
    def lemmas(self):
        """
        Returns the lemmas section
        """
        return self._lemmas

    @property
    def rules(self):
        """
        Returns the rules section
        """
        return self._rules

    @property
    def prelude(self):
        """
        Returns the prelude
        """
        return self._prelude

    @property
    def body(self):
        """
        Returns the body
        """
        return self._body

    def append(self, section, vertex, **settings):
        """
        decorates a given vertex using the settings.

        section: str
            one of {"body", ""notation"}

        vertex: dict
            one of {"definitions", "notations", "lemmas", "rules"}

        settings: dict
            some input attributs.
        """
        if section == "prelude":
            for key, txt in vertex.items():
                self._prelude.append(assign(key, txt))
        elif section == "body":
            self._body += vertex
        else:
            print ("> Wrong section name.")
            raise()

    def _str_prelude(self, mode='plain'):
        """
        Prints the prelude

        mode: str
            printing mode as used in sympy latex
        """
        # ...
        txt  = ""
        txt += "Definitions:" + "\n"
        # ...

        # ...
        if mode is "plain":
            txt += r"\begin{align}" + "\n"
        # ...

        # ...
        for line in self.prelude[:-1]:
            txt += line + "\n"
            txt += r"\\"
            txt += "\n"
        txt += self.prelude[-1] + "\n"
        # ...

        # ...
        if mode is "plain":
            txt += r"\end{align}" + "\n"
        # ...

        return txt

    def _str_body(self):
        """
        Prints the body
        """
        txt = ""

#        if mode is "plain":
#            txt += r"\begin{align}" + "\n"

        for line in self.body:
            txt += line + "\n"

#        if mode is "plain":
#            txt += r"\end{align}" + "\n"

        return txt

    def __str__(self):
        """
        Returns the whole document as a string
        """
        # ...
        txt = ""
        # ...

        # ...
        try:
            txt += self._str_prelude()
        except:
            pass
        # ...

        # ...
        try:
            txt += self._str_body()
        except:
            pass
        # ...

        return txt
