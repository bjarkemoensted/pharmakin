from __future__ import annotations
import networkx as nx
import numpy as np
import sympy


# A	Instantaneous Input	Drug appears immediately in systemic circulation.
# B	Zero-Order Absorption	Drug enters bloodstream at a constant rate over time.
# C	First-Order Absorption	Rate of absorption is proportional to drug remaining at site.
# D	Delayed Absorption	Time lag before absorption begins. Then behaves like first-order or zero-order.
# E	Controlled Release / Sustained	Complex release: could be pseudo-zero-order or multiphasic.



# Symbol for time variable
_t = sympy.Symbol("t")


def _make_symbol(name: str):
    res = sympy.Symbol(name, real=True, positive=True)
    return res


class Quantity:
    def __init__(self, name, model: Model):
        self.name = name
        self._model = model
        self.sym = sympy.Function(name)(_t)
    #


class Relation:
    def __init__(self, *eqs: sympy.Eq):
        self.equations = eqs




class Model:
    def __init__(self):
        self.G = nx.DiGraph()
        self._nodes = dict()
    
    def add_quantity(self, name: str):
        if name in self._nodes:
            return
        self.G.add_node(name)
        q = Quantity(name=name, model=self)
        self._nodes[name] = q
    
    def add_quantities(self, *names: str):
        for name in names:
            self.add_quantity(name=name)
    
    def add_transition(self, u: str, v: str, t_half: float):
        self.add_quantities(u, v)
        # TODO expand to allow non-linear metabolism, but play around with first-order first !!!
        self.G.add_edge(u, v)
        
        k = float(np.log(2))/t_half
        gamma = Rate(u=self[u], v=self[v], t_half=t_half)
        print(gamma.out, type(gamma.eq), isinstance(gamma.eq, sympy.Equality))
        
        a = sympy.Symbol("A")
        b = sympy.Symbol("B")
        
        print(gamma.eq)
        e1 = sympy.Eq()
        wat = Relation()
    
    def __getitem__(self, name: str):
        return self._nodes[name]
    
if __name__ == '__main__':
    pass
    # m = Model()
    # m.add_quantity("LDX")
    # m.add_transition("LDX", "AMP", t_half=1.5)
    
    # print(m["AMP"])

