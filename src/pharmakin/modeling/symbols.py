import sympy

# Time variable
t = sympy.Symbol("t")

# Amount
A = sympy.Function("A")
Ap = sympy.Derivative(A(t), t)
