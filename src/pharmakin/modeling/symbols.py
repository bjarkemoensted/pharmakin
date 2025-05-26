import sympy

# Time variable
t = sympy.Symbol("t")

# Amount
A = sympy.Function("A")
A_t = A(t)
Ap = sympy.Derivative(A(t), t)
