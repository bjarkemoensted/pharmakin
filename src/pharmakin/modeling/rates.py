import numpy as np

# from sympy import Function, dsolve, Eq, Derivative, symbols
# from sympy.abc import t

# # Define the functions x(t) and y(t)
# x = Function('x')
# y = Function('y')

# # Define the system of differential equations
# eq1 = Eq(Derivative(x(t), t), x(t) + y(t))
# eq2 = Eq(Derivative(y(t), t), x(t) - y(t))

# # Solve the system with initial conditions
# sol = dsolve([eq1, eq2], ics={x(0): 1, y(0): 0})

# # Print the solutions
# for s in sol:
#     print(s)

import sympy
from sympy import Function, dsolve, Derivative, Eq, symbols


class Decay:
    def __init__(self, rate):
        pass



# Define the function y(x)
A = Function('A')

t = sympy.Symbol('t')
k = sympy.Symbol('k')

# Define the differential equation: dy/dx = y
ode = Eq(Derivative(A(t), t), k*A(t))

# Solve the ODE with initial condition y(0) = 1
solution = dsolve(ode, A(t), ics={A(0): 1})

print(solution)



def first_order_rate(t_half):
    res = np.log(2)/t_half
    return res




if __name__ == '__main__':
    pass
