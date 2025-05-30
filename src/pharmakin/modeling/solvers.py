from __future__ import annotations
import logging
logger = logging.getLogger(__name__)
import numpy as np
import sympy

from typing import Any, Iterable

from pharmakin.modeling import symbols


def solve_numerical(
        funcs: Iterable[sympy.core.function.AppliedUndef],
        gradients: dict[sympy.core.function.Derivative, sympy.Expr],
        ics: dict[sympy.core.function.Application, float|int],
        t_vals: np.ndarray,
        max_relative_step_size: float=0.01) -> dict[sympy.core.function.AppliedUndef, np.ndarray]:
    """Solves numerically a system of differential equations.
    This takes a simple approach and computes at each time step:
        * The current rate of change given the current values.
        * How quickly said rates are changing, given the current rates.
    The above steps are equivalent to computing a velocity given a current position, and an acceleration given
    the current velocities. For this reason, standard symbols representing position (x), velocity (v), and
    acceleration (a) are used in the code.
    
    funcs: Iterable of functions in the system. Must be function of the time parameter (symbols.t).
    gradients: dict mapping symbols for (time) derivative of functions to expressions for them.
    ics: dict representing initial conditions. Keys can be f(0) where f is a function from the funcs parameter.
    t_vals: numpy array of time values. This must begin at 0.0, as the solver assumes the initial values
        stored under .A_0 in each compund are initially valid. The reason for explicitly specifying an
        array of time values (rather than e.g. a stop-time, and step size), is to make it easier to
        compare with analytical solutions, where function values are already available as a function of time,
        and must be passed a time array to obtain a sequence of values.
    max_relative_step_size (float) - The maximum allowed change relative to current values before
        a warning is raised. For example, if this value is 0.01, and if at any time step sum of absolute changes
        in values exceed 1% of the current concentration, the logger will issue a warning.
        The logger will only warn the first time this threshold is exceeded."""

    # Check we start at t=0 (because we initialize with the initial values A_0 for each compound)
    assert t_vals[0] == 0.0
    warned = False

    # Make sure functions are ordered and assign a coordinate to each
    conc_t = list(funcs)
    x_vars = tuple(funcs)

    # Determine the rate of change in concentration of each compound (dA/dt), as a function of current concentrations
    dadt = [sympy.Derivative(a, symbols.t) for a in conc_t]
    slopes = [gradients[Ap] for Ap in dadt]
    fp = sympy.lambdify(x_vars, slopes, modules="numpy")
    
    # Express acceleration of change (d^2A/dt^2) as a function of current rate of change
    v_vars = tuple(dadt)
    second_derivatives = [sympy.Derivative(s, symbols.t).doit() for s in slopes]
    fpp = sympy.lambdify(v_vars, second_derivatives, modules="numpy")

    # Make a matrix for holding the results, and set the first row to initial values
    x = np.array([ics[A_t.subs(symbols.t, 0)] for A_t in conc_t])
    m = np.empty(shape=(len(t_vals), len(conc_t)))
    m.fill(np.nan)
    ind = 0
    m[ind, :] = x

    # Fill up the matrix using the gradients to interpolate from previous points
    for dt in np.diff(t_vals):
        ind += 1
        v = np.array(fp(*x))
        a = np.array(fpp(*v))
        step = dt*v +0.5*a*dt**2

        total_change = np.sum(np.abs(step))
        total_conc = np.sum(np.abs(x))
        toofast = total_change >= total_conc*max_relative_step_size
        if toofast and not warned:
            logger.warning(f"Step size exceeded threshold {max_relative_step_size}. Try greater temporal resolution.")
            warned = True
        x += step
        m[ind, :] = x

    # Map the label for each compound to an array of its concentrations at the input times
    d = {A_t: col for A_t, col in zip(conc_t, m.T, strict=True)}

    return d