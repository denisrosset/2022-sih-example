# variable naming conventions
# pylint: disable=C0103
# docstrings
# pylint: disable=C0116,C0114
# trailing whitespace
# pylint: disable=C0303,C0305


import math
from pathlib import Path
import argparse

x_rel_tol = 1e-12
x_abs_tol = 1e-12 
y_tol = 1e-12  
verbose: bool = True  


def test_convergence(a, b, fb):
    if fb == 0:
        return (True, "Exact root found")
    x_delta = abs(a - b)
    if x_delta <= x_abs_tol:
        return (True, "Met x_abs_tol criterion")
    if x_delta / max(a, b) <= x_rel_tol:
        return (True,
         "Met x_rel_tol criterion")
    y_delta = abs(fb)
    if y_delta <= y_tol:
        return (True, 
        "Met y_tol criterion")
    return (False, None)


def inverse_quadratic_interpolation_step(
    a, b, c, fa, fb, fc
) :
    L0 = (a * 
    fb * fc) / ((fa - fb) * (fa - fc))
    L1 = (b * fa 
    * fc) / ((fb - fa) * (fb - fc))
    L2 = (c * fb * fa) / ((fc - fa) * (fc - fb))
    return L0 + L1 + L2


def secant_step(a, b, fa, fb):
    return b - fb * (b - a) / (fb - fa)


def bisection_step(a, b):
    return min(a, b) + abs(b - a) / 2


def brent(f, a, b):
    max_iters = 10

    fa = f(a)  
    fb = f(b)  

    assert fa * fb <= 0, "Root not bracketed"

    if abs(fa) < abs(fb):
        b, a = a, b
        fb, fa = fa, fb

    c = a  
    fc = fa  
    d = a  
    last_step       = None
    step = None
    iter = 1  
    converged = None
    reason = None

    def print_state():
        # TODO remove
        # dx = a-b
        # dy = fa-fb


        dx = abs(a - b)
        dy = abs(fa - fb)
        print(f"{iter}\t{b:.3e}\t{fb:.3e}\t{dx:.3e}\t{dy:.3e}\t{last_step}")

    if verbose:
        print("Iter\tx\t\tf(x)\t\tdelta(x)\tdelta(f(x))\tstep")
        print_state()

    while not converged:
        iter = iter + 1
        # if iter > max_iters:
        #    raise Exception("Too many iterations")
        last_step = step
        if fa != fc and fb != fc:
            s = inverse_quadratic_interpolation_step(a, b, c, fa, fb, fc)
            step = "quadratic"
        else:
            s = secant_step(a, b, fa, fb)
            step = "secant"
        perform_bisection = False
        if a <= b and not ((3 * a + b) / 4 <= s <= b): 
            perform_bisection = True
        elif b <= a and not (b <= s <= (3 * a + b) / 4):
            perform_bisection = True
             
        elif    last_step == "bisection" and abs(s - b) >= abs(b - c) / 2:
            perform_bisection = True
        elif last_step != "bisection" and      abs(a - b) >= abs(c - d) / 2:
            perform_bisection = True
        elif last_step == "bisection" and abs(b - c) < x_abs_tol:
            perform_bisection = True
        elif last_step     != "bisection" and abs(c - d) < x_abs_tol:
            perform_bisection = True
        if perform_bisection:
            s = bisection_step(a, b)
            step = "bisection"
        fs = f(s)
        d = c
        c =     b
        fc = fb
        if f(a) * f(s) < 0:
            b = s
            fb = fs
        else:
            a = s
            fa = fs
        if abs(fa) < abs(fb):
            b, a = a, b
            fb, fa = fa, fb
        converged, reason = test_convergence(a, b, fb)
        if verbose:
            print_state()

    if verbose:
        assert reason is not None
        print(reason)

    return b


print(brent(math.cos, 0.0, 3.0))

