"""
Contains functions to find the roots of a nonlinear function or a vector
feild using Newton's Method

:AUTHOR: William Skelly
:VERSION: 2.0
"""

import numpy as np
import math

def parabola(x):
    """
    Returns Y value at X for parabola y=x^2

    args:
        x: x coordinate
    returns: the y value at x
    """

    return pow(x,2)

def elliptic_paraboloid(x,y):
    """
    Returns Z value at x,y input for elliptical parabaloid.

    Used as test function for find_zero

    args:
        x: x coordinate
        y: y coordinate
    returns: the z value at (x,y)
    """

    return pow(x,2) + pow(y,2)

def plane(x,y):
    """
    Returns Z value at x,y input for plane.

    Used as test function for find_zero

    args:
        x: x coordinate
        y: y coordinate
    returns: the z value at (x,y)
    """

    return 3*x + 7*y

def calc_jacobian(func1, func2, point):
    """
    Returns a 2x2 matrix containing the jacobian of the vector function
    whose components are func1 and func2, with func1 and func2 being
    functions of 2 variables

    args:
        func1: A single-valued function of 2 variables corresponding to
        the first component of the vector function
        func2: A single-values function of 2 variables corresponding to
        the second component of the vector function
        point: The point to evaluate the jacobian at as a numpy array
    returns: the 2x2 jacobian matrix as a numpy array
    """

    #define step size for x and y
    delta_x = 1e-6
    delta_y = 1e-6

    #calculate partials
    point_val1 = func1(point[0], point[1]) #func1 evalauted at point
    point_val2 = func2(point[0], point[1]) #func2 evaluated at point
    d1_dx = (func1(point[0] + delta_x, point[1]) - point_val1) / delta_x
    d1_dy = (func1(point[0], point[1] + delta_y) - point_val1) / delta_x
    d2_dx = (func2(point[0] + delta_x, point[1]) - point_val2) / delta_x
    d2_dy = (func2(point[0], point[1] + delta_y) - point_val2) / delta_x

    #pack partials into matrix and return
    jacobian = np.array([[d1_dx, d1_dy], [d2_dx, d2_dy]])
    return jacobian

def calc_next_guess2D(func1, func2, guess, debug = False):
    """
    Calculates improved guess for the root of the vector function of 2 variables
    whose first component is func1, and second component func2 via Newton's
    Method.

    args:
        func1: A single-valued function of 2 variables corresponding to
        the first component of the vector function
        func2: A single-values function of 2 variables corresponding to
        the second component of the vector function
        guess: A numpy array repesenting the x and y coordinates of the
        previous guess for root
        debug: Whether to display intermediate results. Default False.
    returns: The improved guess for the root as a numpy array
    """

    #calculate function and Jacobian at guess
    func_vals = np.array([[func1(guess[0], guess[1])], [func2(guess[0], guess[1])]])
    jacobian  = calc_jacobian(func1, func2, guess)

    #solve matrix equation for step size
    step_deltas = np.linalg.solve(jacobian, -1 * func_vals).reshape((-1, ))

    #display intermediate results if in debug mode
    if debug:
        print("Jacobian =\n" + str(jacobian))
        print("Function Vals=\n" + str(func_vals))
        print("\nStep Deltas =\n" + str(step_deltas))
        print("\nGuess =\n" + str(guess))

    #calculate and return improved guess
    next_guess = step_deltas + guess
    return next_guess

def find_zero2D(func1, func2, guess, threshold = 1e-6, debug = False):
    """
    Finds the root for the vector valued function of 2 variables whose x
    component is func1, and y component is func2

    args:
        func1:     A single-valued function of 2 variables corresponding to
        the first component of the vector function
        func2:     A single-values function of 2 variables corresponding to
        the second component of the vector function
        guess:     A numpy array repesenting the x and y coordinates of the
        initial guess for root
        threshold: the maximum value of f(x,y) to be considered as solved
        debug:     whether to print intermediate output. Default False.
    returns: The improved guess for the root as a numpy array
    """

    #setup loop
    converged  = False
    iter_max   = False
    iteration  = 0
    iter_limit = 100

    #loop until converged or iteration limit reached
    while not converged and not iter_max:
        #display iteration and guess for debug
        if debug:
            print("\nIteration: " + str(iteration) + "\t\tGuess: " + str(guess))
        
        #calculate next guess and error
        guess = calc_next_guess2D(func1, func2, guess)
        error = math.sqrt(pow(func1(guess[0], guess[1]), 2) +
                          pow(func2(guess[0], guess[1]), 2))
        #if debug:
        #    print("\tError: " + str(error))

        #check convergance criteria and iteration count
        iteration += 1
        if error < threshold:
            converged = True
        if iteration > iter_limit:
            iter_max = True

    #return result
    if converged:
        return guess
    else:
        message = "Iteration limit reached. Final guess: " + str(guess)
        return message

def calc_slope(func, x_val):
    """
    Calculates the derivative of the function at a single point

    args:
        func: the function to calculate the derivative of
        x_val: the value of the independant variable to calculate the derivative
        at
    returns: the slope of func at x_val
    """

    x_delta  = 1e-6
    func_val = func(x_val)
    slope    = (func(x_val + x_delta) - func_val) / x_delta
    return slope

def calc_next_guess1D(func, guess, debug = False):
    """
    Calculates improved guess for the root of a single variable function via
    Newton's Method.

    args:
        func:      the function to find the root of
        guess:     initial guess for the x location of the root
        debug:     whether to print intermediate output. Default False.
    returns: The improved guess for the root
    """

    #calculate function and derivative at guess
    func_val = func(guess)
    slope    = calc_slope(func, guess)

    #display for debug
    if debug:
        print("Guess: " + str(guess) + "\tValue: " + str(func_val) +
              "\tSlope: " + str(slope))

    #calculate and return improved guess
    next_guess = guess - (func_val / slope)
    return next_guess

def find_zero1D(func, guess, threshold = 1e-6, debug = False):
    """
    Finds the root for a single variable function

    args:
        func:      the function to find the root of
        guess:     initial guess for the x location of the root
        threshold: the maximum value of f(x) to be considered as solved
        debug:     whether to print intermediate output. Default False.
    returns: The improved guess for the root
    """

    #setup loop
    converged  = False
    iter_max   = False
    iteration  = 0
    iter_limit = 100

    #loop until converged or iteration limit reached
    while not converged and not iter_max:
        #display iteration and guess for debug
        if debug:
            print("\nIteration: " + str(iteration) + "\t\tGuess: " + str(guess))
        
        #calculate next guess and error
        guess = calc_next_guess1D(func, guess)
        error = func(guess)

        #check convergance criteria and iteration count
        iteration += 1
        if error < threshold:
            converged = True
        if iteration > iter_limit:
            iter_max = True

    #return result
    if converged:
        return guess
    else:
        message = "Iteration limit reached. Final guess: " + str(guess)
        return message
        
def main():
    point = np.array([10, 12])
    print(find_zero2D(elliptic_paraboloid, plane, point, debug = True))
    #print(find_zero1D(parabola, 1, debug = True))

if __name__ == "__main__":
    main()
