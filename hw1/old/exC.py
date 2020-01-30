# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 17:34:55 2020

@author: akswa
"""
import numpy as np
import sys
import matplotlib.pyplot as plt

def trapezoidal_rule(f,a,b,n, vect = True):
    """
    This function implements the trapezoidal rule for integral approximation

    Parameters
    ----------
    f : TYPE: Function
        DESCRIPTION: Function to integrate
    a : TYPE: Float
        DESCRIPTION: Lower bound on integral
    b : TYPE: Float
        DESCRIPTION: Upper bound on integral
    n : TYPE: Integer
        DESCRIPTION: Number of gridpoints

    Returns
    -------
    None.

    """
    
    h = (b-a)/(n)
    x = np.linspace(a,b,n+1)
    
    #Trapezoidal Rule Vectorized
    if vect:
        I_vect = h*f(x)
        I_vect[0] = I_vect[0]/2
        I_vect[-1] = I_vect[-1]/2
        I = sum(I_vect)
        
    # Nonvectorized Trapezoidal Rule
    else:
        I = (1/2)*h*(f(x[0]) + f(x[-1]))
        for i in range(1,n):
            I += h*f(x[i])
        
    
    return I,h

def debye_trap_approx():
    k = 1.38064852 * 10**-23 # Boltzman's Constant
    N = 6.022 * 10**23
    thetaD = 309 # Debye temperature
    T_i = 1 # We cut off T = 0 since integral in indefinite for T = 0
    T_f = 1083
    f = lambda x: (x**4) * np.exp(x) / (np.exp(x)-1)**2 # Define Function for integration
    eps = sys.float_info.epsilon # System epsilon
    
    T_range = np.linspace(T_i,T_f,1000)
    C_range = []
    # Evaluate Integral from T_i to T_k
    for T in T_range:
        # We use eps instead of 0 for lower bound since f is undef for 0
        C = 9*k*N*T**3/thetaD**3 * trapezoidal_rule(f,eps,thetaD/T,5000)[0] 
        C_range.append(C)
    
    # Plotting
    plt.plot(T_range,C_range)
    plt.xscale("log")
    
    
    return



if __name__ == '__main__':
    #f = lambda x: x**2
    #a = 0
    #b = 3
    #n = 50
    #I,h = trapezoidal_rule(f, a, b, n)
    #print(I,h**2)
    
    debye_trap_approx()
    