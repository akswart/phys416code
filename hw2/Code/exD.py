# -*- coding: utf-8 -*-
"""
Phys 416
HW 1
Alexander Swart
"""
import numpy as np
import matplotlib.pyplot as plt

def interpf(xi,x,y):
    """
    Function to interpolate between data points
    using Lagrange polynomial (quadratic)
    Inputs
        x    Vector of x coordinates of data points (3 values)
        y    Vector of y coordinates of data points (3 values)
        xi   The x value where interpolation is computed
    Output
        yi   The interpolation polynomial evaluated at xi

    Calculate yi = p(xi) using Lagrange polynomial  
    """
    
    yi = (xi-x[1]) * (xi-x[2]) / ( (x[0]-x[1]) * (x[0]-x[2]) ) * y[0] \
    + (xi-x[0]) * (xi-x[2]) / ( (x[1]-x[0]) * (x[1]-x[2]) ) * y[1] \
    + (xi-x[0]) * (xi-x[1]) / ( (x[2]-x[0]) * (x[2]-x[1]) ) * y[2]
       
    return yi


def interpf_mod(xi,x,y):
    n = len(x)
    yi = 0
    for i in range(n):
        temp_product = 1
        for j in range(n):
            if i != j:
                temp_product *= (xi - x[j]) / (x[i] - x[j])
        temp_product *= y[i]
        yi += temp_product
    return yi


if __name__ == "__main__":
    # part a
    fig1,ax1 = plt.subplots(2,1,sharex=True,sharey=True )
    a = 0
    b = np.pi*2
    x = np.linspace(a,b,100)
    x1 = np.linspace(a,b,5)
    x2 = np.linspace(a,b,7)
    f = lambda x: np.sin(x)
    y1 = f(x1)
    y2 = f(x2)
    xi = np.pi
    
    ax1[0].plot(x,f(x))
    ax1[1].plot(x,f(x))
    ax1[0].scatter(x1,y1) 
    ax1[1].scatter(x2,y2)
    yi1 = []
    yi2 = []
    
    for i in x:
        yi1_t = interpf_mod(i,x1,y1)
        yi2_t = interpf_mod(i,x2,y2)
        yi1.append(yi1_t)
        yi2.append(yi2_t)
    ax1[0].plot(x,yi1)
    ax1[1].plot(x,yi2)
    ax1[0].set_title('5 data points')
    ax1[1].set_title('7 data points')
    fig1.suptitle("Part A: Sine Function")

    # part b
    fig2,ax2 = plt.subplots(2,1,sharex=True,sharey=True )
    a = -2
    b = 2
    x = np.linspace(a,b,100)
    x1 = np.linspace(a,b,5)
    x2 = np.linspace(a,b,7)
    f = lambda x: .5*(1 + np.tanh(10*x))
    y1 = f(x1)
    y2 = f(x2)
    xi = np.pi
    
    ax2[0].plot(x,f(x))
    ax2[1].plot(x,f(x))
    ax2[0].scatter(x1,y1) 
    ax2[1].scatter(x2,y2)
    yi1 = []
    yi2 = []
    
    for i in x:
        yi1_t = interpf_mod(i,x1,y1)
        yi2_t = interpf_mod(i,x2,y2)
        yi1.append(yi1_t)
        yi2.append(yi2_t)
    ax2[0].plot(x,yi1)
    ax2[1].plot(x,yi2)
    ax2[0].set_title('5 data points')
    ax2[1].set_title('7 data points')
    fig2.suptitle("Part B: Heavyside Function")
    

