#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 13:55:15 2020

@author: aswart
"""
import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import warnings

g = 9.81
k1 = 10
k2 = 20
L1 = .1
L2 = .1
d = .1
m = .1

V = lambda x,y :  .5*k1*(np.sqrt(x**2+y**2)-L1)**2 + .5*k2*(np.sqrt((x - D)**2+y**2)-L1)**2 - m*g*y
# newtn - Program to solve a system of nonlinear equations
# using Newtons method.  Equations defined by function fnewt.
# clear all;  help newtn;  # Clear memory and print header

def  fnewt(x,a):
    """
      Function used by the N-variable Newtons method
      Inputs
        x     State vector [x y]
        a     Parameters [k1,k2,L1,L2,D,m,g]
      Outputs
        f     Lorenz model r.h.s. [dx/dt dy/dt]
        D     Jacobian matrix, D(i,j) = df(j)/dx(i)
    """
    k1 = a[0]
    k2 = a[1]
    L1 = a[2]
    L2 = a[3]
    d = a[4]
    m = a[5]
    g = a[6]

    
    f = np.zeros(2)
    D = np.matrix(np.zeros([2,2]))
    # Evaluate f(i)
    # force in x dir
    f[0] = -( (k1*x[0]*(np.sqrt(x[0]**2 + x[1]**2) - L1))/np.sqrt(x[0]**2 + x[1]**2) + \
           (k2*(x[0]-d)*(np.sqrt((x[0] - d)**2 + x[1]**2) - L2))/np.sqrt((x[0] - d)**2 + x[1]**2) )
    # force in y dir
    f[1] = -( (k1*x[1]*(np.sqrt(x[0]**2 + x[1]**2) - L1))/np.sqrt(x[0]**2 + x[1]**2) + \
           (k2*x[1]*(np.sqrt((x[0] - d)**2 + x[1]**2) - L2))/np.sqrt((x[0] - d)**2 + x[1]**2) - m*g)
    
    
    # Evaluate D(i,j) (Thanks to wolfram)
    D[0,0] = -( -(k1*L1*x[1]**2)/(x[0]**2+x[1]**2)**(3/2) - (k2*L2*x[1]**2)/((d-x[0])**2+x[1]**2)**(3/2) + (k1+k2) )      # df(1)/dx(1)
    
    D[0,1] = -( x[1]*((k1*L1*x[0])/(x[0]**2 + x[1]**2)**(3/2) - (k2*L2*(d - x[0]))/(d**2 - 2*d*x[0] + x[0]**2 + x[1]**2)**(3/2)) )   # df(2)/dx(1)
    
    D[1,0] = -( (k1*L1*x[0]*x[1])/(x[0]**2 + x[1]**2)**(3/2) + k2*x[1]*(d - x[0])*((np.sqrt(d**2 - 2*d*x[0] + x[0]**2 + x[1]**2) - L2) \
               /(d**2 - 2*d*x[0] + x[0]**2 + x[1]**2)**(3/2) - 1/(d**2 - 2*d*x[0] + x[0]**2 + x[1]**2)) )        # df(1)/dx(2)
        
    D[1,1] = -( k1*(1-(L1*x[0]**2)/(x[0]**2+x[1]**2)**(3/2)) + k2*(((d**2 - 2*d*x[0] +x[0]**2) * \
           (np.sqrt(d**2 - 2*d*x[0] + x[0]**2 + x[1]**2) - L2))/(d**2 - 2*d*x[0] + x[0]**2 + x[1]**2)**(3/2) + \
                                                  x[1]**2/(d**2-2*d*x[0] + x[0]**2 + x[1]**2)) )        # df(2)/dx(2)
        
    # Test if D is singular
    det_D = np.linalg.det(D)
    if abs(det_D) <= np.finfo(np.float32).eps:
        warnings.warn("Jacobian Matrix is singular")
    
    # print('f',f)
    # print('D',D)
    return f,D

if __name__ == "__main__":
    # Set initial guess and parameters
    x0,y0 = 1,1
    x = np.array([float(x0),float(y0)])
    xp = np.copy(x)  # Record initial guess for plotting

    # Set params
    a=np.array([k1,k2,L1,L2,d,m,g])
    
    # Loop over desired number of steps 
    nStep = 10   # Number of iterations before stopping
    for iStep in range(0,nStep):
    	
      # Evaluate function f and its Jacobian matrix D
        [f, D] = fnewt(x,a)      # fnewt returns value of f and D
      # Find dx by Gaussian elimination
        dx = np.squeeze(np.asarray(f*np.linalg.inv(D)))
      # Update the estimate for the root  
        x = x - dx             # Newton iteration for new x
        xp = np.vstack([xp,x]) # Save current estimate for plotting
    
    #* Print the final estimate for the root
        print('After %d iterations the root is'%iStep)
        print(x)
    
    print(xp)
    
    fin_posx,fin_posy = xp[-1,0],xp[-1,1]
    print("The final position is x: %f, y: %f" % (fin_posx,fin_posy))