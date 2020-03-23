#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 13:55:15 2020

@author: aswart
"""
import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3


g = 9.81
k1 = 10
k2 = 20
L1 = .1
L2 = .1
D = .1
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
    D = a[4]
    m = a[5]
    g = a[6]

    
    f = np.zeros(2)
    D = np.matrix(np.zeros([2,2]))
    # Evaluate f(i)
    f[0] = (k1*x[0]*(np.sqrt(x[0]**2 + x[1]**2) - L1))/np.sqrt(x[0]**2 + x[1]**2) + \
           (k2*(x[0]-D)*(np.sqrt((x[0] - D)**2 + x[1]**2) - L2))/np.sqrt((x[0] - D)**2 + x[1]**2)
    
    f[1] = (k1*x[1]*(np.sqrt(x[0]**2 + x[1]**2) - L1))/np.sqrt(x[0]**2 + x[1]**2) + \
           (k2*x[1]*(np.sqrt((x[0] - D)**2 + x[1]**2) - L2))/np.sqrt((x[0] - D)**2 + x[1]**2) - m*g
    
    # Evaluate D(i,j)
    D[0,0] = -a[1];        # df(1)/dx(1)
    D[0,1] = a[0]-x[2];    # df(2)/dx(1)
    D[1,0] = a[1];         # df(1)/dx(2)
    D[1,1] = -1;           # df(2)/dx(2)

    # print('f',f)
    # print('D',D)
    return f,D

if __name__ == "__main__":
    x0,y0 = 1,1
    x = np.array([float(x0),float(y0)])
    xp = np.copy(x)  # Record initial guess for plotting





#"""    
#* Set initial guess and parameters
x0,y0,z0 = input('Enter the initial guess (x,y,z) (row vector): ').split(',')
x = np.array([float(x0),float(y0),float(z0)])  # Copy initial guess
xp = np.copy(x)  # Record initial guess for plotting
# a1,a2,a3 = raw_input('Enter the 3 parameters a: ').split(',')
# a=np.array([float(a1),float(a2),float(a3)])
a1 = float(input('Enter the r: '))
a2 = 10.0
a3 = 8./3.
a=np.array([a1,a2,a3])

#* Loop over desired number of steps 
nStep = 10;   # Number of iterations before stopping
for iStep in range(0,nStep):
	
  #* Evaluate function f and its Jacobian matrix D
    [f, D] = fnewt(x,a);      # fnewt returns value of f and D
  #* Find dx by Gaussian elimination
    dx = np.squeeze(np.asarray(f*np.linalg.inv(D)))
  #* Update the estimate for the root  
    x = x - dx             # Newton iteration for new x
    xp = np.vstack([xp,x]) # Save current estimate for plotting

#* Print the final estimate for the root
    print('After %d iterations the root is'%iStep)
    print(x)

print(xp)
#* Plot the iterations from initial guess to final estimate
#  Clear figure 1 window and bring forward
fig=plt.figure(1)
ax=p3.Axes3D(fig)
# Left plot
ax = fig.add_subplot(1, 2, 1, projection='3d')
ax.plot3D(xp[0:,0],xp[0:,1],xp[0:,2],'o-')#,[x[0],x[1],x[2]],'*')
ax.set_xlabel('x');  ax.set_ylabel('y'); ax.set_zlabel('z');
ax.grid(True)
# Right plot
ax = fig.add_subplot(1, 2, 2, projection='3d')
ax.plot3D(xp[0:,0],xp[0:,1],xp[0:,2],'o-') #,x[0],x[1],x[2],'*');
ax.set_xlabel('x');  ax.set_ylabel('y'); ax.set_zlabel('z');
# view([-127.5, 30]);  # Viewing angle
fig.suptitle('After '+str(nStep)+' iterations, the root is '+str(x))
plt.grid(True)

# # Plot data from lorenz (if available). To write lorenz data, use:
# # >>save xplot; save yplot; save zplot;
# # after running the lorenz program.
flag = int(input('Plot data from lorenz program? (1=Yes/0=No): '))
if( flag == 1 ):
#   figure(2); clf;  # Clear figure 1 window and bring forward
    fig=plt.figure(2)
    ax=p3.Axes3D(fig)
    xplot=np.load('xplot.npy')
    yplot=np.load('yplot.npy')
    zplot=np.load('zplot.npy')
    ax.plot3D(xplot,yplot,zplot,'-')
    ax.plot3D(xp[0:,0],xp[0:,1],xp[0:,2],'o--');
    ax.set_xlabel('x');  ax.set_ylabel('y'); ax.set_zlabel('z');
#   view([40 10]);  # Rotate to get a better view
    plt.grid(True)          # Add a grid to aid perspective

plt.show()
#"""